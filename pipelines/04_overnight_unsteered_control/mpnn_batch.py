"""Stage 3: LigandMPNN for every selected design, mirroring Maxi's prep_mpnn.py.

Maxi's settings (prep_mpnn.py), kept as is:
  model ligand_mpnn v_32_010_25, temperature 0.1, 4 sequences, seed 111,
  --ligand_mpnn_use_atom_context 1 (side-chain context left at its default),
  motif residues FIXED, and the residue types K R H D E N Q S T Y C FORBIDDEN at
  every non-motif residue within 6 A of the base carboxylate O (OD1/OD2).

One difference, on purpose: all designs go through ONE LigandMPNN process using
--pdb_path_multi / --fixed_residues_multi / --omit_AA_per_residue_multi, so the
checkpoint is loaded once instead of once per design (our per-design runs spent
most of their time reloading it).

Verification (we do not trust the flags): after the run every sequence is
checked for (a) motif residues unchanged and (b) no forbidden residue type at a
forbidden position. Violations are listed and the sequence is dropped.

Writes: <out>/mpnn/in/<design>.pdb, fixed/omit/pdb json, <out>/mpnn/seqs/*.fa,
        <out>/sequences.csv  (design, seq_id, source, sequence, fixed, n_omit)

Usage: python -I mpnn_batch.py <selection.csv> <out dir> [--prepare-only]
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
from biotite.structure import concatenate
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
ENV_ROOT = os.environ.get("ENV_ROOT", "/mnt/ref/software/envs")
MPNN = os.environ.get("LIGANDMPNN_ROOT", "/mnt/ref/software/LigandMPNN")
CKPT = os.environ.get("LIGANDMPNN_CHECKPOINT",
                      "/mnt/ref/data_weights/ligandmpnn/ligandmpnn_v_32_010_25.pt")
OMIT = "KRHDENQSTYC"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def read_fa(p: Path):
    out, name, buf = [], None, []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def main(sel_csv: Path, out: Path, prepare_only: bool = False) -> None:
    rows = list(csv.DictReader(open(sel_csv)))
    base = out / "mpnn"
    for d in ("in", "seqs"):
        (base / d).mkdir(parents=True, exist_ok=True)

    pdb_multi, fixed_multi, omit_multi, meta, seqs_rfd4 = {}, {}, {}, {}, {}
    for r in rows:
        name = r["design"]
        a = get_structure(CIFFile.read(r["cif"]), model=1)
        a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]
        lig, prot = a[a.hetero], a[~a.hetero]
        prot.chain_id[:] = "A"; lig.chain_id[:] = "B"; lig.res_id[:] = 1
        pdb = base / "in" / f"{name}.pdb"
        f = PDBFile(); f.set_structure(concatenate([prot, lig])); f.write(pdb)

        fixed_ids = [int(float(r[k])) for k in ("base_resid", "partner_resid", "wall_resid")]
        bo = prot.coord[(prot.res_id == fixed_ids[0]) & np.isin(prot.atom_name, ("OD1", "OD2"))]
        tree = cKDTree(prot.coord)
        near = {int(prot.res_id[i]) for i in
                np.unique(np.concatenate(tree.query_ball_point(bo, 6.0)))} - set(fixed_ids)
        key = str(pdb)
        pdb_multi[key] = ""
        fixed_multi[key] = " ".join(f"A{i}" for i in fixed_ids)
        omit_multi[key] = {f"A{i}": OMIT for i in sorted(near)}
        ca = prot[prot.atom_name == "CA"]
        seqs_rfd4[name] = "".join(AA3.get(n, "X") for n in ca.res_name)
        meta[name] = dict(fixed=fixed_multi[key], omit=sorted(near), ids=fixed_ids,
                          res_ids=[int(x) for x in ca.res_id])

    for fn, obj in (("pdb_multi.json", pdb_multi), ("fixed_multi.json", fixed_multi),
                    ("omit_multi.json", omit_multi)):
        (base / fn).write_text(json.dumps(obj))
    print(f"prepared {len(rows)} designs; mean forbidden-zone size "
          f"{np.mean([len(m['omit']) for m in meta.values()]):.1f} residues")
    if prepare_only:
        return

    cmd = [f"{ENV_ROOT}/env-ligandmpnn/bin/python", f"{MPNN}/run.py",
           "--model_type", "ligand_mpnn", "--checkpoint_ligand_mpnn", CKPT,
           "--pdb_path_multi", str(base / "pdb_multi.json"),
           "--fixed_residues_multi", str(base / "fixed_multi.json"),
           "--omit_AA_per_residue_multi", str(base / "omit_multi.json"),
           "--out_folder", str(base), "--number_of_batches", "4", "--batch_size", "1",
           "--temperature", "0.1", "--seed", "111",
           "--ligand_mpnn_use_atom_context", "1", "--save_stats", "0"]
    print("running one LigandMPNN process for all designs ...", flush=True)
    p = subprocess.run(cmd, capture_output=True, text=True)
    (base / "ligandmpnn.log").write_text(p.stdout + "\n--- stderr ---\n" + p.stderr)
    if p.returncode != 0:
        raise SystemExit(f"LigandMPNN failed rc={p.returncode}: {p.stderr[-500:]}")

    out_rows, bad = [], 0
    for name, m in meta.items():
        out_rows.append(dict(design=name, seq_id=f"{name}__rfd4", source="rfd4",
                             sequence=seqs_rfd4[name], fixed=m["fixed"], n_omit=len(m["omit"])))
        fa = base / "seqs" / f"{name}.fa"
        if not fa.exists():
            print(f"  [FAIL] no output for {name}"); bad += 1; continue
        recs = read_fa(fa)
        native = recs[0][1].split(":")[0]
        idx = {rid: i for i, rid in enumerate(m["res_ids"])}
        for k, (_, s) in enumerate(recs[1:], 1):
            s = s.split(":")[0]
            viol = []
            for rid in m["ids"]:                      # motif residues unchanged
                if s[idx[rid]] != native[idx[rid]]:
                    viol.append(f"motif {rid}")
            for rid in m["omit"]:                     # forbidden types absent
                if s[idx[rid]] in OMIT:
                    viol.append(f"forbidden {s[idx[rid]]}@{rid}")
            if viol:
                bad += 1; print(f"  [drop] {name}_mpnn{k}: {viol[:3]}"); continue
            out_rows.append(dict(design=name, seq_id=f"{name}__mpnn{k}", source=f"mpnn{k}",
                                 sequence=s, fixed=m["fixed"], n_omit=len(m["omit"])))
    with open(out / "sequences.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0])); w.writeheader(); w.writerows(out_rows)
    n_m = sum(r["source"].startswith("mpnn") for r in out_rows)
    print(f"sequences: {len(out_rows)} ({n_m} MPNN + {len(out_rows) - n_m} RFD4); "
          f"dropped/failed {bad}  ->  {out / 'sequences.csv'}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), "--prepare-only" in sys.argv)
