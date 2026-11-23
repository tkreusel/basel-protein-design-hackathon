"""LigandMPNN re-design of stage-1 success backbones, prepared for refolding.

Produces, in one pass:
  <work>/ligandmpnn/<design_id>/   raw LigandMPNN output
  <work>/armB_sequences.fasta      every kept sequence
  <work>/boltz_in/*.yaml           one Boltz input per kept sequence
  <work>/designs.json              catalytic positions, for scoring the refold

Rules enforced (plan 6.3):
  * only designs that PASSED the stage-1 gate are redesigned -- consult the
    gate, not just the crash flag
  * catalytic residues are fixed, and the output is VERIFIED; any sequence
    that altered one is dropped with the reason recorded
  * 6NT is supplied as ligand context

Usage: python -I run_ligandmpnn.py <root> <stage1 run> <work dir> [n_seq]
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp.record import read_manifest  # noqa: E402

GATE = "stage1_motif"
EXPECTED = {"base": "D", "donor": "Q", "wall": "W"}
ENV_ROOT = os.environ.get("ENV_ROOT", "/mnt/ref/software/envs")
MPNN_ROOT = os.environ.get("LIGANDMPNN_ROOT", "/mnt/ref/software/LigandMPNN")
CKPT = os.environ.get("LIGANDMPNN_CHECKPOINT",
                      "/mnt/ref/data_weights/ligandmpnn/ligandmpnn_v_32_010_25.pt")


def read_fasta(p: Path):
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


def main(root: Path, run: Path, work: Path, n_seq: int = 4) -> None:
    recs = list(read_manifest(run / "manifest.jsonl"))
    ok = [r for r in recs if r.status == "ok" and r.gates.get(GATE) == "pass"]
    print(f"stage-1: {len(recs)} designs, {len(ok)} passed the gate -> redesigning these")

    mpnn_dir = work / "ligandmpnn"
    boltz_in = work / "boltz_in"
    for d in (mpnn_dir, boltz_in):
        d.mkdir(parents=True, exist_ok=True)

    fasta_lines, designs = [], {}
    n_kept = n_dropped = n_failed = 0

    for rec in ok:
        did = rec.design_id
        pdb = Path(rec.paths.get("backbone_pdb", ""))
        if not pdb.exists():
            print(f"  [FAIL] {did}: backbone pdb missing")
            n_failed += 1
            continue
        cat = rec.catalytic_residues
        wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat["wall"]
        fixed = f"A{cat['base']} A{cat['donor']} A{wall}"
        dest = mpnn_dir / did

        cmd = [f"{ENV_ROOT}/env-ligandmpnn/bin/python", f"{MPNN_ROOT}/run.py",
               "--model_type", "ligand_mpnn",
               "--checkpoint_ligand_mpnn", CKPT,
               "--pdb_path", str(pdb), "--out_folder", str(dest),
               "--fixed_residues", fixed,
               "--ligand_mpnn_use_atom_context", "1",
               "--ligand_mpnn_use_side_chain_context", "1",
               "--batch_size", "1", "--number_of_batches", str(n_seq),
               "--pack_side_chains", "0", "--seed", "1", "--verbose", "0"]
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            print(f"  [FAIL] {did}: LigandMPNN rc={p.returncode} {p.stderr[-200:]}")
            n_failed += 1
            continue

        fas = sorted((dest / "seqs").glob("*.fa")) if (dest / "seqs").exists() else []
        if not fas:
            print(f"  [FAIL] {did}: no .fa output")
            n_failed += 1
            continue

        entries = read_fasta(fas[0])
        kept_here = 0
        # entry 0 is the input/native row; designs follow.
        for idx, (hdr, seq) in enumerate(entries[1:], start=1):
            bad = []
            for role, pos in (("base", cat["base"]), ("donor", cat["donor"]),
                              ("wall", wall)):
                got = seq[pos - 1] if pos - 1 < len(seq) else "?"
                if got != EXPECTED[role]:
                    bad.append(f"{role}@{pos}={got}!={EXPECTED[role]}")
            if bad:
                n_dropped += 1
                print(f"  [drop] {did}_s{idx}: catalytic changed: {';'.join(bad)}")
                continue
            name = f"{did}_s{idx}"
            fasta_lines.append(f">{name} len={len(seq)}\n{seq}")
            (boltz_in / f"{name}.yaml").write_text(yaml.safe_dump({
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": seq, "msa": "empty"}},
                    {"ligand": {"id": "B", "ccd": "6NT"}},
                ]}, sort_keys=False))
            designs[name] = {
                "design_id": name, "parent_design": did,
                "stage1_rank": len(designs) + 1, "length": len(seq),
                "catalytic": {"base": cat["base"], "donor": cat["donor"],
                              "wall": wall},
                "stage1_scores": {k: v for k, v in rec.scores.items()
                                  if k.startswith("geom.")},
                "source_run": run.name, "arm": "armB_ligandmpnn",
            }
            kept_here += 1
            n_kept += 1
        print(f"  [ok]   {did}: kept {kept_here}/{max(len(entries)-1,0)}")

    (work / "armB_sequences.fasta").write_text("\n".join(fasta_lines) + "\n")
    (work / "designs.json").write_text(json.dumps(designs, indent=2))
    print(f"\nkept={n_kept}  dropped_catalytic_changed={n_dropped}  failed={n_failed}")
    print(f"boltz inputs: {len(designs)} in {boltz_in}")
    if not designs:
        raise SystemExit("no sequences survived -- nothing to refold")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]),
         int(sys.argv[4]) if len(sys.argv) > 4 else 4)
