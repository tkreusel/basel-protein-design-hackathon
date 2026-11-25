"""Stage 3: LigandMPNN sequences for the gated designs.

usage: mpnn.py --sel sel.csv --out DIR [--shard k/n] [--limit N]
       mpnn.py --merge --out DIR         (combine the shard tables)

sel.csv columns: design, cif, motif, arm. Per design (motif residues fixed in every variant, located by superposing the
reference theozyme on the design's ligand):
  rfd4        the RFD4-native sequence (kept as a control)
  mpnn1..4    full redesign, ligand context, amino acids in `omit_aa` forbidden at non-motif residues within 6 A of the base carboxylate O
  shell1      steered arms only: motif + catalytic shell fixed at the RFD4 sequence, the rest redesigned. Keeps the environment the
              reward-guided search built (TMol, nitro donor, polar exclusion).
Writes DIR/mpnn/<design>/{in.pdb,omit.json,full/,shell/} and DIR/sequences[_shardKofN].csv with
design, seq_id, source, variant, sequence, fixed, n_omit, n_shell_fixed.
Environment: LIGANDMPNN_DIR (clone), LIGANDMPNN_PY (interpreter), LIGANDMPNN_WEIGHTS (directory with the checkpoint).
"""
import argparse
import glob
import json
import os
import subprocess

import numpy as np
import pandas as pd
import biotite.structure as bs
import biotite.structure.io.pdb as pdb
from scipy.spatial import cKDTree

from . import config, structure

COLUMNS = ["design", "seq_id", "source", "variant", "sequence", "fixed", "n_omit", "n_shell_fixed"]


def split_design(path):
    """Protein (chain A) and ligand (chain B, residue 1) of an RFD4 design CIF."""
    atoms = structure.load_cif(path)
    ligand, protein = atoms[atoms.hetero], atoms[~atoms.hetero]
    protein.chain_id[:] = "A"
    ligand.chain_id[:] = "B"
    ligand.res_id[:] = 1
    return protein, ligand


def design_context(protein, ligand, motif, is_steered, cfg=None, reference=None):
    """Residue sets for one design.

    Returns dict(motif_ids, omit_ids, shell_ids): motif residues; residues with omit_aa forbidden (within base_shell_A of the base O,
    excluding the motif); residues fixed in the shell variant (None for unsteered arms).
    Raises ValueError if a motif residue cannot be located.
    """
    cfg = cfg or config.load("sequence_design")
    located = structure.locate_motif(protein, ligand, motif, reference)
    motif_ids = [rid for _, rid, _ in located]
    if any(r is None for r in motif_ids):
        raise ValueError(f"motif not located: {located}")
    base_id = next(rid for res, rid, _ in located if res == "ASP")
    base_o = protein[(protein.res_id == base_id) & np.isin(protein.atom_name, config.load("motifs")["base_oxygens"])]
    near_base = {int(protein.res_id[i]) for i in np.unique(np.concatenate(cKDTree(protein.coord).query_ball_point(base_o.coord, cfg["base_shell_A"])))}
    omit = sorted(near_base - set(motif_ids))
    shell = None
    if is_steered:
        d_lig = cKDTree(ligand.coord).query(protein.coord)[0]
        near_lig = {int(r) for r, d in zip(protein.res_id, d_lig) if d <= cfg["shell_variant"]["ligand_shell_A"]}
        shell = sorted(near_lig | near_base | set(motif_ids))
    return dict(motif_ids=motif_ids, omit_ids=omit, shell_ids=shell)


def run_ligandmpnn(out, tag, pdb_path, fixed_ids, n_batches, seed, omit_json=None, cfg=None):
    cfg = cfg or config.load("sequence_design")
    lm_dir = config.env_path("LIGANDMPNN_DIR")
    cmd = [str(config.env_path("LIGANDMPNN_PY")), str(lm_dir / "run.py"), "--model_type", "ligand_mpnn",
           "--checkpoint_ligand_mpnn", str(config.env_path("LIGANDMPNN_WEIGHTS") / cfg["checkpoint"]),
           "--pdb_path", pdb_path, "--out_folder", f"{out}/{tag}", "--number_of_batches", str(n_batches), "--batch_size", "1",
           "--temperature", str(cfg["temperature"]), "--seed", str(seed), "--fixed_residues", " ".join(f"A{i}" for i in fixed_ids),
           "--ligand_mpnn_use_atom_context", "1", "--save_stats", "0"]
    if omit_json:
        cmd += ["--omit_AA_per_residue", omit_json]
    if not glob.glob(f"{out}/{tag}/seqs/*.fa"):                 # resumable
        result = subprocess.run(cmd, cwd=lm_dir, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        if result.returncode:
            raise RuntimeError("LigandMPNN failed: " + result.stderr[-1500:])
    records = open(glob.glob(f"{out}/{tag}/seqs/*.fa")[0]).read().split(">")[1:]
    return [rec.split("\n", 1)[1].replace("\n", "").split(":")[0] for rec in records[1:]]       # record 0 is the input sequence


def sequences_for_design(row, out_root, cfg=None):
    """All sequence rows of one design (runs LigandMPNN)."""
    cfg = cfg or config.load("sequence_design")
    out = f"{out_root}/mpnn/{row.design}"
    os.makedirs(out, exist_ok=True)
    protein, ligand = split_design(row.cif)
    f = pdb.PDBFile()
    f.set_structure(bs.concatenate([protein, ligand]))
    f.write(f"{out}/in.pdb")
    ctx = design_context(protein, ligand, row.motif, config.is_steered(row.arm), cfg)
    json.dump({f"A{i}": cfg["omit_aa"] for i in ctx["omit_ids"]}, open(f"{out}/omit.json", "w"))
    fixed = " ".join(f"A{i}" for i in ctx["motif_ids"])
    rows = [dict(design=row.design, seq_id=f"{row.design}__rfd4", source="rfd4", variant="native", sequence=structure.sequence(protein),
                 fixed=fixed, n_omit=len(ctx["omit_ids"]), n_shell_fixed=0)]
    for j, seq in enumerate(run_ligandmpnn(out, "full", f"{out}/in.pdb", ctx["motif_ids"], cfg["n_sequences"], cfg["seed_full"], f"{out}/omit.json", cfg), 1):
        rows.append(dict(design=row.design, seq_id=f"{row.design}__mpnn{j}", source=f"mpnn{j}", variant="full", sequence=seq, fixed=fixed,
                         n_omit=len(ctx["omit_ids"]), n_shell_fixed=0))
    if ctx["shell_ids"] is not None:
        sv = cfg["shell_variant"]
        for j, seq in enumerate(run_ligandmpnn(out, "shell", f"{out}/in.pdb", ctx["shell_ids"], sv["n_sequences"], sv["seed"], cfg=cfg), 1):
            rows.append(dict(design=row.design, seq_id=f"{row.design}__shell{j}", source=f"shell{j}", variant="shell", sequence=seq, fixed=fixed,
                             n_omit=0, n_shell_fixed=len(ctx["shell_ids"])))
    return rows


def merge_shards(out):
    """Concatenate sequences_shard*.csv into sequences.csv (shards without rows are skipped)."""
    parts = [pd.read_csv(f) for f in sorted(glob.glob(f"{out}/sequences_shard*.csv")) if len(open(f).read().strip().splitlines()) > 1]
    merged = pd.concat(parts, ignore_index=True)
    merged.to_csv(f"{out}/sequences.csv", index=False)
    print("sequences:", len(merged), merged.groupby("variant").size().to_dict())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sel")
    ap.add_argument("--out", required=True)
    ap.add_argument("--shard", help="k/n: process every n-th design starting at k")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--merge", action="store_true", help="merge the shard tables in --out into sequences.csv and exit")
    a = ap.parse_args()
    if a.merge:
        return merge_shards(a.out)
    if not a.sel:
        ap.error("--sel is required unless --merge")
    sel = pd.read_csv(a.sel)
    if a.shard:
        k, n = map(int, a.shard.split("/"))
        sel = sel.iloc[k::n]
    if a.limit:
        sel = sel.head(a.limit)
    rows = []
    for row in sel.itertuples():
        try:
            rows += sequences_for_design(row, a.out)
            print(row.design, "ok", flush=True)
        except ValueError as exc:
            print("MOTIF NOT LOCATED", row.design, exc, flush=True)
    name = f"sequences_shard{a.shard.replace('/', 'of')}.csv" if a.shard else "sequences.csv"
    pd.DataFrame(rows, columns=COLUMNS).to_csv(f"{a.out}/{name}", index=False)
    print("wrote", f"{a.out}/{name}", len(rows))


if __name__ == "__main__":
    main()
