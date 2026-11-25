"""Write Boltz-2 inputs (no MSA; protein chain A, ligand chain B as CCD code) for a table of sequences.

usage: boltz_prep.py --seq sequences.csv --ligand 6NT|H5J|6VP|apo --out DIR [--ids ids.txt] [--chunk 40]

File name: <seq_id>__<LIGAND>.yaml, so the evaluation suite resolves entry = seq_id and arm = ligand.
Chunk directories DIR/chunk_NNN/ allow array tasks; the campaign uses one chunk per run.
"""
import argparse
import os

import pandas as pd
import yaml


def write_inputs(df, ligand, out, chunk=40):
    n = 0
    for k in range(0, len(df), chunk):
        chunk_dir = f"{out}/chunk_{k // chunk:03d}"
        os.makedirs(chunk_dir, exist_ok=True)
        for r in df.iloc[k:k + chunk].itertuples():
            seqs = [{"protein": {"id": "A", "sequence": r.sequence, "msa": "empty"}}]
            if ligand != "apo":
                seqs.append({"ligand": {"id": "B", "ccd": ligand}})
            with open(f"{chunk_dir}/{r.seq_id}__{ligand}.yaml", "w") as fh:
                yaml.safe_dump({"version": 1, "sequences": seqs}, fh, sort_keys=False)
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seq", required=True)
    ap.add_argument("--ligand", required=True, choices=["6NT", "H5J", "6VP", "apo"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--ids")
    ap.add_argument("--chunk", type=int, default=40)
    a = ap.parse_args()
    df = pd.read_csv(a.seq)
    df = df[df.sequence.str.fullmatch("[ACDEFGHIKLMNPQRSTVWY]+")]
    if a.ids:
        df = df[df.seq_id.isin(set(open(a.ids).read().split()))]
    n = write_inputs(df, a.ligand, a.out, a.chunk)
    print(n, "inputs ->", a.out, f"({-(-len(df) // a.chunk)} chunks)")


if __name__ == "__main__":
    main()
