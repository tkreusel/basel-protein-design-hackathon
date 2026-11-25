"""Select the extended evaluation set of a finished run.

usage: extended_set.py RUN [--k 20]

Sequences that did not pass the strict Boltz gate but fold (pLDDT and CA RMSD pass) are scored on how many relaxed criteria they meet
(ligand RMSD <= 2.5 A, motif tips <= 1.5 A, refold geometry, ligand ipTM). At most `max_per_design` per design, `k` in total, best first.
Writes RUN/survivors_x.txt and RUN/extended_selection.csv.
"""
import argparse

import pandas as pd

from . import config


def select_extended(boltz, strict_ids, k=None, cfg=None):
    cfg = cfg or config.load("gates")["extended_set"]
    k = cfg["max_sequences"] if k is None else k
    rel = cfg["relaxed"]
    x = boltz[~boltz.seq_id.isin(strict_ids) & boltz.ok_plddt & boltz.ok_ca].copy()
    x["relax"] = ((x.ligand_rmsd_ca_aligned <= rel["ligand_rmsd_max"]).astype(int) + (x.motif_tip_rmsd_ca_aligned <= rel["motif_tip_rmsd_max"]).astype(int)
                  + x.geom_ok_refold.astype(int) + x.ok_liptm.astype(int))
    x = x.sort_values(["relax", "ligand_rmsd_ca_aligned"], ascending=[False, True])
    x["rk"] = x.groupby("design").cumcount()
    return x[x.rk < cfg["max_per_design"]].head(k)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run")
    ap.add_argument("--k", type=int)
    a = ap.parse_args()
    boltz = pd.read_csv(f"{a.run}/boltz1_eval.csv")
    strict = set(open(f"{a.run}/survivors.txt").read().split())
    x = select_extended(boltz, strict, a.k)
    open(f"{a.run}/survivors_x.txt", "w").write("\n".join(x.seq_id) + "\n")
    x.to_csv(f"{a.run}/extended_selection.csv", index=False)
    print("extended set:", len(x), "sequences from", x.design.nunique(), "designs")


if __name__ == "__main__":
    main()
