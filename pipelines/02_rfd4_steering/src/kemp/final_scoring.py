"""Final filter and composite score.

Filter: ipSAE2(TSA, mean of the Boltz samples) above `final_selection.ipsae_tsa_min` (0.8). Score: equal-weight composite of seven metrics,
all oriented so that lower is better (their correlations with activity in the literature set are negative), each z-scored over the filtered
population; the mean percentile rank is kept as a robustness check.

  S01  tsa_pk.int.lk_ball_iso - prod.int.lk_ball_iso
  S03  sub.int.fa_elec - prod.int.fa_elec
  S04  tsa_pk.int.total
  S05  tsa_pk.E_total - prod.E_total
  S06  tsa_pk.buns.n_unsat_ligand - prod.buns.n_unsat_ligand
  S07  polar-residue ligand interaction (S,T,N,Q,Y,H,K,R,D,E; summed tmol per-residue energies), TSA minus substrate
  S12  mean Boltz PAE from catalytic residues to pocket residues (within 8 A of the ligand, catalytic excluded), substrate minus product

usage: final_scoring.py [--camp NAME]   reads $KEMP_RESULTS/tables/evaluated_sequences_full_suite.csv,
                                        writes final_scoring.csv and evaluated_with_final_metrics.csv next to it
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

from . import config

POLAR = set("STNQYHKRDE")
SCORE_COLUMNS = ["S01_lkballiso_TSA_minus_PROD", "S03_faelec_SUB_minus_PROD", "S04_interaction_TSA", "S05_Etotal_TSA_minus_PROD",
                 "S06_unsat_ligand_polars_TSA_minus_PROD", "S07_tmol_polar_residue_int_TSA_minus_SUB", "S12_boltz_cat_to_pocket_PAE_SUB_minus_PROD"]


def add_energy_terms(e):
    """The five metrics that are differences or copies of columns already in the evaluation table."""
    e = e.copy()
    e["S01_lkballiso_TSA_minus_PROD"] = e["tsa_pk.int.lk_ball_iso"] - e["prod.int.lk_ball_iso"]
    e["S03_faelec_SUB_minus_PROD"] = e["sub.int.fa_elec"] - e["prod.int.fa_elec"]
    e["S04_interaction_TSA"] = e["tsa_pk.int.total"]
    e["S05_Etotal_TSA_minus_PROD"] = e["tsa_pk.E_total"] - e["prod.E_total"]
    e["S06_unsat_ligand_polars_TSA_minus_PROD"] = e["tsa_pk.buns.n_unsat_ligand"] - e["prod.buns.n_unsat_ligand"]
    return e


def per_run_terms(x, root, seqs):
    """S07 (tmol per-residue energies) and S12 (Boltz PAE) need the run directory."""
    run = f"{root}/{x.motif}/{x.arm}/{x.batch}"
    ev_dir, forest = (f"{run}/eval_x", f"{run}/eval_forest_x") if x.eval_set == "extended" else (f"{run}/eval", f"{run}/eval_forest")
    out = {"S07": np.nan, "S12": np.nan}
    jf = f"{ev_dir}/work/tmol/{x.seq_id}_boltz_tsa_boltz_m0.json"
    if not os.path.exists(jf):
        return out
    detail, sq = json.load(open(jf))["detail"], seqs.get(x.seq_id, "")

    def polar(tag):
        per_res = detail.get(tag, {}).get("int_per_res")
        if not per_res or not sq:
            return np.nan
        return float(sum(v for k, v in per_res.items() if 0 < int(k) <= len(sq) and sq[int(k) - 1] in POLAR))
    out["S07"] = polar("tsa_pk") - polar("sub")

    site = detail.get("site", {})
    pocket = set(site.get("pocket_resids_8A") or [])
    table = pd.read_csv(f"{run}/mpnn/sequences.csv")
    fixed = [int(t[1:]) for t in str(table.loc[table.seq_id == x.seq_id, "fixed"].iloc[0]).split()]
    cat, pk = [i - 1 for i in fixed], [i - 1 for i in pocket if i not in fixed]
    if cat and pk:
        pae = {}
        for lig in ("H5J", "6VP"):
            vals = []
            for f in glob.glob(f"{forest}/{x.seq_id}__{lig}/pae_*_model_*.npz"):
                z = np.load(f)
                vals.append(float(z[z.files[0]][np.ix_(cat, pk)].mean()))
            pae[lig] = np.mean(vals) if vals else np.nan
        out["S12"] = pae["H5J"] - pae["6VP"]
    return out


def composite(e, ipsae_min, min_metrics):
    """Filter on ipSAE2 and add z-score and rank composites. Returns the ranked filtered table."""
    f = e[e["ipsae.tsa"] > ipsae_min].copy()
    z = pd.DataFrame({c: -(f[c] - f[c].mean()) / (f[c].std() + 1e-12) for c in SCORE_COLUMNS})       # lower raw = better = higher z
    r = pd.DataFrame({c: f[c].rank(ascending=False, pct=True) for c in SCORE_COLUMNS})
    f["n_metrics"] = z.notna().sum(axis=1)
    f["composite_z"] = np.where(f.n_metrics >= min_metrics, z.mean(axis=1), np.nan)
    f["composite_rank"] = np.where(f.n_metrics >= min_metrics, r.mean(axis=1), np.nan)
    for c in SCORE_COLUMNS:
        f["z_" + c] = z[c]
    f = f.sort_values("composite_z", ascending=False)
    f["final_rank"] = np.arange(1, len(f) + 1)
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--camp")
    ap.add_argument("--ipsae", type=float)
    a = ap.parse_args()
    cfg = config.load("gates")["final_selection"]
    ipsae_min = a.ipsae if a.ipsae is not None else cfg["ipsae_tsa_min"]
    root, tables = config.runs_root(a.camp), config.results_dir() / "tables"
    e = pd.read_csv(tables / "evaluated_sequences_full_suite.csv")
    e = e[~e.arm.str.contains("sheet")].copy()
    seqs = {}
    for p in glob.glob(f"{root}/*/*/b*/mpnn/sequences.csv"):
        s = pd.read_csv(p)
        seqs.update(dict(zip(s.seq_id, s.sequence)))
    terms = [per_run_terms(x, root, seqs) for x in e.itertuples()]
    e["S07_tmol_polar_residue_int_TSA_minus_SUB"] = [t["S07"] for t in terms]
    e["S12_boltz_cat_to_pocket_PAE_SUB_minus_PROD"] = [t["S12"] for t in terms]
    e = add_energy_terms(e)
    e["pass_ipsae"] = e["ipsae.tsa"] > ipsae_min
    f = composite(e, ipsae_min, cfg["min_composite_metrics"])
    f.to_csv(tables / "final_scoring.csv", index=False)
    e.to_csv(tables / "evaluated_with_final_metrics.csv", index=False)
    print("evaluated (non-sheet):", len(e), "| strict:", int((e.boltz_tier == 1).sum()), f"| ipSAE2>{ipsae_min}:", len(f),
          "| scored:", int(f.composite_z.notna().sum()))
    print("metric availability among ipSAE-passing:", {c: int(f[c].notna().sum()) for c in SCORE_COLUMNS})


if __name__ == "__main__":
    main()
