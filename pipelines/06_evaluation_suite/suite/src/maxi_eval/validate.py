"""Validation table: do the in-silico metrics order the enzymes the way measured activity does?

Reads work/report/summary_replicates.csv (aggregate.py output). Writes work/report/validation.md + validation_*.csv.

Three views, each deliberately modest about what n permits:
  A  per-entry table of headline metrics next to kcat/KM
  B  LADDER   HG3 (7K4Q+5RGA, 146) -> HG3.7 (7K4X, 27,000) -> HG3.17 (7K4Z, 126,000): is each metric monotone with activity?
              3 enzymes: this is a sanity check on direction, not a statistical test (no p-values are quoted)
  C  CLASS    low (HG3: no Gln50 donor) vs evolved (HG3.7, HG3.17, HG3.R5, HG3.17 E47N/N300D) using every crystal copy:
              mean difference, and whether the classes are completely separated. The class labels are inferred from lineage
              (see CLASSES) and are NOT measured values for 8RD5 / 4BS0.
Expected direction ('+' higher is better) is declared per metric BEFORE looking at data, in EXPECT.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd

# metric -> expected direction for a MORE ACTIVE enzyme (+1 higher, -1 lower)
EXPECT = {
    "rlx.geo.n_donors_3p5": +1, "rlx.geo.n_sidechain_donors_3p5": +1,
    "placer.holo.frac_pass_core": +1, "placer.apo.frac_pass_core": +1, "placer.holo.frac_donor": +1,
    "placer.holo.lig_rmsd_spread": -1, "placer.holo.chi_sd_base": -1, "placer.holo.chi_sd_donor": -1, "placer.holo.prmsd_mean": -1,
    "placer.holo.field_CH_sd": -1, "placer.holo.field_LG_sd": -1,
    "rlx.cst.n_pass": +1, "rlx.cst.rms_z": -1, "strain.E_noncst": -1,
    "delta.tsa_minus_sub.total": -1, "delta.prod_minus_sub.total": +1, "delta.tsa_minus_sub.hbond": -1,
    "rlx.int.total": -1, "rlx.int.hbond": -1, "rlx.rep_intra_pocket": -1, "rlx.E_total_per_res": -1, "apo.reorg_E": -1,
    "rlx.buns.n_unsat": -1, "rlx.void.volume_A3": -1,
    "rlx.field.CH.env_nobase.vac": +1, "rlx.field.LG.env.vac": +1,
    "rlx.geo.wsite_open_total": -1,
    "ros.tsa.dG_sep": -1, "ros.delta.tsa_minus_sub.dG_sep": -1, "ros.delta.prod_minus_sub.dG_sep": +1,
    "delta_native.tsa_minus_sub.total": -1, "delta_native.prod_minus_sub.total": +1, "ros.delta_native.tsa_minus_sub.dG_sep": -1,
    "ros.delta_native.prod_minus_sub.dG_sep": +1, "ipsae.delta_tsa_sub": +1,
    "ros.tsa.sc": +1, "ros.tsa.packstat_int": +1, "ros.tsa.filter.buns_all": -1, "ros.tsa.delta_unsat_hb": -1,
}
LADDER = [("HG3", ["7K4Q", "5RGA"], 146), ("HG3.7", ["7K4X"], 27000), ("HG3.17", ["7K4Z"], 126000)]
CLASSES = {"low (HG3)": ["7K4Q", "5RGA"], "evolved": ["7K4X", "7K4Z", "8RD5", "4BS0"]}
KEY = ["rlx.geo.n_donors_3p5", "rlx.geo.d_O_X", "rlx.geo.ang_O_H_X", "rlx.cst.n_pass", "placer.holo.frac_pass_core", "placer.apo.frac_pass_core",
       "placer.holo.lig_rmsd_spread", "delta.tsa_minus_sub.total", "ros.delta.tsa_minus_sub.dG_sep", "delta.prod_minus_sub.total",
       "rlx.rep_intra_pocket", "rlx.E_total_per_res", "rlx.field.LG.env.vac", "ros.tsa.dG_sep"]


def fmt(x):
    return "-" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="work/report")
    a = ap.parse_args()
    R = Path(a.report)
    rep = pd.read_csv(R / "summary_replicates.csv")
    ent = pd.read_csv(R / "summary_entries.csv", index_col=0)
    metrics = [m for m in EXPECT if m in rep.columns]
    out = ["# Validation against measured activity\n",
           "Heuristic features; 'agree' = the metric moves in the direction declared in advance (EXPECT) for a more active enzyme.\n"]

    # ---- A: per-entry
    cols = [c for c in KEY if c in ent.columns]
    A = ent[cols + ["kcat_KM"]].copy()
    A = A.loc[[e for e in A.index if not e.endswith("_boltz")] + [e for e in A.index if e.endswith("_boltz")]]
    out += ["## A. Per-entry means (crystal/NMR copies, then Boltz refolds)\n", A.round(2).to_string(), "\n"]
    A.to_csv(R / "validation_entries.csv")

    # ---- B: ladder
    rows = []
    for m in metrics:
        vals = []
        for name, ents, k in LADDER:
            v = rep[rep.entry.isin(ents)][m].astype(float)
            vals.append(v.mean() if v.notna().any() else np.nan)
        if np.isnan(vals).any():
            rows.append([m, *vals, EXPECT[m], "n/a (missing)"]); continue
        inc = vals[0] <= vals[1] <= vals[2]
        dec = vals[0] >= vals[1] >= vals[2]
        mono = "increasing" if inc and not dec else "decreasing" if dec and not inc else "flat" if inc and dec else "non-monotone"
        agree = "AGREE" if (mono == "increasing" and EXPECT[m] > 0) or (mono == "decreasing" and EXPECT[m] < 0) else \
                ("n/a" if mono == "flat" else "disagree")
        rows.append([m, *vals, EXPECT[m], f"{mono} -> {agree}"])
    B = pd.DataFrame(rows, columns=["metric", "HG3 (146)", "HG3.7 (27,000)", "HG3.17 (126,000)", "expected", "verdict"])
    out += ["## B. Activity ladder: HG3 -> HG3.7 -> HG3.17 (3 enzymes; direction check only)\n", B.round(2).to_string(index=False), "\n"]
    B.to_csv(R / "validation_ladder.csv", index=False)
    full = B[B.verdict.str.contains("->")]
    out.append(f"Ladder: {int(full.verdict.str.endswith('AGREE').sum())} of {len(full)} metrics with complete data are monotone in the expected direction.\n")

    # ---- C: class separation using every copy
    low = rep[rep.entry.isin(CLASSES["low (HG3)"])]
    hi = rep[rep.entry.isin(CLASSES["evolved"])]
    rows = []
    for m in metrics:
        l, h = low[m].astype(float).dropna(), hi[m].astype(float).dropna()
        if len(l) < 2 or len(h) < 2:
            continue
        d = h.mean() - l.mean()
        sep = (h.min() > l.max()) if d > 0 else (h.max() < l.min())
        good = np.sign(d) == EXPECT[m]
        rows.append([m, len(l), len(h), l.mean(), h.mean(), d, "yes" if sep else "no", "AGREE" if good else "disagree"])
    C = pd.DataFrame(rows, columns=["metric", "n_low", "n_evolved", "mean_low", "mean_evolved", "diff", "fully_separated", "direction"])
    out += ["## C. Low (HG3) vs evolved (HG3.7, HG3.17, HG3.R5, HG3.17 E47N/N300D), every copy\n",
            "(class membership of 8RD5 and 4BS0 is inferred from lineage, not a measured kcat/KM)\n", C.round(2).to_string(index=False), "\n"]
    C.to_csv(R / "validation_classes.csv", index=False)
    if len(C):
        out.append(f"Classes: {int((C.direction == 'AGREE').sum())} of {len(C)} metrics point the expected way; {int((C.fully_separated == 'yes').sum())} separate the classes completely.\n")

    # ---- D: crystal vs Boltz refold agreement (does the predicted pose behave like the crystal?)
    pairs = [(e, e + "_boltz") for e in ["7K4Q", "5RGA", "8RD5", "4BS0", "29SB"] if e + "_boltz" in ent.index and e in ent.index]
    if pairs:
        rows = []
        for m in [c for c in KEY if c in ent.columns]:
            rows.append([m] + [fmt(ent.loc[x, m]) for p in pairs for x in p])
        D = pd.DataFrame(rows, columns=["metric"] + [x for p in pairs for x in p])
        out += ["## D. Crystal vs Boltz refold of the same entry\n", D.to_string(index=False), "\n"]
    out += ["## How to read this\n",
            "* B has 3 enzymes: it can show that a metric is *compatible* with the activity ordering, not that it predicts it.",
            "* C uses 6 enzymes (2 copies for most) but class labels for two are inferred; copies of one enzyme are not independent.",
            "* tmol relax is stochastic (several REU run to run): differences smaller than that are not evidence.",
            "* A metric that AGREEs in B and C and is not a trivial proxy for 'has a Gln50 donor' is the interesting one. Check `rlx.geo.n_donors_3p5` first:",
            "  if most metrics only reproduce it, they add nothing beyond the known donor determinant."]
    (R / "validation.md").write_text("\n".join(out) + "\n")
    print("\n".join(out[:14])); print("->", R / "validation.md")

if __name__ == "__main__":
    main()
