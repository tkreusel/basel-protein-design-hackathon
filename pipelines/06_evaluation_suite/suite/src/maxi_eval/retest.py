"""Test-retest of the tmol protocol: run 1 (work/tmol) vs an independent run 2 (work/tmol_rep2) on the same Boltz structures.

  python -m maxi_eval.retest [--a work/tmol] [--b work/tmol_rep2] [--labels config/labels_manual_dataset.csv]

For each protein: the key quantities (swap dE_total/dint TSA-minus-substrate, TSA-minus-product, native interaction differences, strain, ...)
are computed from each run. Reports (1) test-retest rho / ICC across proteins, (2) rho vs kcat/KM for run 1, run 2 and their MEAN.
A metric whose test-retest rho is low cannot carry a real correlation; a metric that is reproducible AND correlates in both runs is trustworthy.
"""
from __future__ import annotations
import argparse, glob, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr


def load(d: str) -> pd.DataFrame:
    rows = {}
    for f in glob.glob(f"{d}/*.json"):
        s = json.load(open(f)).get("scalars", {})
        if s.get("status") == "ok":
            rows[s["entry"].replace("_boltz", "")] = {k: v for k, v in s.items() if isinstance(v, (int, float))}
    return pd.DataFrame(rows).T


def quantities(t: pd.DataFrame) -> pd.DataFrame:
    q = pd.DataFrame(index=t.index)
    for term in ("E_total", "int.total", "E_protein_in_complex", "int.fa_elec", "int.hbond", "int.fa_lk", "int.lk_ball_iso", "pocket_E_per_res"):
        for lab, a, b in (("tsa-sub", "tsa_pk", "sub"), ("tsa-prod", "tsa_pk", "prod")):
            if f"{a}.{term}" in t and f"{b}.{term}" in t:
                q[f"swap {term} {lab}"] = t[f"{a}.{term}"] - t[f"{b}.{term}"]
    for term in ("E_total_noncst", "int.total", "int.hbond", "int.fa_elec"):
        for lab, b in (("tsa-sub", "native_sub"), ("tsa-prod", "native_prod")):
            if f"rlx.{term}" in t and f"{b}.{term}" in t:
                q[f"native {term} {lab}"] = t[f"rlx.{term}"] - t[f"{b}.{term}"]
    for c in ("strain.E_noncst", "apo.reorg_E", "rlx.E_total_per_res", "rlx.int.total"):
        if c in t:
            q[c] = t[c]
    return q


def icc(a: np.ndarray, b: np.ndarray) -> float:   # ICC(3,1), consistency
    x = np.c_[a, b]; n = len(x); gm = x.mean(); ms_r = 2 * ((x.mean(1) - gm) ** 2).sum() / (n - 1)
    ms_e = ((x - x.mean(1, keepdims=True) - x.mean(0, keepdims=True) + gm) ** 2).sum() / (n - 1)
    return float((ms_r - ms_e) / (ms_r + ms_e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", default="work/tmol"); ap.add_argument("--b", default="work/tmol_rep2"); ap.add_argument("--labels", default="config/labels_manual_dataset.csv")
    ap.add_argument("--out", default="work/report/retest.csv")
    a = ap.parse_args()
    A, B = quantities(load(a.a)), quantities(load(a.b))
    common = A.index.intersection(B.index)
    lab = pd.read_csv(a.labels); k = {r.entry: np.log10(float(r.kcat_KM)) for r in lab.itertuples()}
    y = pd.Series({p: k.get(p, np.nan) for p in common})
    rows = []
    for c in A.columns.intersection(B.columns):
        x1, x2 = A.loc[common, c], B.loc[common, c]; ok = x1.notna() & x2.notna() & y.notna()
        if ok.sum() < 8:
            continue
        rt = spearmanr(x1[ok], x2[ok])[0]
        rows.append(dict(quantity=c, n=int(ok.sum()), retest_rho=rt, ICC=icc(x1[ok].to_numpy(), x2[ok].to_numpy()), mean_abs_diff=float((x1 - x2).abs().mean()),
                         sd_across_proteins=float(x1.std()), rho_kcat_run1=spearmanr(x1[ok], y[ok])[0], rho_kcat_run2=spearmanr(x2[ok], y[ok])[0],
                         rho_kcat_mean=spearmanr(((x1 + x2) / 2)[ok], y[ok])[0]))
    R = pd.DataFrame(rows).sort_values("rho_kcat_mean", key=abs, ascending=False)
    R.to_csv(a.out, index=False); pd.set_option("display.width", 220)
    print(f"{len(common)} proteins in both runs"); print(R.round(2).to_string(index=False))

if __name__ == "__main__":
    main()
