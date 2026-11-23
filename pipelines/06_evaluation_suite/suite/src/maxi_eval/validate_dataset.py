"""Discriminative power on a larger labelled dataset (one row per protein).

  python -m maxi_eval.validate_dataset [--group boltz|crystal|all] [--active-threshold 1000] [--control rlx.geo.n_donors_3p5]

Unit of analysis = protein (entry): replicates/copies are averaged first, so copies of one enzyme never count as independent points.
For every numeric metric with enough labelled, non-constant data:
  * Spearman rho vs log10(kcat/KM), permutation-free p, Benjamini-Hochberg q across ALL metrics tested, bootstrap 95% CI of rho
  * AUROC for 'active' (kcat/KM >= threshold, or class == active) vs the rest
  * partial Spearman controlling for a nuisance metric (default: number of oxyanion donors) -> does it add anything beyond the known determinant?
  * LEAVE-ONE-OUT COMPOSITE: for each protein, pick the top-k metrics on the OTHER proteins (selection inside the loop), z-score and
    sign them, average -> predicted rank; Spearman of these held-out predictions vs truth is the honest estimate of what a metric-selected
    ranking can do (selection from hundreds of metrics inflates in-sample rho badly; this does not).
With few proteins (< 12) or many metrics the report says so and refuses to rank. No p-value is meaningful without the BH column.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata


def bh(p: np.ndarray) -> np.ndarray:
    n = len(p); o = np.argsort(p); q = np.empty(n)
    prev = 1.0
    for i, idx in enumerate(o[::-1]):
        rank = n - i
        prev = min(prev, p[idx] * n / rank); q[idx] = prev
    return q


def auroc(x: np.ndarray, y: np.ndarray) -> float:
    r = rankdata(x); n1 = int(y.sum()); n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def partial_spearman(x, y, z) -> float:
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    def res(a, b):
        A = np.c_[np.ones_like(b), b]
        return a - A @ np.linalg.lstsq(A, a, rcond=None)[0]
    ex, ey = res(rx, rz), res(ry, rz)
    if ex.std() == 0 or ey.std() == 0:
        return float("nan")
    return float(np.corrcoef(ex, ey)[0, 1])


def screen(X: pd.DataFrame, y: np.ndarray, active: np.ndarray | None, control: np.ndarray | None, n_boot=500, seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed); rows = []
    for m in X.columns:
        v = X[m].to_numpy(float); ok = np.isfinite(v) & np.isfinite(y)
        if ok.sum() < 8 or np.nanstd(v[ok]) == 0:
            continue
        rho, p = spearmanr(v[ok], y[ok])
        idx = np.where(ok)[0]; bs = []
        for _ in range(n_boot):
            b = rng.choice(idx, len(idx))
            if np.std(v[b]) > 0 and np.std(y[b]) > 0:
                bs.append(spearmanr(v[b], y[b])[0])
        lo, hi = (np.percentile(bs, [2.5, 97.5]) if bs else (np.nan, np.nan))
        row = dict(metric=m, n=int(ok.sum()), rho=rho, p=p, rho_lo=lo, rho_hi=hi)
        if active is not None and np.isfinite(active[ok]).all() and 0 < active[ok].sum() < ok.sum():
            row["auroc"] = auroc(v[ok], active[ok].astype(int))
        if control is not None:
            c = np.isfinite(control) & ok
            row["partial_rho"] = partial_spearman(v[c], y[c], control[c]) if c.sum() >= 8 else np.nan
        rows.append(row)
    d = pd.DataFrame(rows)
    if len(d):
        d["q_BH"] = bh(d.p.to_numpy())
        d = d.reindex(d.rho.abs().sort_values(ascending=False).index)
    return d


def loo_composite(X: pd.DataFrame, y: np.ndarray, k: int = 3) -> tuple[float, np.ndarray]:
    """Honest held-out estimate: metric selection happens inside the loop."""
    n = len(y); pred = np.full(n, np.nan)
    Z = (X - X.mean()) / X.std(ddof=0).replace(0, np.nan)
    for i in range(n):
        tr = np.array([j for j in range(n) if j != i])
        scores = {}
        for m in X.columns:
            v = X[m].to_numpy(float)[tr]; ok = np.isfinite(v)
            if ok.sum() >= 8 and np.nanstd(v[ok]) > 0:
                scores[m] = spearmanr(v[ok], y[tr][ok])[0]
        top = sorted(scores, key=lambda m: -abs(scores[m]) if np.isfinite(scores[m]) else 0)[:k]
        if top:
            vals = [np.sign(scores[m]) * Z[m].iloc[i] for m in top]
            pred[i] = np.nanmean(vals)
    ok = np.isfinite(pred)
    return (float(spearmanr(pred[ok], y[ok])[0]) if ok.sum() >= 5 else float("nan")), pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="work/report")
    ap.add_argument("--config", default="config/entries.yaml")
    ap.add_argument("--group", default="boltz", choices=["boltz", "crystal", "all"])
    ap.add_argument("--active-threshold", type=float, default=None, help="kcat/KM >= this counts as active (default: median)")
    ap.add_argument("--panel", default="curated", choices=["curated", "all"],
                    help="curated = ~40 metrics declared in advance (validate.EXPECT + ipSAE/site columns); all = every numeric column (exploratory, huge multiplicity)")
    ap.add_argument("--control", default="rlx.geo.n_donors_3p5")
    ap.add_argument("--min-proteins", type=int, default=12)
    ap.add_argument("--topk", type=int, default=3)
    a = ap.parse_args()
    R = Path(a.report)
    ent = pd.read_csv(R / "summary_entries.csv", index_col=0)
    if a.group == "boltz":
        ent = ent[ent.index.str.endswith("_boltz")]
    elif a.group == "crystal":
        ent = ent[~ent.index.str.endswith("_boltz")]
    num = ent.select_dtypes("number").drop(columns=[c for c in ent.columns if c.endswith("|sd") or c in ("kcat_KM", "log10_kcat_KM")], errors="ignore")
    if a.panel == "curated":
        from .validate import EXPECT
        keep = set(EXPECT) | {c for c in num.columns if c.startswith(("ipsae.", "site.pae", "site.plddt", "boltz.")) and not c.endswith("|sd")} \
               | {c for c in num.columns if c.startswith(("placer_sub.", "placer_prod.")) and c.endswith(("frac_pass_core", "prmsd_mean", "lig_rmsd_spread"))}
        num = num[[c for c in num.columns if c in keep]]
    lab = ent[ent.log10_kcat_KM.notna()]
    out = [f"# Discriminative power ({a.group} proteins)\n", f"{len(ent)} proteins, {len(lab)} with a kcat/KM label, {num.shape[1]} metrics.\n"]
    if len(lab) < a.min_proteins:
        out.append(f"**Not ranked: only {len(lab)} labelled proteins (need >= {a.min_proteins}).** With this n a Spearman rho needs |rho| > ~0.7 to reach p<0.05 for ONE metric "
                   f"and ~{num.shape[1]} metrics are being screened, so any 'winner' would be noise. Add labels to config/labels.csv.\n")
        (R / "validation_dataset.md").write_text("\n".join(out)); print("\n".join(out)); return
    y = lab.log10_kcat_KM.to_numpy(float); X = num.loc[lab.index]
    thr = np.log10(a.active_threshold) if a.active_threshold else np.median(y)
    active = (y >= thr).astype(float)
    ctrl = X[a.control].to_numpy(float) if a.control in X else None
    d = screen(X, y, active, ctrl)
    d.to_csv(R / "validation_dataset.csv", index=False)
    nsig = int((d.q_BH < 0.1).sum()) if len(d) else 0
    out += [f"Active = log10 kcat/KM >= {thr:.2f} ({int(active.sum())} active / {int((1-active).sum())} inactive). Control for partial rho: `{a.control}`.\n",
            f"Metrics with BH q < 0.1: **{nsig}** of {len(d)} screened (expect ~{0.1*len(d):.0f} false positives if none were real at q=0.1... interpret q, not p).\n",
            "## Top metrics\n", d.head(25).round(3).to_string(index=False), "\n"]
    r, pred = loo_composite(X, y, a.topk)
    out += [f"## Leave-one-out composite (top-{a.topk} metrics chosen inside each fold)\n",
            f"held-out Spearman = **{r:.2f}** (n={len(y)}). This, not the in-sample table above, is the honest estimate of ranking power from selected metrics.\n"]
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        top = d.head(6).metric.tolist(); fig, ax = plt.subplots(2, 3, figsize=(11, 6))
        for axx, m in zip(ax.ravel(), top):
            axx.scatter(X[m], y, s=18); axx.set_title(m[:34], fontsize=8); axx.set_ylabel("log10 kcat/KM", fontsize=7)
        fig.tight_layout(); fig.savefig(R / "validation_dataset_top.png", dpi=130)
    except Exception as exc:  # noqa: BLE001
        out.append(f"(plot skipped: {exc})")
    (R / "validation_dataset.md").write_text("\n".join(out) + "\n"); print("\n".join(out[:12]))

if __name__ == "__main__":
    main()
