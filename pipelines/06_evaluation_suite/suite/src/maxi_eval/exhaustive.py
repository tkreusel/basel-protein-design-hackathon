"""Exhaustive, selection-aware screen of ALL per-arm metrics and ALL arm differences against activity.

  python -m maxi_eval.exhaustive --labels config/labels_manual_dataset.csv --out work/report/exhaustive [--family HG3]

Candidates = for every metric that exists for all three arms in a tool (ipSAE2/Boltz, PLACER holo, tmol, Rosetta; native refolds and same-pocket swap):
  the value in each arm (tsa, sub, prod)  +  the three differences (tsa-sub, tsa-prod, sub-prod).
Per candidate: Spearman rho vs log10 kcat/KM, bootstrap CI, BH-FDR (over ALL candidates), within-family rho, AUROC at several activity thresholds,
leave-one-protein-out range of rho, rho under the OLD label set (robustness to label source).
Selection-aware tests: the null distribution of the best |rho| (labels permuted) is computed for the whole candidate set AND separately for each
tool and for each arm-type, so "is there signal in the PLACER TSA-product differences?" has its own honest answer.
Also: leave-one-out composites with metric selection inside the fold.
"""
from __future__ import annotations
import argparse, re, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
warnings.filterwarnings("ignore")
from .analyze_full import load_master, within_family_rho
from .validate_dataset import bh, auroc, loo_composite

FAMS = {  # tool-family -> (tsa prefix, sub prefix, prod prefix)
    "ipSAE2/Boltz": ("ip.tsa.", "ip.sub.", "ip.prod."),
    "PLACER holo": ("placer.holo.", "placer_sub.holo.", "placer_prod.holo."),
    "tmol native": ("rlx.", "native_sub.", "native_prod."),
    "tmol swap": ("tsa_pk.", "sub.", "prod."),
    "Rosetta native": ("ros.tsa.", "ros.native_sub.", "ros.native_prod."),
    "Rosetta swap": ("ros.tsa.", "ros.sub.", "ros.prod."),
}
TOOL = lambda fam: fam.split("/")[0].split(" ")[0]          # ipSAE2 / PLACER / tmol / Rosetta
SKIP = re.compile(r"(^|\.)(n_models|place_rmsd|status|n_int_res|n_res|time)|per_res_energy|atoms$|resid$")
TYPES = ["tsa", "sub", "prod", "tsa-sub", "tsa-prod", "sub-prod"]


def build_candidates(ent: pd.DataFrame) -> pd.DataFrame:
    cols, meta = {}, {}
    for fam, (pt, ps, pp) in FAMS.items():
        suffixes = [c[len(pt):] for c in ent.columns if c.startswith(pt)]
        for s in dict.fromkeys(suffixes):
            if SKIP.search(s):
                continue
            ct, cs, cp = pt + s, ps + s, pp + s
            if not all(c in ent.columns for c in (ct, cs, cp)):
                continue
            t, su, pr = (ent[c].astype(float) for c in (ct, cs, cp))
            if min(t.notna().sum(), su.notna().sum(), pr.notna().sum()) < 15:
                continue
            for typ, v in (("tsa", t), ("sub", su), ("prod", pr), ("tsa-sub", t - su), ("tsa-prod", t - pr), ("sub-prod", su - pr)):
                name = f"{fam}|{s}|{typ}"
                if v.nunique() > 1:
                    cols[name] = v; meta[name] = (fam, TOOL(fam), s, typ)
    X = pd.DataFrame(cols)
    X.attrs["meta"] = meta
    return X


def rank_cols(X: np.ndarray) -> np.ndarray:
    R = np.apply_along_axis(rankdata, 0, X); R = R - R.mean(0)
    return R / np.maximum(np.linalg.norm(R, axis=0), 1e-12)


def perm_null(R: np.ndarray, y: np.ndarray, groups: dict[str, np.ndarray], nperm=10000, seed=0):
    """observed + null of max|rho| for the whole set and each group (indices into R's columns)."""
    ry = rankdata(y); ry0 = (ry - ry.mean()) / np.linalg.norm(ry - ry.mean())
    obs_all = np.abs(ry0 @ R)
    rng = np.random.default_rng(seed)
    P = np.stack([rng.permutation(ry) for _ in range(nperm)]); P = P - P.mean(1, keepdims=True); P = P / np.linalg.norm(P, axis=1, keepdims=True)
    A = np.abs(P @ R)                              # nperm x M
    res = {"ALL": (obs_all.max(), float((A.max(1) >= obs_all.max()).mean()), float(np.percentile(A.max(1), 95)))}
    for g, idx in groups.items():
        if len(idx):
            o = obs_all[idx].max(); nm = A[:, idx].max(1)
            res[g] = (float(o), float((nm >= o).mean()), float(np.percentile(nm, 95)))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="work/report"); ap.add_argument("--work", default="work")
    ap.add_argument("--labels", default="config/labels_manual_dataset.csv"); ap.add_argument("--old-labels", default="config/labels.csv")
    ap.add_argument("--out", default="work/report/exhaustive"); ap.add_argument("--family", default=None)
    ap.add_argument("--nperm", type=int, default=10000)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    ent = load_master(Path(a.report), Path(a.work), a.labels)
    old = pd.read_csv(a.old_labels) if Path(a.old_labels).exists() else pd.DataFrame(columns=["entry", "kcat_KM"])
    old_k = {r.entry.lower(): np.log10(float(r.kcat_KM)) for r in old.itertuples() if pd.notna(r.kcat_KM)}
    lab = ent[ent.log10_kcat_KM.notna()]
    if a.family:
        lab = lab[lab.family == a.family]
    X = build_candidates(ent).loc[lab.index]; meta = build_candidates(ent).attrs["meta"]
    y = lab.log10_kcat_KM.to_numpy(float); fam = lab.family.to_numpy(); n = len(y)
    # --- univariate screen (columns with <=2 missing: pairwise rho; permutation uses complete-case columns)
    rows = []
    for c in X.columns:
        v = X[c].to_numpy(float); ok = np.isfinite(v)
        if ok.sum() < max(10, n - 3) or np.std(v[ok]) == 0:
            continue
        rho, p = spearmanr(v[ok], y[ok])
        # leave-one-protein-out range
        idx = np.where(ok)[0]; lo = [spearmanr(np.delete(v[ok], j), np.delete(y[ok], j))[0] for j in range(len(idx))]
        f, t_, s_, ty = meta[c]
        row = dict(candidate=c, tool=t_, family=f, metric=s_, type=ty, n=int(ok.sum()), rho=rho, p=p, rho_loo_min=np.nanmin(lo), rho_loo_max=np.nanmax(lo),
                   rho_within_family=within_family_rho(v[ok], y[ok], fam[ok]) if not a.family else np.nan)
        for thr in (1e4, 1e5, 5e5):
            act = (lab.kcat_KM.to_numpy(float)[ok] >= thr).astype(int)
            row[f"auroc_{thr:g}"] = auroc(v[ok], act) if 0 < act.sum() < len(act) else np.nan
        # robustness to the OLD label set (subset of proteins that have an old value)
        yo = np.array([old_k.get(p_.lower(), np.nan) for p_ in lab.index]); oo = ok & np.isfinite(yo)
        row["rho_old_labels"] = spearmanr(v[oo], yo[oo])[0] if oo.sum() >= 6 and np.std(v[oo]) > 0 else np.nan
        row["n_old"] = int(oo.sum())
        rows.append(row)
    S = pd.DataFrame(rows)
    S["q_BH"] = bh(S.p.to_numpy())
    S = S.reindex(S.rho.abs().sort_values(ascending=False).index)
    S.to_csv(out / "candidates_all.csv", index=False)
    # --- selection-aware nulls: whole set, per tool, per arm-type, per (tool, type)
    Xc = X.dropna(axis=1); cc = list(Xc.columns)
    R = rank_cols(Xc.to_numpy(float))
    ix = lambda f: np.array([i for i, c in enumerate(cc) if f(meta[c])])
    groups = {f"tool:{t}": ix(lambda m, t=t: m[1] == t) for t in ("ipSAE2", "PLACER", "tmol", "Rosetta")}
    groups.update({f"type:{ty}": ix(lambda m, ty=ty: m[3] == ty) for ty in TYPES})
    groups.update({f"{t}|{ty}": ix(lambda m, t=t, ty=ty: m[1] == t and m[3] == ty) for t in ("ipSAE2", "PLACER", "tmol", "Rosetta") for ty in TYPES})
    nul = perm_null(R, y, groups, a.nperm)
    N = pd.DataFrame([dict(group=g, n_candidates=(len(cc) if g == "ALL" else len(groups[g])), best_abs_rho=v[0], p_best_selection_aware=v[1], null95=v[2]) for g, v in nul.items()])
    # best metric of each group
    best = {}
    for g, idx in groups.items():
        if len(idx):
            k = idx[np.argmax(np.abs(rankdata(y) @ R[:, idx]))] if False else idx[np.argmax(np.abs(((rankdata(y) - rankdata(y).mean()) / np.linalg.norm(rankdata(y) - rankdata(y).mean())) @ R[:, idx]))]
            best[g] = cc[k]
    best["ALL"] = cc[int(np.argmax(np.abs(((rankdata(y) - rankdata(y).mean()) / np.linalg.norm(rankdata(y) - rankdata(y).mean())) @ R)))]
    N["best_candidate"] = N.group.map(best)
    N.to_csv(out / "selection_aware_by_group.csv", index=False)
    # --- composites with selection inside the fold
    comp = {}
    for k in (1, 3, 5):
        comp[k] = loo_composite(Xc, y, k)[0]
    pd.Series(comp, name="loo_spearman").to_csv(out / "loo_composites.csv")
    # --- report
    md = [f"# Exhaustive screen ({'family ' + a.family if a.family else 'all proteins'}; n={n}; {len(S)} candidates screened, {len(cc)} complete-case for the permutation tests)\n",
          f"BH q < 0.1: **{int((S.q_BH < 0.1).sum())}**; q < 0.2: {int((S.q_BH < 0.2).sum())}.  Leave-one-out composites (metric selection inside the fold): "
          + ", ".join(f"top-{k}: {v:+.2f}" for k, v in comp.items()) + "\n",
          "## Selection-aware tests (labels permuted; p = P(noise best |rho| >= observed best in that group))\n", "```", N.round(3).to_string(index=False), "```\n",
          "## Top 30 candidates overall\n", "```", S.head(30)[["candidate", "n", "rho", "p", "q_BH", "rho_within_family", "rho_loo_min", "rho_loo_max", "auroc_100000", "rho_old_labels"]].round(2).to_string(index=False), "```\n"]
    for ty in ("tsa-prod", "tsa-sub", "sub-prod"):
        sub = S[S.type == ty].head(8)
        md += [f"## Best {ty} difference candidates\n", "```", sub[["candidate", "n", "rho", "p", "q_BH", "rho_within_family", "rho_loo_min", "rho_loo_max", "rho_old_labels"]].round(2).to_string(index=False), "```\n"]
    for t in ("ipSAE2", "PLACER", "tmol", "Rosetta"):
        sub = S[S.tool == t].head(6)
        md += [f"## Best {t} candidates (any arm / difference)\n", "```", sub[["candidate", "n", "rho", "p", "q_BH", "rho_within_family", "rho_old_labels"]].round(2).to_string(index=False), "```\n"]
    (out / "exhaustive.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[:6])); print("->", out)

if __name__ == "__main__":
    main()
