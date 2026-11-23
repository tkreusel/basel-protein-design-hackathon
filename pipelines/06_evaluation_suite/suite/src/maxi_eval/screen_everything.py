"""Screen EVERYTHING against activity, with selection-aware nulls per feature family.

  python -m maxi_eval.screen_everything --labels config/labels_manual_dataset.csv --out work/report/everything [--family HG3]

Candidate sets: ARM (per-arm + TSA/sub/prod differences from all four tools), EXTRA (features_extra.csv: PLACER ensemble distributions, tmol per-residue
energies, Boltz pocket PAE/pLDDT, composition), LEGACY (every other numeric column of the master table).
For each candidate: Spearman rho, BH-FDR over ALL candidates, bootstrap-free robustness (leave-one-protein-out range, within-family rho, rho under the OLD labels),
AUROC at 1e4/1e5/5e5, and (HG3 only) partial rho controlling for number of mutations vs HG3 (does it merely track how evolved the variant is?).
Selection-aware permutation p (best |rho| of the family vs labels permuted) for: everything, each set, each EXTRA family, and each tool.
"""
from __future__ import annotations
import argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
warnings.filterwarnings("ignore")
from .analyze_full import load_master, within_family_rho
from .exhaustive import build_candidates, rank_cols, perm_null
from .validate_dataset import bh, auroc, loo_composite


def set_of(c: str, meta: dict) -> str:
    if c in meta:
        return f"ARM-{meta[c][1]}"
    if c.startswith(("pl2.", "pr.", "bz.", "cmp.")):
        return f"EXTRA-{c.split('.')[0]}"
    return "LEGACY"


def perm_null_nan(X: pd.DataFrame, y: np.ndarray, groups: dict[str, np.ndarray], nperm=5000, seed=0):
    """Selection-aware null of max|rho| per group, with PAIRWISE-complete rho (columns may miss a few proteins)."""
    V = X.to_numpy(float); n, M = V.shape; mask = np.isfinite(V); full = mask.all(0)
    rng = np.random.default_rng(seed); perms = np.stack([rng.permutation(n) for _ in range(nperm)])     # indices into y
    A = np.zeros((nperm, M)); obs = np.zeros(M)
    fi = np.where(full)[0]
    if len(fi):
        R = rank_cols(V[:, fi]); ry = rankdata(y); r0 = (ry - ry.mean()) / np.linalg.norm(ry - ry.mean()); obs[fi] = np.abs(r0 @ R)
        Y = np.stack([rankdata(y[pm]) for pm in perms]); Y = Y - Y.mean(1, keepdims=True); Y = Y / np.linalg.norm(Y, axis=1, keepdims=True); A[:, fi] = np.abs(Y @ R)
    for j in np.where(~full)[0]:
        idx = np.where(mask[:, j])[0]
        if len(idx) < 8 or np.std(V[idx, j]) == 0:
            continue
        r = rankdata(V[idx, j]); r = (r - r.mean()) / np.linalg.norm(r)
        yy = y[idx]; ry = rankdata(yy); obs[j] = abs(((ry - ry.mean()) / np.linalg.norm(ry - ry.mean())) @ r)
        for k, pm in enumerate(perms):           # permute the labels among ALL proteins, then restrict to those with data
            yp = rankdata(y[pm][idx]); yp = yp - yp.mean(); nn = np.linalg.norm(yp)
            A[k, j] = abs((yp / nn) @ r) if nn > 0 else 0.0
    res = {}
    for g, idx in groups.items():
        if len(idx):
            o = obs[idx].max(); nm = A[:, idx].max(1)
            res[g] = (float(o), float((nm >= o).mean()), float(np.percentile(nm, 95)), cols_best(idx, obs))
    return res, obs


def cols_best(idx, obs):
    return int(idx[np.argmax(obs[idx])])


def partial(x, y, z):
    d = pd.concat([x, y, z], axis=1).dropna()
    if len(d) < 6 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 2].nunique() < 2:
        return np.nan
    r = [rankdata(d.iloc[:, i]) for i in range(3)]; A = np.c_[np.ones(len(d)), r[2]]
    ex = r[0] - A @ np.linalg.lstsq(A, r[0], rcond=None)[0]; ey = r[1] - A @ np.linalg.lstsq(A, r[1], rcond=None)[0]
    return float(np.corrcoef(ex, ey)[0, 1]) if ex.std() > 0 and ey.std() > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="work/report"); ap.add_argument("--work", default="work")
    ap.add_argument("--labels", default="config/labels_manual_dataset.csv"); ap.add_argument("--old-labels", default="config/labels.csv")
    ap.add_argument("--extra", default="work/features_extra.csv"); ap.add_argument("--out", default="work/report/everything"); ap.add_argument("--family", default=None)
    ap.add_argument("--nperm", type=int, default=4000)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    ent = load_master(Path(a.report), Path(a.work), a.labels)
    ARM = build_candidates(ent); meta = ARM.attrs["meta"]
    EX = pd.read_csv(a.extra, index_col=0) if Path(a.extra).exists() else pd.DataFrame(index=ent.index)
    skip = set(ARM.columns)
    LEG = ent.select_dtypes("number").drop(columns=[c for c in ent.columns if c in ("kcat_KM", "log10_kcat_KM") or c.endswith("|sd")], errors="ignore")
    LEG = LEG.loc[:, LEG.notna().sum() >= 12]
    ALL = pd.concat([ARM, EX, LEG], axis=1)
    ALL = ALL.loc[:, ~ALL.columns.duplicated()]
    lab = ent[ent.log10_kcat_KM.notna()]
    if a.family:
        lab = lab[lab.family == a.family]
    X = ALL.loc[lab.index]; y = lab.log10_kcat_KM.to_numpy(float); fam = lab.family.to_numpy(); n = len(y)
    old = pd.read_csv(a.old_labels) if Path(a.old_labels).exists() else pd.DataFrame(columns=["entry", "kcat_KM"])
    old_k = {r.entry.lower(): np.log10(float(r.kcat_KM)) for r in old.itertuples() if pd.notna(r.kcat_KM)}
    yo = np.array([old_k.get(p.lower(), np.nan) for p in lab.index])
    nmut = EX["cmp.n_mut_vs_HG3"].reindex(lab.index) if "cmp.n_mut_vs_HG3" in EX else pd.Series(np.nan, index=lab.index)
    rows = []
    for c in X.columns:
        v = X[c].to_numpy(float); ok = np.isfinite(v)
        if ok.sum() < max(10, n - 3) or np.std(v[ok]) == 0:
            continue
        rho, p = spearmanr(v[ok], y[ok]); idx = np.where(ok)[0]
        lo = [spearmanr(np.delete(v[ok], j), np.delete(y[ok], j))[0] for j in range(len(idx))]
        row = dict(candidate=c, set=set_of(c, meta), n=int(ok.sum()), rho=rho, p=p, rho_loo_min=np.nanmin(lo), rho_loo_max=np.nanmax(lo),
                   rho_within_family=within_family_rho(v[ok], y[ok], fam[ok]) if not a.family else np.nan)
        for thr in (1e4, 1e5, 5e5):
            act = (lab.kcat_KM.to_numpy(float)[ok] >= thr).astype(int)
            row[f"auroc_{thr:g}"] = auroc(v[ok], act) if 0 < act.sum() < len(act) else np.nan
        oo = ok & np.isfinite(yo)
        row["rho_old_labels"] = spearmanr(v[oo], yo[oo])[0] if oo.sum() >= 6 and np.std(v[oo]) > 0 else np.nan
        row["partial_given_nmut"] = partial(pd.Series(v, index=lab.index), pd.Series(y, index=lab.index), nmut)
        rows.append(row)
    S = pd.DataFrame(rows); S["q_BH"] = bh(S.p.to_numpy())
    S = S.reindex(S.rho.abs().sort_values(ascending=False).index); S.to_csv(out / "candidates_all.csv", index=False)
    ALL.to_csv(out / "candidate_matrix.csv")
    # ------------- selection-aware nulls
    Xv = X.loc[:, X.notna().sum() >= max(10, n - 3)]
    Xv = Xv.loc[:, Xv.nunique() > 1]; cc = list(Xv.columns)
    sets = pd.Series([set_of(c, meta) for c in cc])
    groups = {"ALL": np.arange(len(cc))}
    for s_ in sets.unique():
        groups[f"set:{s_}"] = np.where(sets == s_)[0]
    groups["set:EXTRA(all)"] = np.where(sets.str.startswith("EXTRA"))[0]
    groups["set:ARM(all)"] = np.where(sets.str.startswith("ARM"))[0]
    nul, obs = perm_null_nan(Xv, y, groups, a.nperm)
    N = pd.DataFrame([dict(group=g, n_candidates=len(groups[g]), best_abs_rho=v[0], p_selection_aware=v[1], null95=v[2], best_candidate=cc[v[3]]) for g, v in nul.items()])
    N.to_csv(out / "selection_aware_by_set.csv", index=False)
    Xc = X.dropna(axis=1); cc2 = list(Xc.columns)
    comp = {k: loo_composite(Xc, y, k)[0] for k in (1, 3, 5)}
    cpool = {k: loo_composite(Xc[[c for c in cc2 if not c.startswith(("pl2.", "pr.", "bz.", "cmp."))]], y, k)[0] for k in (1, 3)}
    md = [f"# Screen of everything ({'family ' + a.family if a.family else 'all proteins'}, n={n}; {len(S)} candidates screened; {len(cc)} valid for the permutation tests)\n",
          f"BH q<0.1: **{int((S.q_BH < 0.1).sum())}**; q<0.2: {int((S.q_BH < 0.2).sum())}. LOO composite with selection inside the fold: " + ", ".join(f"top-{k}: {v:+.2f}" for k, v in comp.items()) +
          " | excluding the new EXTRA features: " + ", ".join(f"top-{k}: {v:+.2f}" for k, v in cpool.items()) + "\n",
          "## Selection-aware test by feature set\n```", N.round(3).to_string(index=False), "```\n"]
    cols = ["candidate", "n", "rho", "p", "q_BH", "rho_within_family", "rho_loo_min", "rho_loo_max", "auroc_100000", "rho_old_labels", "partial_given_nmut"]
    md += ["## Top 40 overall\n```", S.head(40)[cols].round(2).to_string(index=False), "```\n"]
    for s in sorted(S["set"].unique()):
        md += [f"## Top 8 in {s}\n```", S[S["set"] == s].head(8)[cols].round(2).to_string(index=False), "```\n"]
    (out / "everything.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[:4])); print("->", out)

if __name__ == "__main__":
    main()
