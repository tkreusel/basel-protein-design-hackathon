"""Full analysis of a labelled protein set: (A) what correlates with kcat/KM, (B) how TSA / substrate / product complexes differ.

  python -m maxi_eval.analyze_full [--report work/report] [--out work/report/analysis]

Unit = protein (replicates averaged).  Inputs (all optional except summary_entries.csv): work/report/summary_entries.csv, work/ipsae/ipsae_models.csv,
work/apo_holo.csv (apo vs holo), work/arm_poses.csv (pose differences between arms), config/labels.csv.

PART A  correlation with log10 kcat/KM
  * Spearman rho, BH-FDR q over the pre-declared PANEL, bootstrap CI, AUROC (active vs inactive at --active-threshold)
  * WITHIN-FAMILY rho: ranks are taken inside each scaffold family (HG3 / 1A53 / KABLE) and pooled, so scaffold identity cannot create a correlation
  * SELECTION-AWARE permutation test: the null distribution of the BEST |rho| over the whole panel (labels permuted), i.e. 'how large a best rho
    does pure noise give when you screen this many metrics?'  Quote this, not the raw p of the winner.
  * leave-one-out composite (metric selection inside the loop)
PART B  arms
  * arm-resolved metrics table (TSA / substrate / product) with paired Wilcoxon tests across proteins
  * does the DIFFERENCE between arms (TSA-sub, TSA-prod, sub-prod) track kcat/KM?
  * pose differences between arms (arm_poses.csv) and apo-vs-holo (apo_holo.csv) vs kcat/KM
"""
from __future__ import annotations
import argparse, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata, wilcoxon
warnings.filterwarnings("ignore")
from .validate_dataset import screen, loo_composite, bh
from .validate import EXPECT

FAMILY = lambda n: "HG3" if n.startswith("HG3") else ("1A53" if n.startswith("1A53") else ("KABLE" if n.startswith("KABLE") else "other"))
IPCOLS = ["ipsae.p2l", "ipsae.l2p", "ipsae.iPAE12", "boltz.iptm", "boltz.ligand_iptm", "boltz.complex_iplddt", "site.pae_cat_to_lig", "site.pae_lig_to_cat", "site.plddt_cat", "site.plddt_lig"]

# arm-resolved metrics: label -> {arm: column in the master table}
ARM_METRICS = {
    "ipSAE2 protein->ligand": {a: f"ip.{a}.ipsae.p2l" for a in ("tsa", "sub", "prod")},
    "ligand ipTM": {a: f"ip.{a}.boltz.ligand_iptm" for a in ("tsa", "sub", "prod")},
    "PAE catalytic residues->ligand (A)": {a: f"ip.{a}.site.pae_cat_to_lig" for a in ("tsa", "sub", "prod")},
    "pLDDT ligand": {a: f"ip.{a}.site.plddt_lig" for a in ("tsa", "sub", "prod")},
    "PLACER holo prmsd (A)": {"tsa": "placer.holo.prmsd_mean", "sub": "placer_sub.holo.prmsd_mean", "prod": "placer_prod.holo.prmsd_mean"},
    "PLACER holo ligand RMSD to input (A)": {"tsa": "placer.holo.lig_rmsd_to_input", "sub": "placer_sub.holo.lig_rmsd_to_input", "prod": "placer_prod.holo.lig_rmsd_to_input"},
    "PLACER holo pose spread (A)": {"tsa": "placer.holo.lig_rmsd_spread", "sub": "placer_sub.holo.lig_rmsd_spread", "prod": "placer_prod.holo.lig_rmsd_spread"},
    "tmol interaction E (REU)": {"tsa": "rlx.int.total", "sub": "native_sub.int.total", "prod": "native_prod.int.total"},
    "tmol H-bond E (REU)": {"tsa": "rlx.int.hbond", "sub": "native_sub.int.hbond", "prod": "native_prod.int.hbond"},
    "tmol electrostatics (REU)": {"tsa": "rlx.int.fa_elec", "sub": "native_sub.int.fa_elec", "prod": "native_prod.int.fa_elec"},
    "tmol repulsion lig-prot (REU)": {"tsa": "rlx.int.fa_ljrep", "sub": "native_sub.int.fa_ljrep", "prod": "native_prod.int.fa_ljrep"},
    "tmol ligand burial (rel. SASA)": {"tsa": "rlx.geo.lig_rel_sasa", "sub": "native_sub.geo.lig_rel_sasa", "prod": "native_prod.geo.lig_rel_sasa"},
    "Rosetta dG_separated (REU)": {"tsa": "ros.tsa.dG_sep", "sub": "ros.native_sub.dG_sep", "prod": "ros.native_prod.dG_sep"},
    "Rosetta interface packstat": {"tsa": "ros.tsa.packstat_int", "sub": "ros.native_sub.packstat_int", "prod": "ros.native_prod.packstat_int"},
    "Rosetta H-bond E (REU)": {"tsa": "ros.tsa.hbond_E", "sub": "ros.native_sub.hbond_E", "prod": "ros.native_prod.hbond_E"},
}


def load_master(report: Path, work: Path, labels_csv: str) -> pd.DataFrame:
    ent = pd.read_csv(report / "summary_entries.csv", index_col=0)
    ent = ent[ent.index.str.endswith("_boltz")].copy()
    ent.index = ent.index.str.replace("_boltz", "", regex=False)
    ent = ent.drop(columns=[c for c in ent.columns if c.endswith("|sd")], errors="ignore")
    ipf = work / "ipsae" / "ipsae_models.csv"
    if ipf.exists():
        ip = pd.read_csv(ipf); ip = ip[ip.status == "ok"].copy(); ip["protein"] = ip.entry.str.replace("_boltz", "", regex=False)
        for arm, g in ip.groupby("arm"):
            m = g.groupby("protein")[[c for c in IPCOLS if c in g]].mean()
            m.columns = [f"ip.{arm}.{c}" for c in m.columns]
            ent = ent.join(m)
    for fn in ("apo_holo.csv", "arm_poses.csv"):
        f = work / fn
        if f.exists():
            d = pd.read_csv(f); d["protein"] = d.entry.str.replace("_boltz", "", regex=False)
            num = d.drop(columns=["rep_id", "entry"], errors="ignore").select_dtypes("number")
            ent = ent.join(num.groupby(d.protein).mean(), rsuffix="_dup")
    # labels (override whatever the aggregate attached)
    lab = pd.read_csv(labels_csv) if Path(labels_csv).exists() else pd.DataFrame(columns=["entry", "kcat_KM"])
    kc = {r.entry.lower(): r.kcat_KM for r in lab.itertuples() if pd.notna(r.kcat_KM)}
    ent["kcat_KM"] = [kc.get(p.lower(), np.nan) for p in ent.index]
    ent["log10_kcat_KM"] = np.log10(ent.kcat_KM.astype(float))
    ent["family"] = [FAMILY(p) for p in ent.index]
    # arm differences
    for label, cols in ARM_METRICS.items():
        key = label.split(" (")[0].replace(" ", "_").replace("->", "2").replace("/", "_")
        t, s, p = (ent.get(cols[a]) for a in ("tsa", "sub", "prod"))
        if t is not None and s is not None:
            ent[f"diff.{key}.tsa_minus_sub"] = t - s
        if t is not None and p is not None:
            ent[f"diff.{key}.tsa_minus_prod"] = t - p
        if s is not None and p is not None:
            ent[f"diff.{key}.sub_minus_prod"] = s - p
    return ent


def panel_columns(ent: pd.DataFrame) -> list[str]:
    keep = set(EXPECT)
    keep |= {c for c in ent.columns if c.startswith(("ip.", "ah.", "arm.", "diff.", "delta_native.", "ros.delta_native."))}
    keep |= {c for c in ent.columns if c.startswith(("placer_sub.", "placer_prod.", "placer_apoarm.")) and c.endswith(("frac_pass_core", "prmsd_mean", "lig_rmsd_to_input", "lig_rmsd_spread"))}
    keep |= {c for c in ent.columns if c.startswith("placer.") and c.endswith(("frac_pass_core", "prmsd_mean", "lig_rmsd_to_input", "lig_rmsd_spread", "chi_sd_base", "frac_donor"))}
    keep |= {c for c in ent.columns if c.startswith(("native_sub.int.", "native_prod.int.", "rlx.int.", "ros.native_", "ros.tsa.", "ros.sub.", "ros.prod.")) and c.count(".") <= 3}
    return [c for c in ent.columns if c in keep and ent[c].dtype.kind in "fi"]


def within_family_rho(x: np.ndarray, y: np.ndarray, fam: np.ndarray) -> float:
    rx, ry = np.zeros(len(x)), np.zeros(len(y))
    for f in np.unique(fam):
        i = fam == f
        if i.sum() < 2:
            continue
        rx[i] = rankdata(x[i]) - rankdata(x[i]).mean(); ry[i] = rankdata(y[i]) - rankdata(y[i]).mean()
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def best_rho_null(X: pd.DataFrame, y: np.ndarray, nperm: int = 10000, seed: int = 0) -> tuple[float, float, np.ndarray]:
    """Null distribution of max_m |rho_m| when y is permuted (complete-case metrics only). Returns (observed max, p, null)."""
    Xc = X.dropna(axis=1)
    if Xc.shape[1] == 0:
        return float("nan"), float("nan"), np.array([])
    Rx = np.apply_along_axis(rankdata, 0, Xc.to_numpy(float)); Rx = Rx - Rx.mean(0)
    Rx = Rx / np.maximum(np.linalg.norm(Rx, axis=0), 1e-12)
    ry = rankdata(y); rng = np.random.default_rng(seed)
    perms = np.stack([rng.permutation(ry) for _ in range(nperm)]); perms = perms - perms.mean(1, keepdims=True)
    perms = perms / np.linalg.norm(perms, axis=1, keepdims=True)
    null = np.abs(perms @ Rx).max(1)
    obs = float(np.abs((ry - ry.mean()) / np.linalg.norm(ry - ry.mean()) @ Rx).max())
    return obs, float((null >= obs).mean()), null


def arm_table(ent: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, cols in ARM_METRICS.items():
        d = {a: ent[cols[a]] for a in cols if cols[a] in ent}
        if len(d) < 2:
            continue
        row = {"metric": label, **{f"mean_{a}": d[a].mean() for a in d}, **{f"n_{a}": int(d[a].notna().sum()) for a in d}}
        for a, b in (("tsa", "sub"), ("tsa", "prod"), ("sub", "prod")):
            if a in d and b in d:
                x = pd.concat([d[a], d[b]], axis=1).dropna()
                if len(x) >= 6 and (x.iloc[:, 0] - x.iloc[:, 1]).abs().sum() > 0:
                    row[f"{a}-{b} median diff"] = float((x.iloc[:, 0] - x.iloc[:, 1]).median())
                    row[f"{a}-{b} wilcoxon p"] = float(wilcoxon(x.iloc[:, 0], x.iloc[:, 1]).pvalue)
        rows.append(row)
    return pd.DataFrame(rows)


def figures(ent, scr, arm_t, out: Path):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    col = {"HG3": "#2a6fdb", "1A53": "#e08a00", "KABLE": "#2aa05a", "other": "#888"}
    lab = ent[ent.log10_kcat_KM.notna()]
    top = scr.head(9).metric.tolist()
    if top:
        fig, ax = plt.subplots(3, 3, figsize=(12, 10))
        for a, m in zip(ax.ravel(), top):
            for f, g in lab.groupby("family"):
                a.scatter(g[m], g.log10_kcat_KM, c=col[f], s=40, label=f)
                for n, r in g.iterrows():
                    a.annotate(n, (r[m], r.log10_kcat_KM), fontsize=6, xytext=(3, 3), textcoords="offset points")
            rr = scr.set_index("metric").loc[m]
            a.set_title(f"{m[:38]}\nrho={rr.rho:+.2f}  q={rr.q_BH:.2f}", fontsize=8); a.set_ylabel("log10 kcat/KM", fontsize=7)
        ax[0, 0].legend(fontsize=7); fig.tight_layout(); fig.savefig(out / "fig_top_correlations.png", dpi=130); plt.close(fig)
    metrics = [l for l in ARM_METRICS if l in set(arm_t.metric)][:12]
    if metrics:
        n = len(metrics); fig, ax = plt.subplots((n + 3) // 4, 4, figsize=(14, 3.2 * ((n + 3) // 4)))
        for a, l in zip(ax.ravel(), metrics):
            cols = ARM_METRICS[l]
            for p, r in ent.iterrows():
                ys = [r.get(cols.get(arm), np.nan) for arm in ("tsa", "sub", "prod")]
                a.plot([0, 1, 2], ys, "-o", c=col[r.family], ms=3, lw=0.8, alpha=0.8)
            a.set_xticks([0, 1, 2]); a.set_xticklabels(["TSA\n6NT", "substrate\nH5J", "product\n6VP"], fontsize=7); a.set_title(l, fontsize=8)
        fig.tight_layout(); fig.savefig(out / "fig_arm_comparison.png", dpi=130); plt.close(fig)
    if len(scr):
        s = scr.head(30).iloc[::-1]
        fig, a = plt.subplots(figsize=(8, 9))
        a.barh(s.metric.str[:48], s.rho, color=["#c0392b" if q < 0.1 else "#999" for q in s.q_BH]); a.axvline(0, c="k", lw=0.5)
        a.set_xlabel("Spearman rho vs log10 kcat/KM  (red = BH q < 0.1)"); fig.tight_layout(); fig.savefig(out / "fig_correlation_bars.png", dpi=130); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", default="work/report"); ap.add_argument("--work", default="work")
    ap.add_argument("--labels", default="config/labels.csv"); ap.add_argument("--out", default="work/report/analysis")
    ap.add_argument("--active-threshold", type=float, default=10000.0)
    ap.add_argument("--control", default="rlx.geo.n_donors_3p5")
    ap.add_argument("--family", default=None, help="restrict the correlation analysis to one scaffold family (HG3 / 1A53 / KABLE)")
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    ent = load_master(Path(a.report), Path(a.work), a.labels)
    ent.to_csv(out / "master_table.csv")
    lab = ent[ent.log10_kcat_KM.notna()]
    if a.family:
        lab = lab[lab.family == a.family]
    y = lab.log10_kcat_KM.to_numpy(float); fam = lab.family.to_numpy()
    cols = panel_columns(ent)
    X = lab[cols]
    md = [f"# Full analysis: {len(ent)} proteins, {len(lab)} labelled ({', '.join(f'{k}: {v}' for k, v in lab.family.value_counts().items())})\n",
          f"Panel = {len(cols)} metrics declared in advance (see panel_columns); n={len(lab)}.  Active = kcat/KM >= {a.active_threshold:g}.\n"]
    # ---------------- PART A
    active = (lab.kcat_KM.to_numpy(float) >= a.active_threshold).astype(float)
    ctrl = X[a.control].to_numpy(float) if a.control in X else None
    scr = screen(X, y, active, ctrl)
    if len(scr):
        wf = []
        for m in scr.metric:
            v = X[m].to_numpy(float); ok = np.isfinite(v)
            wf.append(within_family_rho(v[ok], y[ok], fam[ok]) if ok.sum() >= 6 else np.nan)
        scr["rho_within_family"] = wf
        scr.to_csv(out / "correlations_panel.csv", index=False)
        obs, pnull, null = best_rho_null(X, y)
        md += ["## A1. Correlation with log10 kcat/KM (panel)\n",
               f"Metrics with BH q < 0.1: **{int((scr.q_BH < 0.1).sum())}** of {len(scr)}.\n",
               f"**Selection-aware test:** best |rho| over the complete-case panel = **{obs:.2f}**; permuting the labels gives a best |rho| >= that in **{pnull:.1%}** of "
               f"10,000 permutations (95th percentile of the null best |rho|: {np.percentile(null, 95):.2f}). If this p is not small, the top row below is what noise produces when many metrics are screened.\n",
               "```", scr.head(25).round(3).to_string(index=False), "```\n",
               "`rho_within_family` = rank correlation after ranking inside each scaffold family (HG3 / 1A53 / KABLE), i.e. with scaffold effects removed; `partial_rho` controls for oxyanion-donor presence.\n"]
        r_loo, _ = loo_composite(X.dropna(axis=1, thresh=int(0.8 * len(X))), y, 3)
        md += [f"## A2. Leave-one-out composite (top-3 metrics chosen inside each fold)\nheld-out Spearman = **{r_loo:.2f}** (n={len(y)}).\n"]
    # ---------------- PART B
    arm_t = arm_table(ent)
    arm_t.to_csv(out / "arm_comparison.csv", index=False)
    md += ["## B1. TSA vs substrate vs product (all proteins; paired Wilcoxon across proteins)\n", "```", arm_t.round(3).to_string(index=False), "```\n"]
    dcols = [c for c in ent.columns if c.startswith("diff.") and ent.loc[lab.index, c].notna().sum() >= 8]
    if dcols:
        sd = screen(lab[dcols], y, active, None)
        if len(sd):
            sd["rho_within_family"] = [within_family_rho(lab[m].to_numpy(float)[np.isfinite(lab[m])], y[np.isfinite(lab[m])], fam[np.isfinite(lab[m])]) for m in sd.metric]
            sd.to_csv(out / "arm_difference_vs_kcat.csv", index=False)
            md += ["## B2. Does the DIFFERENCE between arms track kcat/KM?\n", "```", sd.head(15).round(3).to_string(index=False), "```\n"]
    pose = [c for c in ent.columns if c.startswith(("arm.", "ah."))]
    if pose:
        pt = ent[pose].describe().T[["count", "mean", "50%", "min", "max"]]
        md += ["## B3. Pose differences between arms (A) and apo-vs-holo, all proteins\n", "```", pt.round(2).to_string(), "```\n"]
    figures(ent, scr if len(scr) else pd.DataFrame(columns=["metric", "rho", "q_BH"]), arm_t, out)
    (out / "analysis.md").write_text("\n".join(md) + "\n")
    print("\n".join(md[:8])); print("->", out)

if __name__ == "__main__":
    main()
