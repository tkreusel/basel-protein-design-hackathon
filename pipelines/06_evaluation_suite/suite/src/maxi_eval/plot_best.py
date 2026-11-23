"""Publication-style scatter plots of the best correlations with activity.

  python -m maxi_eval.plot_best [--screen work/report/everything] [--out work/report/figures]

Figures:
  fig1_best_correlations.png      the 9 best candidates (all 19 proteins), one theme per panel
  fig2_preorganisation.png        the preorganisation-type metrics (apo pre-formed, apo/holo motion, base burial, substrate NAC, constraint penalty, ...)
  fig3_signal_by_feature_set.png  best |rho| per feature set vs the selection-aware noise level (95th percentile of the null)
  fig4_HG3_only.png               the best candidates inside the HG3 family (n=12)
Each scatter: points coloured by scaffold family and labelled, least-squares line with bootstrap 95% band, Spearman rho / BH q / leave-one-protein-out range /
within-family rho in the corner.  Colour-blind-safe palette.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd

COL = {"HG3": "#0072B2", "1A53": "#E69F00", "KABLE": "#009E73"}
NICE = {
    "tmol swap|E_total|tsa-sub": ("tmol: total-energy difference TSA − substrate (REU)", "more favourable for the TSA →"),
    "tmol swap|E_total|tsa-prod": ("tmol: total-energy difference TSA − product (REU)", ""),
    "pr.tsa_pk-sub.int_polar": ("tmol: ligand–polar-residue interaction, TSA − substrate (REU)", ""),
    "pr.tsa_pk-sub.pocket_E_sum": ("tmol: pocket-residue energy sum, TSA − substrate (REU)", ""),
    "bz.sub-prod.cat_to_pocket_pae": ("Boltz: catalytic→pocket PAE, substrate − product (Å)", ""),
    "ipSAE2/Boltz|boltz.ligand_iptm|sub-prod": ("Boltz: ligand ipTM, substrate − product", ""),
    "tmol native|int.total|tsa": ("tmol: TSA interaction energy (REU, independent refold)", ""),
    "placer_sub.holo.lig_rmsd_to_input": ("PLACER: substrate ligand RMSD to Boltz pose (Å)", ""),
    "pl2.sub_holo.nac30": ("PLACER: substrate near-attack fraction (d≤3 Å, ∠≥130°)", ""),
    "placer.apo.frac_pass_core": ("PLACER apo: fraction pre-formed (HG3 windows)", ""),
    "placer_apoarm.apo.frac_pass_core": ("PLACER on Boltz-apo: fraction pre-formed", ""),
    "ah.apo.cst_n_pass": ("TSA transplanted into Boltz-apo pocket: constraints met", ""),
    "ah.pocket_bb_rmsd": ("Boltz apo vs holo: pocket backbone RMSD (Å)", ""),
    "tsa_pk.geo.base_sasa": ("catalytic-base SASA with TSA bound (Å²)", ""),
    "sub.cst.penalty": ("theozyme-constraint penalty, substrate complex", ""),
    "strain.E_noncst": ("tmol: strain of imposing the theozyme (REU)", ""),
    "apo.reorg_E": ("tmol: apo reorganisation energy (REU)", ""),
    "rlx.cst.n_pass": ("theozyme constraints met after relax (of 10)", ""),
    "pr.tsa_pk-sub.int_base": ("tmol: base-residue interaction, TSA − substrate (REU)", ""),
    "tmol swap|int.lk_ball_iso|tsa-prod": ("tmol: ligand–pocket lk_ball_iso (isotropic solvation), TSA − product (REU)", ""),
    "tmol swap|int.total|tsa": ("tmol: TSA–pocket interaction energy (REU)", ""),
    "tmol swap|int.fa_lk|tsa-sub": ("tmol: ligand–pocket fa_lk (LK solvation), TSA − substrate (REU)", ""),
    "tmol swap|int.total|tsa-sub": ("tmol: ligand–pocket total interaction, TSA − substrate (REU)", ""),
    "pr.tsa_pk.int_sum": ("tmol: summed per-residue ligand interaction, TSA complex (REU)", ""),
    "tmol swap|int.fa_elec|sub-prod": ("tmol: ligand–pocket electrostatics, substrate − product (REU)", ""),
    "tmol swap|buns.n_unsat_ligand|tsa-prod": ("tmol: buried unsatisfied ligand polars, TSA − product", ""),
    "sub.geo.d_O_H": ("substrate complex: base O···H distance (Å)", ""),
    "placer_sub.holo.ang_C_O_H_mean": ("PLACER substrate: mean C–O···H angle (°)", ""),
    "pl2.tsa_holo.donor_med": ("PLACER TSA: median donor distance (Å)", ""),
    "rlx.geo.n_donors_3p5": ("oxyanion donors within 3.5 Å (TSA)", ""),
}


def fit_band(x, y, rng, n=400):
    xs = np.linspace(x.min(), x.max(), 60); P = []
    for _ in range(n):
        i = rng.integers(0, len(x), len(x))
        if np.std(x[i]) == 0:
            continue
        b = np.polyfit(x[i], y[i], 1); P.append(np.polyval(b, xs))
    P = np.array(P); b0 = np.polyfit(x, y, 1)
    return xs, np.polyval(b0, xs), np.percentile(P, 2.5, 0), np.percentile(P, 97.5, 0)


def panel(ax, name, M, lab, stats, rng):
    x = M[name].reindex(lab.index).to_numpy(float); y = lab.log10_kcat_KM.to_numpy(float); ok = np.isfinite(x)
    ttl, _ = NICE.get(name, (name, ""))
    xs, ys, lo, hi = fit_band(x[ok], y[ok], rng)
    ax.fill_between(xs, lo, hi, color="#bbbbbb", alpha=0.35, lw=0); ax.plot(xs, ys, color="#555555", lw=1.2)
    for f, c in COL.items():
        m = (lab.family.to_numpy() == f) & ok
        ax.scatter(x[m], y[m], s=58, c=c, edgecolor="white", lw=0.8, zorder=3, label=f"{f} (n={int(m.sum())})")
    for xi, yi, n in zip(x[ok], y[ok], lab.index[ok]):
        ax.annotate(n.replace("HG3.", "").replace("1A53-", "") if n.startswith(("HG3.", "1A53-")) else n, (xi, yi), fontsize=5.8, xytext=(3.5, 3), textcoords="offset points", color="#333333")
    s = stats.get(name)
    if s is not None:
        t = f"ρ = {s['rho']:+.2f}   q = {s['q_BH']:.2f}\nLOPO {s['rho_loo_min']:+.2f}…{s['rho_loo_max']:+.2f}"
        if np.isfinite(s.get("rho_within_family", np.nan)):
            t += f"\nwithin-family ρ = {s['rho_within_family']:+.2f}"
        if s.get("extra"):
            t += "\n" + s["extra"]
        ax.text(0.03 if s["rho"] < 0 else 0.97, 0.04, t, transform=ax.transAxes, fontsize=7, va="bottom", ha="left" if s["rho"] < 0 else "right",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", alpha=0.9))
    ax.set_title(ttl, fontsize=8.2, loc="left"); ax.set_xlabel("")
    ax.set_yticks(range(0, 8)); ax.set_yticklabels([f"10$^{{{k}}}$" for k in range(0, 8)], fontsize=7); ax.tick_params(axis="x", labelsize=7)
    ax.grid(True, color="#eeeeee", lw=0.7, zorder=0); [ax.spines[k].set_visible(False) for k in ("top", "right")]


def grid(names, M, lab, stats, out, title, ncols=3, sub=""):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rows = (len(names) + ncols - 1) // ncols
    fig, axs = plt.subplots(rows, ncols, figsize=(5.0 * ncols, 3.9 * rows), facecolor="white"); axs = np.atleast_1d(axs).ravel(); rng = np.random.default_rng(0)
    for a, n in zip(axs, names):
        panel(a, n, M, lab, stats, rng)
    for a in axs[len(names):]:
        a.axis("off")
    for a in axs[::ncols]:
        a.set_ylabel("kcat/KM (M$^{-1}$s$^{-1}$)", fontsize=8)
    h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="upper right", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(title, x=0.01, ha="left", fontsize=12, fontweight="bold"); fig.text(0.01, 0.945, sub, fontsize=8, color="#555555")
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(out, dpi=170); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", default="work/report/everything"); ap.add_argument("--screen-hg3", default="work/report/everything_HG3")
    ap.add_argument("--labels", default="config/labels_manual_dataset.csv"); ap.add_argument("--out", default="work/report/figures"); ap.add_argument("--best", default=None); ap.add_argument("--hg", default=None)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    from .analyze_full import load_master
    ent = load_master(Path("work/report"), Path("work"), a.labels); lab = ent[ent.log10_kcat_KM.notna()]
    M = pd.read_csv(Path(a.screen) / "candidate_matrix.csv", index_col=0)
    S = pd.read_csv(Path(a.screen) / "candidates_all.csv").set_index("candidate"); stats = S.to_dict("index")
    sub = f"Spearman ρ vs log10 kcat/KM, n=19 (manual dataset; tmol = mean of two independent relax runs). q = Benjamini–Hochberg over {len(S)} screened candidates; LOPO = range of ρ when any one protein is dropped. Shaded: bootstrap 95% band."
    best = a.best.split(";") if a.best else ["tmol swap|int.lk_ball_iso|tsa-prod", "tmol swap|int.fa_lk|tsa-sub", "tmol swap|int.total|tsa-sub", "bz.sub-prod.cat_to_pocket_pae",
            "ipSAE2/Boltz|boltz.ligand_iptm|sub-prod", "pr.tsa_pk.int_sum", "tmol swap|E_total|tsa-prod", "tmol native|int.total|tsa", "placer_sub.holo.lig_rmsd_to_input"]
    grid(best, M, lab, stats, out / "fig1_best_correlations.png", "Best correlations with enzyme activity", sub=sub)
    pre = ["placer.apo.frac_pass_core", "placer_apoarm.apo.frac_pass_core", "ah.apo.cst_n_pass", "ah.pocket_bb_rmsd", "tsa_pk.geo.base_sasa", "sub.cst.penalty",
           "pl2.sub_holo.nac30", "strain.E_noncst", "apo.reorg_E"]
    grid(pre, M, lab, stats, out / "fig2_preorganisation.png", "Preorganisation-type metrics vs activity", sub=sub)
    # fig3: best |rho| per set vs null
    N = pd.read_csv(Path(a.screen) / "selection_aware_by_set.csv"); N = N[~N.group.isin(["ALL"])].sort_values("best_abs_rho")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.5, 5.2), facecolor="white")
    y = np.arange(len(N)); ax.barh(y, N.best_abs_rho, color=["#c0392b" if p < 0.05 else "#9aa5b1" for p in N.p_selection_aware], height=0.62)
    ax.scatter(N.null95, y, marker="|", s=420, c="black", zorder=3, label="95th percentile of the noise (labels permuted)")
    for yi, (v, p) in enumerate(zip(N.best_abs_rho, N.p_selection_aware)):
        ax.text(v + 0.01, yi, f"p = {p:.3f}", va="center", fontsize=8)
    ax.set_yticks(y); ax.set_yticklabels([g.replace("set:", "") for g in N.group], fontsize=8); ax.set_xlim(0, 1)
    ax.set_xlabel("best |Spearman ρ| in the feature set"); ax.set_title("Where is the signal? Best correlation per feature set vs selection-aware noise level (red: p < 0.05)", fontsize=10, loc="left")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, -0.13), fontsize=8, frameon=False); [ax.spines[k].set_visible(False) for k in ("top", "right")]; fig.tight_layout(); fig.savefig(out / "fig3_signal_by_feature_set.png", dpi=170); plt.close(fig)
    # fig4 HG3 only
    MH = pd.read_csv(Path(a.screen_hg3) / "candidate_matrix.csv", index_col=0); SH = pd.read_csv(Path(a.screen_hg3) / "candidates_all.csv").set_index("candidate")
    labh = lab[lab.family == "HG3"]
    hg = a.hg.split(";") if a.hg else ["tmol swap|int.lk_ball_iso|tsa-prod", "tmol swap|int.fa_elec|sub-prod", "placer_sub.holo.ang_C_O_H_mean", "tmol swap|buns.n_unsat_ligand|tsa-prod", "sub.geo.d_O_H", "tsa_pk.geo.base_sasa"]
    NH = pd.read_csv(Path(a.screen_hg3) / "selection_aware_by_set.csv").set_index("group").loc["ALL"]
    grid(hg, MH, labh, SH.to_dict("index"), out / "fig4_HG3_only.png", "HG3 family only (n = 12): best candidates", ncols=3,
         sub=f"Within the HG3 scaffold (valid reference windows; tmol = mean of two runs). Selection-aware p over {len(SH)} candidates = {NH.p_selection_aware:.3f} (noise 95th pct {NH.null95:.2f}); n = 12, so treat as a lead to test on new proteins.")
    print("figures ->", out)

if __name__ == "__main__":
    main()
