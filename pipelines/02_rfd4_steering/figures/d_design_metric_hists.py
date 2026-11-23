"""Design-level comparison of unsteered vs steered vs steered+channel designs, one PNG per design-gate metric (two panels: the two motifs, same axes).
Data: every generated design of the finished runs (tables/all_designs_metrics.csv; sheet runs excluded).  Output: ../design_metrics/20261122_dNN_<metric>.png"""
import glob, textwrap
from kstyle import *
D = pd.read_csv(os.path.join(T, "all_designs_metrics.csv"), low_memory=False); D = D[~D.arm.str.contains("sheet")]; D = D[D.arm.isin(ARMS)]
OUTD = os.path.join(OUT, "design_metrics"); os.makedirs(OUTD, exist_ok=True)
NOTE_SRC = f"Source: 20261122_kemp_c1, {len(D)} generated designs (non-sheet runs); all designs are shown before any gate (unsteered runs are gated on geometry, steered on the full core filter)."
def savep(f, k, name): f.savefig(os.path.join(OUTD, f"{DATE}_d{k:02d}_{name}.png"), bbox_inches="tight", pad_inches=0.05); plt.close(f); print("saved", k, name)
def legend_arms(ax):
    hs = [plt.Line2D([], [], color=COL[a], lw=1.5, ls="--" if a == "unsteered" else "-") for a in ARMS]
    lab = [f"{LAB[a]:<22} n={int((D.arm == a).sum()):>4}" for a in ARMS]; legend_mono(ax, hs, lab, title="designs (both motifs)")
def cont(k, name, col, xlabel, lo, hi, nb=26, thr=None, band=None, target=None, clip=True, thr_label=None):
    f, axs = fig(183, 62, ncols=2, sharey=True); f.subplots_adjust(wspace=0.06, right=0.72); bins = np.linspace(lo, hi, nb + 1)
    for ax, m in zip(axs, ("hgnt", "hgwt")):
        if band: ax.axvspan(*band, color="#52606d", alpha=0.10, lw=0)
        if target: ax.axvspan(*target, color="#e2661a", alpha=0.10, lw=0)
        for a in ARMS:
            v = D[(D.motif == m) & (D.arm == a)][col].dropna(); v = v.clip(lo, hi - 1e-9) if clip else v
            if len(v) == 0: continue
            w = np.ones(len(v)) / len(v)
            ax.hist(v, bins=bins, weights=w, histtype="stepfilled", color=COL[a], alpha=0.16, lw=0); ax.hist(v, bins=bins, weights=w, histtype="step", color=COL[a], lw=1.4, ls="--" if a == "unsteered" else "-")
            ax.plot([v.median()], [0], marker="v", ms=4.5, color=COL[a], mec="white", mew=0.5, clip_on=False, zorder=6)
        if thr is not None: ax.axvline(thr, color=INK, lw=0.8, ls=":")
        ax.set_xlim(lo, hi); ax.set_xlabel(textwrap.fill(xlabel, 34)); ax.set_title(MOTLAB[m], fontsize=7, loc="left")
    axs[0].set_ylabel("fraction of designs")
    legend_arms(axs[1]); extra = []
    if thr is not None: extra.append(f"dotted line = gate threshold{thr_label or ''}")
    if band: extra.append("grey band = gate range")
    if target: extra.append("orange band = steering target")
    footnote(f, NOTE_SRC + (" Triangles = medians; " + "; ".join(extra) + "." if extra else " Triangles = medians."), y=-0.13); savep(f, k, name)
def disc(k, name, col, xlabel, cats, thr_note):
    f, axs = fig(183, 62, ncols=2, sharey=True); f.subplots_adjust(wspace=0.06, right=0.72); w = 0.26
    for ax, m in zip(axs, ("hgnt", "hgwt")):
        for i, a in enumerate(ARMS):
            v = D[(D.motif == m) & (D.arm == a)][col].dropna(); n = len(v)
            fr = [((v >= c0) & (v <= c1)).sum() / n if n else 0 for c0, c1, _ in cats]
            ax.bar(np.arange(len(cats)) + (i - 1) * w, fr, w * 0.92, color="white" if a == "unsteered" else COL[a], edgecolor=COL[a], linewidth=1.0 if a == "unsteered" else 0)
        ax.set_xticks(np.arange(len(cats))); ax.set_xticklabels([c[2] for c in cats]); ax.set_xlabel(textwrap.fill(xlabel, 30)); ax.set_title(MOTLAB[m], fontsize=7, loc="left")
    axs[0].set_ylabel("fraction of designs"); legend_arms(axs[1]); footnote(f, NOTE_SRC + " " + thr_note, y=-0.2); savep(f, k, name)
cont(1, "ligand_sasa_fraction", "ligand_sasa_frac", "ligand SASA fraction (exposed fraction of 6NT)", 0, 0.5, thr=0.10, thr_label=" (<= 0.10)")
cont(2, "base_relative_sasa", "base_sasa_rel", "base carboxylate relative SASA", 0, 0.6, thr=0.10, thr_label=" (<= 0.10)")
cont(3, "pocket_volume_ratio", "volume_ratio", "pocket volume / ligand volume", 0, 2.4, band=(0.4, 1.4), target=(0.4, 0.6))
disc(4, "polar_atoms_near_base", "base_sc_polar_4p5", "side-chain polar atoms within 4.5 A of the base O", [(0, 0, "0"), (1, 1, "1"), (2, 2, "2"), (3, 99, ">=3")], "Gate: 0 (hgwt: <= 1, Trp44 NE1).")
disc(5, "cations_near_base", "base_cation_6", "Lys/Arg within 6 A of the base O", [(0, 0, "0"), (1, 1, "1"), (2, 99, ">=2")], "Gate: 0.")
disc(6, "nitro_donor_far_from_base", "nitro_hb_far", "nitro O atoms with a donor >= 6 A from the base", [(0, 0, "0"), (1, 1, "1"), (2, 2, "2")], "Gate: >= 1 (in the gate together with: none with near-base donors only).")
disc(7, "nitro_donor_near_base", "nitro_hb_near_base", "nitro O atoms with only near-base donors", [(0, 0, "0"), (1, 1, "1"), (2, 2, "2")], "Gate: 0.")
disc(8, "ligand_clashes", "lig_clash_atoms", "protein-ligand atom pairs < 2.4 A", [(0, 0, "0"), (1, 2, "1-2"), (3, 6, "3-6"), (7, 999, ">=7")], "Gate: 0.")
cont(9, "motif_base_rmsd", "motif_base_rmsd", "catalytic Asp tip RMSD to theozyme (A)", 0, 1.0, thr=0.3, thr_label=" (<= 0.3 A)")
cont(10, "motif_donor_rmsd", "motif_partner_rmsd", "oxyanion Gln tip RMSD to theozyme (A)", 0, 1.0, thr=0.3, thr_label=" (<= 0.3 A)")
cont(11, "tmol_interaction_energy", "tmol_ddg", "TMol protein-6NT interaction energy (REU)", -35, 10, thr=0.0, thr_label=" (< 0)")
cont(12, "base_apolar_fraction", "base_apolar_frac", "apolar fraction within 6 A of the base O", 0.3, 1.0)
cont(13, "channel_bottleneck_radius", "bottleneck_radius", "channel bottleneck radius (A)", 0.4, 2.4, nb=20, target=(1.7, 1.9))
# pass rates per filter
crit = [("ok_lig_sasa", "ligand SASA"), ("ok_base_buried", "base SASA"), ("ok_pocket", "pocket ratio"), ("ok_clash", "no clashes"), ("ok_motif", "motif RMSD"), ("ok_no_sc_polar", "no polar near base"), ("ok_no_cation", "no cation"), ("ok_nitro", "nitro donor"), ("ok_tmol", "TMol < 0"), ("pass_geom", "all geometry"), ("pass_core", "all (core)")]
f, axs = fig(183, 70, ncols=2, sharey=True); f.subplots_adjust(wspace=0.06, right=0.72); w = 0.26
for ax, m in zip(axs, ("hgnt", "hgwt")):
    for i, a in enumerate(ARMS):
        d = D[(D.motif == m) & (D.arm == a)]; ax.bar(np.arange(len(crit)) + (i - 1) * w, [100 * d[c].mean() for c, _ in crit], w * 0.92, color="white" if a == "unsteered" else COL[a], edgecolor=COL[a], linewidth=1.0 if a == "unsteered" else 0)
    ax.set_xticks(np.arange(len(crit))); ax.set_xticklabels([l for _, l in crit], rotation=40, ha="right", fontsize=6.5); ax.set_title(MOTLAB[m], fontsize=7, loc="left"); ax.set_ylim(0, 105)
axs[0].set_ylabel("designs passing (%)"); legend_arms(axs[1]); footnote(f, NOTE_SRC + " Each bar = share of designs passing that criterion alone; the last two bars combine criteria.", y=-0.2); savep(f, 0, "filter_pass_rates")
