"""Steering objectives: every generated design is scored with the SAME reward terms the beam search used (weight x min(value, clip) / norm; TMol clipped to [-80, 20]); unsteered designs included for comparison.
Reward A = steered arm (13 terms, ligand-SASA weight -0.5); reward B = channel arm (ligand-SASA weight -0.25 + bottleneck band term).  Scores are computed on the FINAL designs (the search also scored coarser structures at steps 30/42/52).
Output: ../steering_objectives/20261122_eNN_<name>.png"""
import glob, textwrap
from kstyle import *
D = pd.read_csv(os.path.join(T, "all_designs_metrics.csv"), low_memory=False); D = D[~D.arm.str.contains("sheet")]; D = D[D.arm.isin(ARMS)].copy()
OUTD = os.path.join(OUT, "steering_objectives"); os.makedirs(OUTD, exist_ok=True)
# term: (column, weight, clip_min, clip_max, norm, label)
T_A = [("base_sc_polar_4p5", -4, None, 3, 1, "polar atoms near base"), ("base_cation_6", -3, None, 2, 1, "Lys/Arg near base"), ("base_sasa_rel", -2, None, 0.5, 0.5, "base SASA"), ("base_apolar_frac", 0.5, None, 1, 1, "apolar fraction at base"),
       ("volume_ratio_band_error", -1, None, 0.6, 0.3, "pocket volume band"), ("ligand_sasa_frac", -0.5, None, 0.2, 0.2, "ligand exposure"), ("nitro_hb_near_base", -1, None, 2, 2, "nitro donor near base"), ("nitro_hb_far", 0.5, None, 1, 1, "nitro donor far from base"),
       ("lig_clash_atoms", -1, None, 3, 3, "ligand clashes"), ("pp_clash_frac", -1, None, 0.01, 0.01, "protein clashes"), ("motif_base_rmsd", -2, None, 1, 0.25, "Asp rigidity (RMSD)"), ("motif_partner_rmsd", -1.5, None, 1, 0.25, "Gln rigidity (RMSD)"),
       ("tmol_ddg", -0.5, -80, 20, 20, "TMol interaction energy")]
T_B = [t if t[0] != "ligand_sasa_frac" else ("ligand_sasa_frac", -0.25, None, 0.2, 0.2, "ligand exposure") for t in T_A] + [("channel_band_error", -2, None, 1, 0.3, "channel bottleneck band")]
def contrib(df, T): return pd.DataFrame({t[5]: t[1] * df[t[0]].clip(lower=t[2], upper=t[3]) / t[4] for t in T}).fillna(0.0)
for nm, T in (("A", T_A), ("B", T_B)):
    C = contrib(D, T); D["reward_" + nm] = C.sum(1)
    for c in C: D[f"{nm}:{c}"] = C[c]
NSRC = f"Source: 20261122_kemp_c1, {len(D)} generated designs, scored on the final structures with the steering reward terms (weight x min(value, clip) / norm)."
def sv(f, k, name): f.savefig(os.path.join(OUTD, f"{DATE}_e{k:02d}_{name}.png"), bbox_inches="tight", pad_inches=0.05); plt.close(f); print("saved", k, name)
def legend_arms(ax):
    hs = [plt.Line2D([], [], color=COL[a], lw=1.5, ls="--" if a == "unsteered" else "-") for a in ARMS]; legend_mono(ax, hs, [f"{LAB[a]:<22} n={int((D.arm == a).sum()):>4}" for a in ARMS], title="designs (both motifs)")
def hist2(k, name, col, xlabel, lo, hi, nb=30, foot=""):
    f, axs = fig(183, 62, ncols=2, sharey=True); f.subplots_adjust(wspace=0.06, right=0.72); bins = np.linspace(lo, hi, nb + 1)
    for ax, m in zip(axs, ("hgnt", "hgwt")):
        for a in ARMS:
            v = D[(D.motif == m) & (D.arm == a)][col].dropna().clip(lo, hi - 1e-9); w = np.ones(len(v)) / max(len(v), 1)
            ax.hist(v, bins=bins, weights=w, histtype="stepfilled", color=COL[a], alpha=0.16, lw=0); ax.hist(v, bins=bins, weights=w, histtype="step", color=COL[a], lw=1.4, ls="--" if a == "unsteered" else "-"); ax.plot([v.median()], [0], marker="v", ms=4.5, color=COL[a], mec="white", mew=0.5, clip_on=False, zorder=6)
        ax.set_xlim(lo, hi); ax.set_xlabel(textwrap.fill(xlabel, 34)); ax.set_title(MOTLAB[m], fontsize=7, loc="left")
    axs[0].set_ylabel("fraction of designs"); legend_arms(axs[1]); footnote(f, NSRC + " Triangles = medians." + foot, y=-0.15); sv(f, k, name)
hist2(0, "total_reward_steered_function", "reward_A", "total steering reward, steered-arm terms (higher is better)", -32, 4, foot=" Maximum possible = +2.5 (13 terms).")
hist2(1, "total_reward_channel_function", "reward_B", "total steering reward incl. channel term (higher is better)", -38, 4, foot=" Channel-arm terms: ligand exposure weight -0.25 plus bottleneck band.")
# mean contribution per term and arm
def breakdown(k, name, T, pref, title_note):
    labs = [t[5] for t in T]; f, axs = fig(183, 92, ncols=2, sharey=True); f.subplots_adjust(wspace=0.08, right=0.72); w = 0.26; y = np.arange(len(labs))
    for ax, m in zip(axs, ("hgnt", "hgwt")):
        for i, a in enumerate(ARMS):
            d = D[(D.motif == m) & (D.arm == a)]; v = [d[f"{pref}:{l}"].mean() for l in labs]
            ax.barh(y + (1 - i) * w, v, w * 0.92, color="white" if a == "unsteered" else COL[a], edgecolor=COL[a], linewidth=1.0 if a == "unsteered" else 0)
        ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=6.8); ax.invert_yaxis(); ax.axvline(0, color=INK, lw=0.6); ax.set_xlabel("mean reward contribution per design"); ax.set_title(MOTLAB[m], fontsize=7, loc="left")
    legend_arms(axs[1]); footnote(f, NSRC + title_note, y=-0.08); sv(f, k, name)
breakdown(2, "reward_breakdown_steered_function", T_A, "A", " Bars to the right of zero are rewards, to the left penalties; mean over designs.")
breakdown(3, "reward_breakdown_channel_function", T_B, "B", " Includes the channel bottleneck band term.")
# gain per term vs unsteered
def gain(k, name, T, pref):
    labs = [t[5] for t in T]; f, axs = fig(183, 92, ncols=2, sharey=True); f.subplots_adjust(wspace=0.08, right=0.72); y = np.arange(len(labs)); w = 0.36
    for ax, m in zip(axs, ("hgnt", "hgwt")):
        base = D[(D.motif == m) & (D.arm == "unsteered")]
        for i, a in enumerate(("steer_tmol", "steer_tmol_ch18")):
            d = D[(D.motif == m) & (D.arm == a)]; ax.barh(y + (i - 0.5) * w, [d[f"{pref}:{l}"].mean() - base[f"{pref}:{l}"].mean() for l in labs], w * 0.92, color=COL[a])
        ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=6.8); ax.invert_yaxis(); ax.axvline(0, color=INK, lw=0.6); ax.set_xlabel("reward gain over unsteered (per design)"); ax.set_title(MOTLAB[m], fontsize=7, loc="left")
    hs = [plt.Rectangle((0, 0), 1, 1, fc=COL[a]) for a in ("steer_tmol", "steer_tmol_ch18")]; legend_mono(axs[1], hs, [LAB[a] for a in ("steer_tmol", "steer_tmol_ch18")], title="arm minus unsteered")
    footnote(f, NSRC + " Positive = the arm scores better on that term than the unsteered designs.", y=-0.08); sv(f, k, name)
gain(4, "reward_gain_over_unsteered", T_B, "B")
# raw plots for terms not yet shown
def hist_raw(k, name, col, xlabel, lo, hi, nb=26, vline=None, foot=""):
    hist2(k, name, col, xlabel, lo, hi, nb, foot=foot)
hist_raw(5, "pocket_volume_band_error", "volume_ratio_band_error", "distance of pocket ratio outside 0.5 +- 0.1", 0, 1.2, foot=" 0 = inside the steering band.")
hist_raw(6, "protein_clash_fraction", "pp_clash_frac", "protein-protein clashes per atom (pairs < 2.2 A)", 0, 0.012, nb=24, foot=" Gate: <= 0.005.")
hist_raw(7, "channel_band_error", "channel_band_error", "bottleneck radius outside 1.7-1.9 A (A)", 0, 1.3, nb=26, foot=" 0 = bottleneck inside the target band.")
