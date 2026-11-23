"""f01: yield funnel per arm (counts per 300 generated designs, log y).  Data: ../../tables/{funnel,all_designs_metrics,candidates,evaluated_sequences_full_suite}.csv; only runs with stage B finished, sheet runs excluded."""
from kstyle import *
fu = pd.read_csv(f"{T}/funnel.csv"); fu = fu[(fu.stageB == True) & ~fu.arm.str.contains("sheet")]
D = designs(); c = cands(); e = ev()
runs = set(zip(fu.motif, fu.arm, fu.batch)); key = lambda df: [(m, a, b) in runs for m, a, b in zip(df.motif, df.arm, df.batch)]
D = D[key(D)]; c = c[key(c)]; e = e[key(e) & (e.eval_set == "main")]
stages = ["designs", "geometry gate", "full core filter", "strict Boltz gate", "+ PLACER signal", "shortlist (all)"]
f, axs = fig(183, 78, ncols=2, sharey=True); f.subplots_adjust(wspace=0.08)
for ax, m in zip(axs, ["hgnt", "hgwt"]):
    for a in ARMS:
        d = D[(D.motif == m) & (D.arm == a)]; n = len(d)
        if n == 0: continue
        s = 300 / n; cc = c[(c.motif == m) & (c.arm == a)]; ee = e[(e.motif == m) & (e.arm == a)]
        v = [300, d.pass_geom.sum() * s, d.pass_core.sum() * s, cc.boltz_gate_strict.sum() * s, (ee.flag_strict_boltz & ee.flag_placer_signal).sum() * s, ee.flag_shortlist.sum() * s]
        x = np.arange(len(v)); v = np.array(v, float); ok = v > 0
        hollow = a == "unsteered"
        ax.plot(x[ok], v[ok], color=COL[a], lw=1.4, zorder=2 if hollow else 3, solid_capstyle="round")
        ax.plot(x[ok], v[ok], "o", ms=4.5, mfc="white" if hollow else COL[a], mec=COL[a], mew=1.0, zorder=4)
    ax.set_yscale("log"); ax.set_xticks(np.arange(len(stages))); ax.set_xticklabels(stages, rotation=35, ha="right", fontsize=6.5); ax.set_title(MOTLAB[m], fontsize=7, loc="left")
    ax.axhline(1, color=SPINE, lw=0.4, ls=":"); ax.axvline(2.5, color=SPINE, lw=0.5, ls="--")
    ax.text(1.2, 650, "designs", fontsize=6.5, color=NOTE, ha="center"); ax.text(4.0, 650, "LigandMPNN sequences", fontsize=6.5, color=NOTE, ha="center")
axs[0].set_ylabel("count per 300 generated designs"); axs[0].set_ylim(0.3, 700)
hs = [plt.Line2D([], [], color=COL[a], marker="o", mfc="white" if a == "unsteered" else COL[a], mec=COL[a], lw=1.4, ms=4.5) for a in ARMS]
nn = {a: int(((D.arm == a)).sum()) for a in ARMS}
legend_mono(axs[1], hs, [f"{LAB[a]:<24} n={nn[a]:>5} designs" for a in ARMS], title="arm (both motifs pooled in n)")
footnote(f, f"Source: 20261122_kemp_c1, {len(fu)} runs with stage B finished. Shortlist = strict Boltz + PLACER (apo>=0.3 or holo>=0.1) + ipSAE2(TSA)>=0.7 + TSA>substrate (z>2) + dG_sep<=-15.\nUnsteered runs are gated on geometry only; steered on the full core filter. Sheet runs excluded.", y=-0.2)
save(f, 1, "funnel")
