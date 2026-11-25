"""f07: design level (all generated non-sheet designs): channel bottleneck radius vs TMol dG per arm.  Shaded band = steering target 1.7-1.9 A; dark points pass the full core filter."""
from kstyle import *
d = designs().dropna(subset=["bottleneck_radius", "tmol_ddg"]); f, axs = fig(183, 70, ncols=3, sharex=True, sharey=True); f.subplots_adjust(wspace=0.06)
rng = np.random.default_rng(0)
for ax, a in zip(axs, ARMS):
    x = d[d.arm == a]; jx = rng.normal(0, 0.012, len(x))
    ax.axvspan(1.7, 1.9, color="#e2661a", alpha=0.10, lw=0)
    ax.scatter(x.bottleneck_radius + jx, x.tmol_ddg.clip(-45, 10), s=3, color=COL[a], alpha=0.18, lw=0, rasterized=True)
    p = x[x.pass_core]; ax.scatter(p.bottleneck_radius + rng.normal(0, 0.012, len(p)), p.tmol_ddg.clip(-45, 10), s=5, color=COL[a], edgecolor=INK, linewidth=0.2, alpha=0.9, rasterized=True)
    ax.set_title(f"{LAB[a]}\nn={len(x)} designs, {len(p)} pass core", fontsize=6.5, loc="left"); ax.set_xlabel("channel bottleneck radius (A)")
axs[0].set_ylabel("TMol dG, protein-6NT (REU)\n(lower is better)")
footnote(f, "Source: 20261122_kemp_c1, all generated designs of finished runs (sheet runs excluded). Bottleneck = largest probe radius that connects the 6NT pocket to bulk (discrete ladder, jittered); band = 1.7-1.9 A target of the channel arm.", y=-0.08)
save(f, 7, "channel_vs_tmol_designs")
