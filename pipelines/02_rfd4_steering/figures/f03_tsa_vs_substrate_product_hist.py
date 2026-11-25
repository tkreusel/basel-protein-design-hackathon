"""f03: ligand discrimination - TSA minus substrate (row 1) and TSA minus product (row 2) for ipTM, ipSAE2, ipDE, tmol energy and PyRosetta binding energy.  Dashed zero line = no preference."""
from kstyle import *
import textwrap
e = ev().copy()
e["d_iptm_s"] = e["6NT.iptm.mean"] - e["H5J.iptm.mean"]; e["d_iptm_p"] = e["6NT.iptm.mean"] - e["6VP.iptm.mean"]
e["d_ipsae_s"] = e["ipsae.tsa"] - e["ipsae.sub"]; e["d_ipsae_p"] = e["ipsae.tsa"] - e["ipsae.prod"]
e["d_ipde_s"] = e["boltz_ipde.tsa"] - e["boltz_ipde.sub"]; e["d_ipde_p"] = e["boltz_ipde.tsa"] - e["boltz_ipde.prod"]
e["d_tm_s"] = e["delta.tsa_minus_sub.total"]; e["d_tm_p"] = e["delta.tsa_minus_sub.total"] - e["delta.prod_minus_sub.total"]
e["d_ros_s"] = e["ros.delta_native.tsa_minus_sub.dG_sep"]; e["d_ros_p"] = e["ros.delta_native.tsa_minus_sub.dG_sep"] - e["ros.delta_native.prod_minus_sub.dG_sep"]
cols = [("iptm", "ipTM (TSA - X)"), ("ipsae", "ipSAE2 (TSA - X)"), ("ipde", "ipDE, A (TSA - X)"), ("tm", "tmol dE, REU (TSA - X)"), ("ros", "PyRosetta dG_sep, REU (TSA - X)")]
f, axs = fig(183, 104, nrows=2, ncols=5); f.subplots_adjust(hspace=0.75, wspace=0.18, right=0.80)
for j, (k, lab) in enumerate(cols):
    for i, s in enumerate(("s", "p")):
        c = f"d_{k}_{s}"; v = e[c].dropna(); lo, hi = np.percentile(v, [1, 99]); bins = np.linspace(min(lo, 0), max(hi, 0), 16)
        hist_panels(axs[i, j], e.assign(**{c: e[c].clip(bins[0], bins[-1])}), c, lab.replace("X", "substrate" if s == "s" else "product"), bins=bins, motifs=("hgnt",) if False else ("hgnt", "hgwt")[:0] or ("hgnt",), zero=True) if False else None
        # pooled over motifs (n is small); motif shown by colour of nothing - arms only
        for a in ARMS:
            vv = e[e.arm == a][c].dropna().clip(bins[0], bins[-1])
            axs[i, j].hist(vv, bins=bins, histtype="stepfilled", color=COL[a], alpha=0.16, lw=0, density=True); axs[i, j].hist(vv, bins=bins, histtype="step", color=COL[a], lw=1.3 if a != "unsteered" else 1.0, ls="--" if a == "unsteered" else "-", density=True)
            axs[i, j].plot([vv.median()], [0], marker="v", ms=4, color=COL[a], mec="white", mew=0.5, clip_on=False, zorder=5)
        axs[i, j].axvline(0, color=INK, lw=0.6, ls=":"); axs[i, j].set_yticks([]); axs[i, j].set_xlabel(textwrap.fill(lab.replace("X", "substrate" if s == "s" else "product"), 16), fontsize=6.5)
axs[0, 0].set_ylabel("TSA vs substrate\ndensity"); axs[1, 0].set_ylabel("TSA vs product\ndensity")
hs = [plt.Line2D([], [], color=COL[a], lw=1.4, ls="--" if a == "unsteered" else "-") for a in ARMS]
legend_mono(axs[0, 4], hs, [f"{LAB[a]:<22} n={int((e.arm == a).sum())}" for a in ARMS], title="evaluated seqs\n(both motifs)")
footnote(f, "Source: 20261122_kemp_c1 evaluated sequences. Boltz values = mean of 10 samples per ligand. Preference for the TSA = ipTM/ipSAE2 > 0 (right), ipDE/energies < 0 (left). tmol dE(TSA-product) = dE(TSA-sub) - dE(product-sub). Triangles = medians; 1-99 % range shown.", y=-0.07)
save(f, 3, "tsa_vs_substrate_product_hist")
