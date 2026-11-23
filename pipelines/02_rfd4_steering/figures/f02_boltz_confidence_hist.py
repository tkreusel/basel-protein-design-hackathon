"""f02: histograms of Boltz-2 / ipSAE2 confidence of the 6NT (TSA) refolds per arm and motif.  Data: ../../tables/evaluated_sequences_full_suite.csv (full-suite sequences, sheet runs excluded)."""
from kstyle import *
e = ev(); cols = [("ipsae.tsa", "ipSAE2 protein->TSA", None), ("boltz_ipde.tsa", "Boltz ipDE, TSA, mean of 10 (A)", None), ("6NT.ligand_iptm.mean", "ligand ipTM, TSA, mean of 10", None), ("6NT.pocket_rmsd.mean", "pocket RMSD to design (A)", None)]
f, axs = fig(183, 98, nrows=2, ncols=4); f.subplots_adjust(hspace=0.55, wspace=0.18, right=0.80)
for j, (c, lab, xl) in enumerate(cols):
    hist_panels(axs[:, j], e, c, lab, xlim=xl)
    for i, m in enumerate(("hgnt", "hgwt")):
        if j == 0: axs[i, j].set_ylabel(MOTLAB[m].split(" ")[0] + "\ndensity", fontsize=7)
arm_legend(axs[0, 3], e); 
footnote(f, "Source: 20261122_kemp_c1, evaluated sequences (strict + fallback survivors and extended set); triangles = medians; histograms normalised per arm. Higher ipSAE2/ipTM and lower ipDE/RMSD are better.", y=-0.07)
save(f, 2, "boltz_confidence_hist")
