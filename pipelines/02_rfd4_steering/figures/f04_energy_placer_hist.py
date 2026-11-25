"""f04: PyRosetta / tmol energies and PLACER catalytic-geometry fractions per arm and motif."""
from kstyle import *
e = ev(); cols = [("ros.tsa.dG_sep", "PyRosetta dG_sep, TSA (REU)"), ("strain.E_noncst", "tmol strain to impose theozyme (REU)"), ("placer.apo.frac_pass_core", "PLACER apo: fraction pre-formed"), ("placer.holo.frac_pass_core", "PLACER holo: fraction in catalytic geometry")]
f, axs = fig(183, 98, nrows=2, ncols=4); f.subplots_adjust(hspace=0.55, wspace=0.18, right=0.80)
for j, (c, lab) in enumerate(cols):
    hist_panels(axs[:, j], e, c, lab, xlim=(0, 1) if "placer" in c else None)
    for i, m in enumerate(("hgnt", "hgwt")):
        if j == 0: axs[i, j].set_ylabel(m + "\ndensity", fontsize=7)
arm_legend(axs[0, 3], e)
footnote(f, "Source: 20261122_kemp_c1 evaluated sequences; PLACER = 20-model ensembles of the best Boltz sample (TSA); triangles = medians. Lower dG_sep / strain are better; PLACER fractions have a ~0.1 binomial uncertainty.", y=-0.07)
save(f, 4, "energy_placer_hist")
