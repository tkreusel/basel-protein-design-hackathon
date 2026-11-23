"""f09: scaffold population of designs passing the full core filter: length vs helix fraction per arm."""
from kstyle import *
d = designs(); p = d[d.pass_core]; f, ax = fig(150, 80); f.subplots_adjust(right=0.64)
rng = np.random.default_rng(2); q = p.assign(len_j=p.length + rng.normal(0, 0.6, len(p)))
scat(ax, q, "len_j", "helix", size=9)
ax.set_xlabel("design length (residues)"); ax.set_ylabel("helix fraction of the design")
scat_legend(ax, q)
footnote(f, "Source: 20261122_kemp_c1, designs passing the full core filter (steered gate) or the core criteria (unsteered, for comparison); helix = biotite P-SEA on CA. Sheet runs excluded.", y=-0.06)
save(f, 9, "scaffold_population")
