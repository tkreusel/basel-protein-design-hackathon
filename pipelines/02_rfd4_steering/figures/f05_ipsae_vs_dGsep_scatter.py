"""f05: raw scatter ipSAE2 (TSA) vs PyRosetta dG_sep (TSA).  Rings = shortlist (strict Boltz + PLACER + ipSAE2>=0.7 + TSA>sub + dG_sep<=-15).  Dashed lines = shortlist thresholds."""
from kstyle import *
e = ev().dropna(subset=["ipsae.tsa", "ros.tsa.dG_sep"]); f, ax = fig(150, 80); f.subplots_adjust(right=0.64)
scat(ax, e, "ipsae.tsa", "ros.tsa.dG_sep", ring=e.flag_shortlist)
ax.axvline(0.7, color=INK, lw=0.6, ls="--"); ax.axhline(-15, color=INK, lw=0.6, ls="--"); ax.set_xlabel("ipSAE2, protein -> TSA (6NT)"); ax.set_ylabel("PyRosetta dG_sep, TSA (REU)\n(lower is better)")
scat_legend(ax, e, f"shortlist (n={int(e.flag_shortlist.sum())})")
footnote(f, "Source: 20261122_kemp_c1 evaluated sequences (sheet runs excluded); each point = one LigandMPNN sequence (best Boltz sample for PyRosetta).", y=-0.06)
save(f, 5, "ipsae_vs_dGsep_scatter")
