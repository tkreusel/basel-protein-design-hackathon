"""f06: raw scatter PLACER apo (pre-formed) vs holo (in catalytic geometry); rings = strict Boltz gate; small jitter (+-0.01) because fractions are multiples of 0.05."""
from kstyle import *
e = ev().dropna(subset=["placer.apo.frac_pass_core", "placer.holo.frac_pass_core"]); f, ax = fig(150, 80); f.subplots_adjust(right=0.64)
scat(ax, e, "placer.holo.frac_pass_core", "placer.apo.frac_pass_core", jitter=0.008, ring=e.flag_strict_boltz)
ax.axvline(0.1, color=INK, lw=0.6, ls="--"); ax.axhline(0.3, color=INK, lw=0.6, ls="--"); ax.set_xlabel("PLACER holo: fraction of 20 models in catalytic geometry"); ax.set_ylabel("PLACER apo: fraction pre-formed")
scat_legend(ax, e, f"strict Boltz gate (n={int(e.flag_strict_boltz.sum())})")
footnote(f, "Source: 20261122_kemp_c1 evaluated sequences; PLACER ensembles of the best Boltz TSA sample; dashed = signal thresholds (holo 0.1, apo 0.3); binomial uncertainty ~0.1.", y=-0.06)
save(f, 6, "placer_holo_vs_apo_scatter")
