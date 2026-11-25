"""f08: which Boltz-1x gate criteria limit the survivors: pass rate per criterion and arm (all LigandMPNN sequences that were folded)."""
from kstyle import *
c = cands(); cr = [("ok_plddt", "pLDDT >= 80"), ("ok_ca", "CA RMSD <= 2 A"), ("ok_liptm", "ligand ipTM >= 0.6"), ("ok_ligand", "ligand RMSD <= 2 A"), ("ok_motif", "motif RMSD <= 1 A"), ("geom_ok_refold", "refold geometry"), ("boltz_gate_strict", "all (strict gate)")]
f, ax = fig(150, 75); f.subplots_adjust(right=0.62); w = 0.26
for i, a in enumerate(ARMS):
    d = c[c.arm == a]; v = [100 * d[k].mean() for k, _ in cr]
    ax.bar(np.arange(len(cr)) + (i - 1) * w, v, w * 0.92, color="white" if a == "unsteered" else COL[a], edgecolor=COL[a], linewidth=1.0 if a == "unsteered" else 0)
ax.set_xticks(np.arange(len(cr))); ax.set_xticklabels([l for _, l in cr], rotation=30, ha="right"); ax.set_ylabel("sequences passing (%)"); ax.set_ylim(0, 105)
hs = [plt.Rectangle((0, 0), 1, 1, fc="white" if a == "unsteered" else COL[a], ec=COL[a]) for a in ARMS]
legend_mono(ax, hs, [f"{LAB[a]:<22} n={int((c.arm == a).sum())}" for a in ARMS], title="sequences into Boltz")
footnote(f, "Source: 20261122_kemp_c1, Boltz-2 1x TSA refolds of all LigandMPNN sequences of finished runs (both motifs pooled; sheet runs excluded).", y=-0.3)
save(f, 8, "boltz_gate_criteria")
