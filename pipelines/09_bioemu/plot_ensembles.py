"""Figures for analyze_ensembles.py output. Run with the ipsae2 env (has matplotlib):
   /home/ec2-user/.conda/envs/ipsae2/bin/python plot_ensembles.py"""
import csv, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H = Path(__file__).resolve().parent / "analysis"
RUNS = {"5_mpnn1": "hgnt__steer_base_tmol__no_trp_0_beam3_5__mpnn1",
        "5_mpnn4": "hgnt__steer_base_tmol__no_trp_0_beam3_5__mpnn4",
        "6_mpnn1": "hgnt__steer_base_tmol__no_trp_0_beam3_6__mpnn1"}
COL = {"5_mpnn1": "#1b6ca8", "5_mpnn4": "#d95f02", "6_mpnn1": "#2a9d55"}
DESIGN = {"b_dist": 2.67, "d_dist": 3.22, "geometry": 0.87}      # measured on the RFD4 design 5 (6: 2.70 / 3.26 / 0.83)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def table(path):
    r = list(csv.DictReader(open(path)))
    return {k: np.array([float(x[k]) if x[k] not in ("", "nan") else np.nan for x in r]) if k not in ("aa", "role") else [x[k] for x in r]
            for k in r[0]}


F = {k: table(H / v / "frames.csv") for k, v in RUNS.items()}
R = {k: table(H / v / "residues.csv") for k, v in RUNS.items()}
S = {k: json.load(open(H / v / "summary.json")) for k, v in RUNS.items()}
out = H / "figures"; out.mkdir(exist_ok=True)

# 1 RMSF + ligand contacts
fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
for k in RUNS:
    ax[0].plot(R[k]["res"], R[k]["rmsf"], color=COL[k], lw=1.4, label=k)
    ax[1].plot(R[k]["res"], R[k]["contact_frac_bound"], color=COL[k], lw=1.2, label=k)
for k, role in (("5_mpnn1", "base"), ("5_mpnn1", "donor")):
    i = R[k]["role"].index(role)
    for a in ax:
        a.axvline(R[k]["res"][i], color="k", ls=":", lw=0.8)
    ax[0].text(R[k]["res"][i] + 0.5, ax[0].get_ylim()[1] * 0.92, f"{role}\n{R[k]['aa'][i]}{int(R[k]['res'][i])}", fontsize=8)
ax[0].set_ylabel("CA RMSF (A), core-aligned"); ax[0].legend(frameon=False, ncol=3)
ax[0].set_title("Backbone flexibility: stable helical core, free C-terminus")
ax[1].set_ylabel("fraction of ligand-bound samples\nwith residue within 4 A of ligand"); ax[1].set_xlabel("residue")
fig.tight_layout(); fig.savefig(out / "fig1_rmsf_contacts.png", dpi=150); plt.close(fig)

# 2 ligand occupancy
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for k in RUNS:
    ax[0].hist(F[k]["lig_pose_rmsd"], bins=np.linspace(0, 40, 81), histtype="step", color=COL[k], lw=1.5, label=f"{k}  bound {S[k]['occupancy']['frac_bound']:.0%}")
ax[0].axvline(3, color="k", ls=":", lw=0.8); ax[0].set_xlabel("ligand heavy-atom RMSD to designed pose (A)"); ax[0].set_ylabel("samples"); ax[0].legend(frameon=False)
ax[0].set_title("Ligand is bimodal: bound or roaming")
for k in RUNS:
    ax[1].scatter(F[k]["lig_com_shift"], F[k]["lig_min_dist"], s=3, alpha=0.35, color=COL[k], label=k)
ax[1].set_xlabel("ligand centre-of-mass shift from design (A)"); ax[1].set_ylabel("closest protein heavy atom (A)")
ax[1].set_title("roaming ligand still touches the surface")
kk = list(RUNS)
ax[2].bar(range(3), [S[k]["occupancy"]["frac_bound"] for k in kk], color=[COL[k] for k in kk])
ax[2].set_xticks(range(3)); ax[2].set_xticklabels(kk); ax[2].set_ylabel("fraction of 1000 samples (RMSD<=3 A)")
ax[2].set_title("pocket occupancy")
fig.tight_layout(); fig.savefig(out / "fig2_ligand_occupancy.png", dpi=150); plt.close(fig)

# 3 catalytic geometry of bound samples
fig, ax = plt.subplots(1, 4, figsize=(14, 3.5))
for j, (key, lab, hi) in enumerate((("b_dist", "base O -> N3 (A)", 14), ("d_dist", "Gln NE2 -> N1 (A)", 10),
                                    ("geometry", "geometry score", 1), ("burial_CG", "CA within 12 A of base CG", 50))):
    for k in RUNS:
        m = F[k]["bound"] == 1
        ax[j].hist(F[k][key][m], bins=np.linspace(0 if key != "burial_CG" else 15, hi, 36), histtype="step", color=COL[k], lw=1.5, label=k, density=True)
    if key in DESIGN:
        ax[j].axvline(DESIGN[key], color="k", ls="--", lw=1)
    if key == "burial_CG":
        ax[j].axvline(41, color="k", ls="--", lw=1)
    if key in ("b_dist", "d_dist"):
        ax[j].axvline(3.5, color="r", ls=":", lw=0.8)
    ax[j].set_xlabel(lab); ax[j].set_yticks([])
ax[0].legend(frameon=False); ax[0].set_title("bound samples only; dashed = RFD4 design")
fig.tight_layout(); fig.savefig(out / "fig3_catalytic_geometry_bound.png", dpi=150); plt.close(fig)

# 4 PCA
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for j, k in enumerate(RUNS):
    sc = ax[j].scatter(F[k]["pc1"], F[k]["pc2"], c=F[k]["burial_CG"], s=5, cmap="viridis")
    b = F[k]["bound"] == 1
    ax[j].scatter(F[k]["pc1"][b], F[k]["pc2"][b], s=14, facecolors="none", edgecolors="r", linewidths=0.4)
    v = S[k]["pca_var_explained"]
    ax[j].set_title(f"{k}  PC1 {v[0]:.0%}, PC2 {v[1]:.0%}")
    ax[j].set_xlabel("PC1 (CA)"); ax[j].set_ylabel("PC2 (CA)")
    plt.colorbar(sc, ax=ax[j], label="CA within 12 A of base")
fig.suptitle("Backbone modes (red rings = ligand bound); PCs include the floppy termini; PC1 correlates with base burial (r -0.4 to -0.6)", fontsize=9)
fig.tight_layout(); fig.savefig(out / "fig4_pca_backbone.png", dpi=150); plt.close(fig)

# 5 where the catalytic atoms move
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
keys = ["shift_N3", "shift_N1", "shift_baseO", "shift_donorNE2"]
w = 0.26
for j, k in enumerate(RUNS):
    ax[0].bar(np.arange(4) + (j - 1) * w, [S[k]["bound_only"]["shift_from_design_median"][x] for x in keys], w, color=COL[k], label=k)
ax[0].set_xticks(range(4)); ax[0].set_xticklabels(["ligand N3\n(H abstracted)", "ligand N1\n(oxyanion)", "Asp Oδ mid", "Gln NE2"])
ax[0].set_ylabel("median displacement from design (A)\nbound samples"); ax[0].legend(frameon=False)
ax[0].set_title("The N3 end and the Asp move; the N1/Gln end holds")
atoms = list(S["5_mpnn1"]["bound_only"]["ligand_rmsf_by_atom"])
for j, k in enumerate(RUNS):
    ax[1].bar(np.arange(len(atoms)) + (j - 1) * w, [S[k]["bound_only"]["ligand_rmsf_by_atom"][a] for a in atoms], w, color=COL[k])
ax[1].set_xticks(range(len(atoms))); ax[1].set_xticklabels(atoms, rotation=60); ax[1].set_ylabel("ligand atom RMSF, bound samples (A)")
ax[1].set_title("ligand mobility inside the pocket")
fig.tight_layout(); fig.savefig(out / "fig5_catalytic_atom_shifts.png", dpi=150); plt.close(fig)

# 6 b_dist vs d_dist
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6), sharex=True, sharey=True)
for j, k in enumerate(RUNS):
    b = F[k]["bound"] == 1
    ax[j].scatter(F[k]["b_dist"][~b], F[k]["d_dist"][~b], s=4, color="#bbbbbb", label="roaming")
    ax[j].scatter(F[k]["b_dist"][b], F[k]["d_dist"][b], s=7, color=COL[k], label="bound")
    ax[j].axvline(3.5, color="r", ls=":", lw=0.8); ax[j].axhline(3.5, color="r", ls=":", lw=0.8)
    ax[j].scatter([DESIGN["b_dist"]], [DESIGN["d_dist"]], marker="*", s=120, color="k", label="design", zorder=5)
    ax[j].set_xlim(2, 20); ax[j].set_ylim(2, 20); ax[j].set_title(k); ax[j].set_xlabel("base O -> N3 (A)")
ax[0].set_ylabel("Gln NE2 -> N1 (A)"); ax[0].legend(frameon=False, loc="upper left")
fig.tight_layout(); fig.savefig(out / "fig6_base_vs_donor_distance.png", dpi=150); plt.close(fig)
print("figures:", sorted(p.name for p in out.glob("*.png")))
