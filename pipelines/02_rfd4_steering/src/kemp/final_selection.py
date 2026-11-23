"""Final selection: geometry checks on the strict + ipSAE2 set, then the top 30 by composite score.

usage: final_selection.py [--camp NAME]     needs tables/final_strict_all_gates_ipsae0.8.csv (export_results.py) and the run directories

For every sequence of final_strict_all_gates_ipsae0.8 the Boltz TSA samples and the RFD4 design are checked against the theozyme:
  oxyanion donor bound   catalytic Gln NE2 .. 6NT leaving-group nitrogen N1 <= 3.5 A in the best TSA sample and in >= 50 % of the samples
                         (hgwt: Trp44 NE1 .. N1 is reported too)
Used to filter, then rank by composite_z and take the first 30. Reported only: base burial (C-beta atoms within 12 A of the base carboxylate carbon),
catalytic-tip RMSD to the theozyme after superposing 6NT onto it, Asp O .. N3 distance, polar atoms within 3.5 A of N1, and sample-to-sample
consistency (ligand and tip RMSD to the best sample after CA superposition).
"""
import argparse
import glob
import os
import shutil

import numpy as np
import pandas as pd

from . import config, structure
from .export_results import all_sequences


def res_xyz(atoms, resid, names):
    """Coordinates of the named protein atoms of one residue."""
    return np.array([atoms.coord[(atoms.res_id == resid) & (atoms.atom_name == n) & ~atoms.hetero][0] for n in names])


def ligand_xyz(atoms, names):
    return np.array([atoms.coord[atoms.hetero & (atoms.atom_name == n)][0] for n in names])


def min_dist(A, B):
    return float(np.linalg.norm(A[:, None] - B[None], axis=-1).min())


def tips_xyz(atoms, resids, kinds):
    tips = structure.residue_tips()
    return np.concatenate([res_xyz(atoms, r, tips[k]) for r, k in zip(resids, kinds)])


def cb_count(atoms, base_resid, radius=12.0):
    """C-beta atoms (CA for Gly) within `radius` of the base carboxylate carbon, excluding the base residue itself."""
    prot = atoms[~atoms.hetero]
    cb = []
    for rid in np.unique(prot.res_id):
        name = "CA" if prot.res_name[prot.res_id == rid][0] == "GLY" else "CB"
        cb.append(prot.coord[(prot.res_id == rid) & (prot.atom_name == name)][0])
    cg = res_xyz(atoms, base_resid, ["CG"])[0]
    return int((np.linalg.norm(np.array(cb) - cg, axis=1) <= radius).sum()) - 1


def geometry_metrics(design, samples, fixed, motif, theozyme):
    """All geometry measurements for one sequence. `samples[0]` is the best TSA sample; `fixed` = motif residue ids (Asp, Gln[, Trp])."""
    motifs = config.load("motifs")
    kinds = motifs["motifs"][motif]["residues"]
    lig_names = motifs["ligand_atoms"]
    th_res = {k: v["resnum"] for k, v in motifs["residues"].items()}
    th_tips = tips_xyz(theozyme, [th_res[k] for k in kinds], kinds)
    th_lig = ligand_xyz(theozyme, lig_names)

    def tips_rmsd_vs_theozyme(atoms):
        R, t = structure.kabsch(ligand_xyz(atoms, lig_names), th_lig)
        return structure.rmsd_xyz(tips_xyz(atoms, fixed, kinds) @ R.T + t, th_tips)

    best = samples[0]
    prot, lig = best[~best.hetero], best[best.hetero]
    asp, gln = fixed[0], fixed[1]
    n1 = lig.coord[lig.atom_name == "N1"]
    out = dict(n_tsa_samples=len(samples), design_tips_rmsd_theozyme=tips_rmsd_vs_theozyme(design), boltz_best_tips_rmsd_theozyme=tips_rmsd_vs_theozyme(best),
               asp_O_to_N3_A=min_dist(res_xyz(best, asp, motifs["base_oxygens"]), lig.coord[lig.atom_name == "N3"]),
               gln_NE2_to_N1_A=min_dist(res_xyz(best, gln, [motifs["donor_atom"]]), n1))
    if motif == "hgwt":
        out["trp_NE1_to_N1_A"] = min_dist(res_xyz(best, fixed[2], ["NE1"]), n1)
    from scipy.spatial import cKDTree
    out["n_polar_atoms_3p5_of_N1"] = int(len(cKDTree(prot.coord[np.isin(prot.element, ("N", "O"))]).query_ball_point(n1[0], 3.5)))
    out["base_CB_within_12A_boltz"], out["base_CB_within_12A_design"] = cb_count(best, asp), cb_count(design, asp)

    ca0 = best[~best.hetero & (best.atom_name == "CA")].coord
    lig_rmsd, tip_rmsd, donor = [], [], []
    for s in samples:
        ca = s[~s.hetero & (s.atom_name == "CA")].coord
        R, t = structure.kabsch(ca, ca0)
        lig_rmsd.append(structure.rmsd_xyz(ligand_xyz(s, lig_names) @ R.T + t, ligand_xyz(best, lig_names)))
        tip_rmsd.append(structure.rmsd_xyz(tips_xyz(s, fixed, kinds) @ R.T + t, tips_xyz(best, fixed, kinds)))
        donor.append(min_dist(res_xyz(s, gln, [motifs["donor_atom"]]), s.coord[s.hetero & (s.atom_name == "N1")]))
    donor = np.array(donor)
    out.update(tsa_samples_ligand_rmsd_median=float(np.median(lig_rmsd)), tsa_samples_ligand_rmsd_max=float(np.max(lig_rmsd)),
               tsa_samples_tips_rmsd_median=float(np.median(tip_rmsd)), tsa_samples_tips_rmsd_max=float(np.max(tip_rmsd)),
               frac_samples_gln_donor_bound_3p5=float(np.mean(donor <= 3.5)), gln_NE2_to_N1_median_A=float(np.median(donor)))
    return out


def sort_key(path):
    return int(path.rsplit("_", 1)[1][:-4])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--camp")
    a = ap.parse_args()
    cfg = config.load("gates")["final_selection"]
    root, results = config.runs_root(a.camp), config.results_dir()
    group = results / "final_strict_all_gates_ipsae0.8"
    strict = pd.read_csv(results / "tables" / "final_strict_all_gates_ipsae0.8.csv")
    theozyme = structure.load_reference()
    rows = []
    for x in strict.itertuples():
        run = f"{root}/{x.motif}/{x.arm}/{x.batch}"
        forest, ev_dir = (f"{run}/eval_forest_x", f"{run}/eval_x") if x.eval_set == "extended" else (f"{run}/eval_forest", f"{run}/eval")
        sq = pd.read_csv(run + "/mpnn/sequences.csv").set_index("seq_id").loc[x.seq_id]
        fixed = [int(t[1:]) for t in str(sq.fixed).split()]
        paths = sorted(glob.glob(f"{forest}/{x.seq_id}__6NT/{x.seq_id}__6NT_model_*.cif"), key=sort_key)
        g = geometry_metrics(structure.load_cif(f"{run}/gen/designs/{x.design}.cif"), [structure.load_cif(p) for p in paths], fixed, x.motif, theozyme)
        rows.append(dict(seq_id=x.seq_id, forest=forest, ev_dir=ev_dir, **g))
    m = strict.merge(pd.DataFrame(rows), on="seq_id")
    m["oxyanion_donor_bound"] = ((m.gln_NE2_to_N1_A <= cfg["donor_bound_distance"]) & (m.frac_samples_gln_donor_bound_3p5 >= cfg["donor_bound_fraction"]))
    m.to_csv(group / "strict_all_gates_with_geometry.csv", index=False)
    print("strict + ipSAE set:", len(m), "| oxyanion donor bound:", int(m.oxyanion_donor_bound.sum()), "| best sample only:", int((m.gln_NE2_to_N1_A <= 3.5).sum()))

    cand = pd.read_csv(results / "tables" / "candidates.csv")[["seq_id", "b_base_sasa_rel", "b_volume_ratio", "b_bottleneck_radius"]].drop_duplicates("seq_id")
    m = m.merge(cand, on="seq_id", how="left")
    m["scaffold"] = m.design.str.replace(r"_beam\d+_\d+$", "", regex=True).str.replace(r"^(hgnt|hgwt)__(steer_tmol_ch18|steer_tmol|unsteered)__", "", regex=True)
    ok = m[m.oxyanion_donor_bound].sort_values("composite_z", ascending=False).copy()
    ok["rank30"] = np.arange(1, len(ok) + 1)
    top = ok.head(cfg["top_n"]).copy()
    cols = [c for c in ["rank30", "motif", "arm", "variant", "scaffold", "composite_z", "ipsae.tsa", "ipsae.tsa|sd", "boltz_ipde.tsa", "6NT.ligand_iptm.mean", "ros.tsa.dG_sep",
                        "delta.tsa_minus_sub.total", "placer.apo.frac_pass_core", "placer.holo.frac_pass_core", "gln_NE2_to_N1_A", "gln_NE2_to_N1_median_A",
                        "frac_samples_gln_donor_bound_3p5", "n_polar_atoms_3p5_of_N1", "asp_O_to_N3_A", "design_tips_rmsd_theozyme", "boltz_best_tips_rmsd_theozyme",
                        "base_CB_within_12A_boltz", "base_CB_within_12A_design", "b_base_sasa_rel", "tsa_samples_ligand_rmsd_median", "tsa_samples_ligand_rmsd_max",
                        "tsa_samples_tips_rmsd_median", "tsa_samples_tips_rmsd_max", "ligand_rmsd_ca_aligned", "motif_tip_rmsd_ca_aligned", "z_tsa_vs_sub", "z_tsa_vs_prod", "seq_id"]
            if c in top]
    out = group / "top30"
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    top[cols].round(3).to_csv(out / "top30_table.csv", index=False)
    ok.round(3).to_csv(out / "all_donor_bound_ranked.csv", index=False)
    seqs = all_sequences(root)
    with open(out / "top30.fasta", "w") as fh:
        for x in top.itertuples():
            fh.write(f">{x.rank30:02d}_{x.seq_id} {x.motif} {x.arm} {x.variant}\n{seqs.get(x.seq_id, '')}\n")
    print("top30:", len(top), "| scaffolds:", top.scaffold.nunique(), "| arms:", top.arm.value_counts().to_dict(), "| donor-bound pool:", len(ok))


if __name__ == "__main__":
    main()
