"""Build the hand-off deliverables: datasets (CSV), correlation table (CSV), data dictionary (CSV) and a PDF report.

  PYTHONPATH=/mnt/work/envs/pdfpkgs python -m maxi_eval.make_deliverables

Inputs (all produced earlier):
  work/report_tmolmean/summary_entries.csv   per-protein table, tmol = mean of the two independent relax runs
  work/features_extra_mean.csv               extra features (PLACER ensembles, tmol per-residue, Boltz PAE, composition), tmol-mean
  work/report/everything_tmolmean{,_HG3}/    selection-aware screens (all proteins / HG3 only) on the mean data
  work/report/everything/                    the same screen on tmol run 1 only,  work/report/everything_tmolrun2/ on run 2 only
  work/report/retest.csv                     tmol test-retest
"""
from __future__ import annotations
import argparse, glob, re, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------ glossary (also rendered into the PDF)
TYPE = {"tsa": "TSA complex (6NT)", "sub": "substrate complex (H5J)", "prod": "product complex (6VP)"}
TOOL = {
    "tmol swap": "tmol, all three ligands placed in the SAME (TSA-relaxed) Boltz pocket (same-pocket comparison)",
    "tmol native": "tmol on the independent Boltz refold made with that ligand",
    "ipSAE2/Boltz": "Boltz-2 confidence / ipSAE2 (protein→ligand)",
    "PLACER holo": "PLACER ensemble with the ligand present (TSA arm)",
    "Rosetta swap": "PyRosetta, same-pocket (swap) ligands", "Rosetta native": "PyRosetta on the independent refolds",
}
TOK = {  # token -> explanation   (used for descriptions of the candidate ids and shown in the PDF glossary)
    # arms / prefixes
    "raw": "unrelaxed Boltz structure with the TSA (no tmol relax)", "rlx": "TSA complex after tmol relax", "rlxcst": "TSA complex after tmol relax WITH the theozyme constraints",
    "tsa_pk": "TSA in the TSA-relaxed pocket (= 'tsa' arm of the same-pocket comparison)", "sub": "substrate (H5J) posed on the TSA frame in the same pocket",
    "prod": "product (6VP) posed on the TSA frame in the same pocket", "sub_ts": "substrate forced into the TS-like geometry (TS-pose control)",
    "prod_ts": "product forced into the TS-like geometry (TS-pose control)", "native_sub": "substrate in its own independent Boltz refold", "native_prod": "product in its own independent Boltz refold",
    "strain": "energy to impose the theozyme geometry (constrained minus free)", "apo": "apo (ligand-free) Boltz refold / PLACER apo ensemble", "ros": "PyRosetta quantity",
    "ip": "Boltz/ipSAE2 confidence quantity", "ipsae": "ipSAE2 quantity", "site": "site-focused Boltz confidence (catalytic residues ↔ ligand)", "ah": "apo-vs-holo comparison (Boltz apo vs TSA holo pocket)",
    "delta": "difference between arms (tmol/Rosetta legacy naming; tsa_minus_sub = TSA − substrate)", "diff": "difference between arms (summary naming)", "delta_native": "arm difference on independent refolds",
    "placer": "PLACER ensemble (TSA arm)", "placer_sub": "PLACER ensemble with the substrate", "placer_prod": "PLACER ensemble with the product", "placer_apoarm": "PLACER ensemble on the Boltz apo model",
    "pl2": "PLACER ensemble distribution features (per model)", "pr": "tmol per-residue ligand-interaction energies", "bz": "Boltz PAE / pLDDT matrices", "cmp": "sequence / pocket composition",
    "tsa_holo": "TSA-bound ensemble", "sub_holo": "substrate-bound ensemble", "prod_holo": "product-bound ensemble", "tsa_apo": "PLACER apo ensemble (TSA-pocket)", "apoarm_apo": "PLACER on the Boltz-apo model",
    # tmol energies
    "int": "ligand ↔ protein pair energies (REU)", "total": "total", "fa_lk": "Lazaridis–Karplus implicit solvation (desolvation penalty) of ligand–pocket pairs", "lk_ball_iso": "isotropic part of the lk_ball solvation model (ligand–pocket)",
    "lk_ball": "anisotropic (water-orientation) lk_ball solvation", "fa_elec": "distance-dependent Coulomb electrostatics", "fa_ljatr": "Lennard-Jones attraction (packing)", "fa_ljrep": "Lennard-Jones repulsion (clashes)",
    "hbond": "hydrogen-bond energy", "rep_over_atr": "repulsion / attraction ratio (packing quality)", "E_total": "total tmol energy of the complex (REU, includes theozyme restraints)", "E_total_noncst": "total energy without the restraint term",
    "E_constraint": "theozyme restraint energy", "E_ligand_intra": "ligand internal strain energy", "E_protein_in_complex": "protein energy in the complex (without restraints)", "E_total_per_res": "total energy per residue",
    "pocket_E_per_res": "mean energy per pocket residue", "pocket_E_max": "worst pocket-residue energy", "res_E_max": "worst residue energy in the protein", "n_res_E_gt5": "number of residues with energy > 5 REU",
    "rep_inter_total": "inter-residue repulsion, whole protein", "rep_intra_total": "intra-residue repulsion, whole protein", "rep_inter_pocket": "inter-residue repulsion, pocket", "rep_intra_pocket": "intra-residue repulsion, pocket",
    # geometry
    "geo": "geometry", "contact_min": "closest ligand–protein heavy-atom contact (Å)", "lig_rel_sasa": "relative solvent-accessible surface of the ligand (burial: low = buried)", "n_shell": "residues in the first shell (≤ 4.5 Å)",
    "n_donors_3p5": "H-bond donors within 3.5 Å of the oxyanion O", "n_contacts": "protein atoms within contact distance of the ligand", "enclosure": "fraction of the ligand's sphere enclosed by protein",
    "base_sasa": "SASA of the catalytic-base carboxylate O (Å²)", "d_O_H": "base O···H3 distance (Å)", "d_O_X": "base O···N3/C3 (donor atom) distance (Å)", "ang_O_H_X": "O···H–X angle (°)", "tors_syn": "carboxylate O′–C–O···H torsion (0 = syn, 180 = anti)",
    "donor_best_d": "best oxyanion-donor distance (Å)", "cst": "theozyme constraints", "penalty": "constraint penalty (flat-bottom)", "n_pass": "number of theozyme constraints satisfied", "rms_z": "RMS constraint deviation in tolerance units",
    "dev": "deviation from the reference value", "z": "deviation in tolerance units",
    "buns": "buried unsatisfied polar atoms (H-bond-unsatisfied)", "n_unsat": "buried unsatisfied polars (pocket)", "n_unsat_ligand": "buried unsatisfied polars on the ligand", "n_buried_polar": "buried polar atoms",
    "void": "pocket void (cavity) proxy", "volume_A3": "void volume (Å³)", "n_points": "void grid points", "field": "electric-field projection from the environment (MV/cm; + = favourable to the developing charge)",
    "CH": "on the transferring C–H bond", "LG": "on the leaving-group O (C7A–O1)", "env": "environment incl. base", "env_nobase": "environment excluding the catalytic base", "vac": "vacuum Coulomb", "rdd4": "distance-dependent dielectric ε = 4r",
    # Boltz / ipSAE
    "p2l": "ipSAE protein→ligand (d0 = 2 Å)", "l2p": "ipSAE ligand→protein", "iPAE12": "interface PAE protein→ligand (Å, low = confident)", "iPAE21": "interface PAE ligand→protein (Å)",
    "iptm": "Boltz ipTM (complex)", "ligand_iptm": "Boltz ligand ipTM", "complex_iplddt": "Boltz complex interface pLDDT", "pae_cat_to_lig": "PAE catalytic residues → ligand (Å)", "pae_lig_to_cat": "PAE ligand → catalytic residues (Å)",
    "plddt_cat": "pLDDT of catalytic residues", "plddt_lig": "pLDDT of ligand atoms", "cat_to_pocket_pae": "mean PAE catalytic residues → pocket residues (Å): how well the model trusts the catalytic geometry relative to its pocket",
    "pocket_pae_int": "mean PAE inside the pocket (Å)", "pocket_to_rest_pae": "PAE pocket → rest of the protein (Å)", "rest_to_pocket_pae": "PAE rest → pocket (Å)", "cat_pae_int": "PAE among the catalytic residues (Å)",
    # PLACER
    "frac_pass_core": "fraction of ensemble models in the catalytic geometry window (HG3-derived; invalid for 1A53/KABLE)", "prmsd_mean": "PLACER predicted RMSD, mean (Å)", "prmsd_top10": "PLACER prmsd of best 10 %", "plddt_mean": "PLACER pLDDT mean",
    "lig_rmsd_to_input": "ligand RMSD of PLACER models to the input Boltz pose (Å)", "lig_rmsd_spread": "ligand-pose spread across the ensemble (Å)", "chi_sd_base": "circular SD of the base chi1 (°)", "chi_sd_wall": "circular SD of wall-residue chi (°)",
    "field_CH_mean": "mean C–H field over the ensemble", "field_CH_sd": "SD of the C–H field over the ensemble", "field_LG_mean": "mean leaving-group field", "field_LG_sd": "SD of the leaving-group field",
    "nac30": "near-attack fraction: base O···H3 ≤ 3.0 Å and angle ≥ 130°", "nac35": "near-attack fraction: ≤ 3.5 Å and angle ≥ 110°", "donor_med": "median oxyanion-donor distance (Å)", "angCOH_med": "median C–O···H angle (°)",
    "ang_C_O_H_mean": "mean C–O···H angle (°) in the ensemble", "cst_n_pass": "constraints met when the TSA is transplanted into the apo pocket", "pocket_bb_rmsd": "pocket backbone RMSD apo vs holo (Å)", "reorg_E": "apo → holo reorganisation energy (REU)",
    # Rosetta
    "dG_sep": "Rosetta dG_separated (REU)", "sc": "shape complementarity (ligand, manual-atom)", "dSASA": "buried SASA on binding (Å²)", "packstat_int": "interface packstat", "n_hb": "interface H-bonds", "hbond_E": "H-bond energy",
    "delta_unsat_hb": "change in unsatisfied H-bonds on binding",
    # composition
    "pocket_frac_pos": "fraction of Arg/Lys/His in the pocket", "pocket_frac_neg": "fraction of Asp/Glu in the pocket", "pocket_charge": "net charge of the pocket", "pocket_frac_arom": "aromatic fraction", "pocket_frac_hydro": "hydrophobic fraction",
    "pocket_frac_polar": "polar fraction", "n_pocket": "number of pocket residues", "pocket_n_gly": "Gly in pocket", "pocket_n_arom": "aromatic residues in pocket", "length": "sequence length", "n_mut_vs_HG3": "number of mutations relative to HG3",
    "int_polar": "ligand interaction energy with polar pocket residues", "int_sum": "summed per-residue ligand interaction energy", "int_base": "ligand interaction energy with the catalytic base residue", "int_neg": "interaction with negatively-charged residues",
    "hbond_sum": "summed per-residue H-bond energy", "pocket_E_sum": "summed energy of pocket residues", "n_res_gt1": "residues with |interaction| > 1 REU", "n_res_lt05": "residues with interaction < −0.5 REU", "n_hbond_res": "residues H-bonding the ligand",
}
DIFF = lambda a, b: f"difference {TYPE.get(a, a)} − {TYPE.get(b, b)}"


def describe(name: str) -> str:
    """Human-readable description of a candidate id."""
    if "|" in name:
        tool, met, ty = name.split("|")
        arm = DIFF(*ty.split("-")) if "-" in ty else TYPE.get(ty, ty)
        return f"{TOOL.get(tool, tool)}: {' · '.join(TOK.get(t, t) for t in met.split('.'))} — {arm}"
    toks = name.split(".")
    out = []
    for t in toks:
        if "-" in t and re.fullmatch(r"[a-z_]+-[a-z_]+", t):
            a, b = t.split("-"); out.append(DIFF(a, b))
        else:
            out.append(TOK.get(t, t))
    return " · ".join(out)


# ------------------------------------------------------------------ data
def load_ipae21(work: Path) -> pd.DataFrame:
    """iPAE ligand→protein (not in the master table) averaged over Boltz models, per protein and arm."""
    from .run_ipsae import prot_lig_row
    rows = []
    for f in glob.glob(str(work / "ipsae" / "*.tsv")):
        m = re.fullmatch(r"(.+)_boltz_(tsa|sub|prod)_boltz_m\d+\.tsv", Path(f).name)
        if not m:
            continue
        try:
            r = prot_lig_row(pd.read_csv(f, sep="\t"))
            rows.append(dict(protein=m.group(1), arm=m.group(2), iPAE21=float(r["iPAE21"])))
        except Exception:
            pass
    d = pd.DataFrame(rows)
    if d.empty:
        return pd.DataFrame()
    P = d.groupby(["protein", "arm"]).iPAE21.mean().unstack()
    out = pd.DataFrame(index=P.index)
    for a in ("tsa", "sub", "prod"):
        out[f"ip.{a}.ipsae.iPAE21"] = P[a]
    out["ipAE21|tsa-sub"] = P.tsa - P["sub"]; out["ipAE21|tsa-prod"] = P.tsa - P["prod"]; out["ipAE21|sub-prod"] = P["sub"] - P["prod"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="work/report/deliverables"); ap.add_argument("--labels", default="config/labels_manual_dataset.csv")
    ap.add_argument("--dataset", default="/mnt/work/manually_created_dataset.csv"); ap.add_argument("--top", type=int, default=300)
    a = ap.parse_args()
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True); work = Path("work")
    from .analyze_full import load_master
    ent = load_master(Path("work/report_tmolmean"), work, a.labels)
    ex = pd.read_csv(work / "features_extra_mean.csv", index_col=0)
    M = pd.read_csv(work / "report/everything_tmolmean/candidate_matrix.csv", index_col=0)
    M1 = pd.read_csv(work / "report/everything/candidate_matrix.csv", index_col=0)
    M2 = pd.read_csv(work / "report/everything_tmolrun2/candidate_matrix.csv", index_col=0)
    C = pd.read_csv(work / "report/everything_tmolmean/candidates_all.csv")
    C1 = pd.read_csv(work / "report/everything/candidates_all.csv").set_index("candidate")
    C2 = pd.read_csv(work / "report/everything_tmolrun2/candidates_all.csv").set_index("candidate")
    CH = pd.read_csv(work / "report/everything_tmolmean_HG3/candidates_all.csv").set_index("candidate")
    N = pd.read_csv(work / "report/everything_tmolmean/selection_aware_by_set.csv"); NH = pd.read_csv(work / "report/everything_tmolmean_HG3/selection_aware_by_set.csv")
    ds = pd.read_csv(a.dataset, sep="\t"); ds.columns = [c.strip() for c in ds.columns]
    ds = ds.rename(columns={"PDB ID": "pdb_ids", "Name": "name", "kcat/km": "kcat_KM_dataset", "source": "kcat_KM_source"})[["name", "pdb_ids", "kcat_KM_source"]]
    ds["protein"] = ds.name
    ipae = load_ipae21(work)

    # ---- identity block
    idb = ent[["family", "kcat_KM", "log10_kcat_KM"]].copy(); idb.insert(0, "protein", idb.index)
    idb = idb.merge(ds.set_index("protein")[["pdb_ids", "kcat_KM_source"]], left_on="protein", right_index=True, how="left"); idb.index = idb.protein
    idb["labelled"] = idb.kcat_KM.notna(); idb["n_mut_vs_HG3"] = ex.get("cmp.n_mut_vs_HG3")

    # ---- full dataset: every master column + extras + all candidate columns (labelled proteins only) + iPAE21
    full = pd.concat([idb, ent.drop(columns=["kcat_KM", "log10_kcat_KM", "family"], errors="ignore"), ex.add_prefix("x:").rename(columns=lambda c: c)], axis=1)
    full = full.loc[:, ~full.columns.duplicated()]
    mc = M.add_prefix("cand:"); full = full.join(mc)
    if not ipae.empty:
        full = full.join(ipae.add_prefix("x:").rename(columns=lambda c: c))
    full.to_csv(out / "kemp_eval_dataset_full.csv", float_format="%.6g")

    # ---- key dataset
    cand_rank = C.reindex(C.rho.abs().sort_values(ascending=False).index)
    top = list(cand_rank.candidate.head(40))
    curated = ["placer.apo.frac_pass_core", "placer_apoarm.apo.frac_pass_core", "ah.apo.cst_n_pass", "ah.pocket_bb_rmsd", "tsa_pk.geo.base_sasa", "sub.cst.penalty", "pl2.sub_holo.nac30", "pl2.tsa_holo.nac30",
               "strain.E_noncst", "apo.reorg_E", "rlx.cst.n_pass", "placer_sub.holo.lig_rmsd_to_input", "placer_sub.holo.ang_C_O_H_mean", "rlx.geo.n_donors_3p5", "rlx.field.CH.env_nobase.vac", "rlx.field.LG.env.vac",
               "ros.tsa.dG_sep", "ros.tsa.sc", "ros.tsa.packstat_int", "cmp.n_mut_vs_HG3", "cmp.pocket_frac_pos", "cmp.pocket_charge", "cmp.length"]
    ipc = [c for c in M.columns if c.startswith("ipSAE2/Boltz|")]
    keycols = list(dict.fromkeys(top + [c for c in curated if c in M.columns] + ipc))
    key = idb[["protein", "family", "kcat_KM", "log10_kcat_KM", "n_mut_vs_HG3", "pdb_ids", "kcat_KM_source"]].copy()
    for c in curated:   # curated single columns may live in master/extras rather than the candidate matrix
        src = M[c] if c in M.columns else (ent[c] if c in ent.columns else (ex[c] if c in ex.columns else None))
        if src is not None and c not in keycols:
            key[c] = src
    kk = pd.DataFrame({c: M[c] for c in keycols if c in M.columns}); key = key.join(kk, how="left")
    for c in top[:12]:   # run-1 / run-2 values of the headline tmol quantities (shows the relax noise)
        if c.startswith("tmol") or c.startswith("pr."):
            if c in M1.columns: key[f"{c}__run1"] = M1[c]
            if c in M2.columns: key[f"{c}__run2"] = M2[c]
    if not ipae.empty:
        key = key.join(ipae)
    key = key.loc[key.kcat_KM.notna()].sort_values("kcat_KM")
    key.to_csv(out / "kemp_eval_dataset_key.csv", float_format="%.6g", index=False)

    # ---- correlation table
    lab = ent[ent.log10_kcat_KM.notna()]; y = lab.log10_kcat_KM
    R = C.copy()
    R["rho_run1"] = R.candidate.map(C1.rho); R["rho_run2"] = R.candidate.map(C2.rho)
    R["rho_HG3_only"] = R.candidate.map(CH.rho); R["q_HG3_only"] = R.candidate.map(CH.q_BH)
    R["description"] = R.candidate.map(describe)
    R = R.rename(columns={"rho": "rho_all19_tmolmean", "p": "p_all19", "q_BH": "q_BH_all19"})
    extra_rows = []
    if not ipae.empty:
        for c in [c for c in ipae.columns if c.startswith("ip.") or c.startswith("ipAE21")]:
            v = ipae[c].reindex(lab.index); ok = v.notna()
            rho, p = spearmanr(v[ok], y[ok]); extra_rows.append(dict(candidate=c, set="ADDED-iPAE21", n=int(ok.sum()), rho_all19_tmolmean=rho, p_all19=p,
                                                                   description="ipSAE2 interface PAE ligand→protein (Å) " + ("(per arm)" if c.startswith("ip.") else "(arm difference)")))
    if extra_rows:
        R = pd.concat([R, pd.DataFrame(extra_rows)], ignore_index=True)
    R = R.reindex(R.rho_all19_tmolmean.abs().sort_values(ascending=False).index)
    cols = ["candidate", "description", "set", "n", "rho_all19_tmolmean", "p_all19", "q_BH_all19", "rho_run1", "rho_run2", "rho_within_family", "rho_loo_min", "rho_loo_max", "auroc_100000", "auroc_10000", "auroc_500000",
            "rho_old_labels", "partial_given_nmut", "rho_HG3_only", "q_HG3_only"]
    R = R[[c for c in cols if c in R.columns]]; R.to_csv(out / "kemp_eval_correlations_all.csv", index=False, float_format="%.4g")

    # ---- data dictionary
    dd = [dict(token=k, meaning=v) for k, v in TOK.items()] + [dict(token=f"type:{k}", meaning=v) for k, v in TYPE.items()] + [dict(token=f"tool:{k}", meaning=v) for k, v in TOOL.items()]
    pd.DataFrame(dd).to_csv(out / "kemp_eval_data_dictionary.csv", index=False)
    print("CSV ->", out, "| key", key.shape, "| full", full.shape, "| correlations", R.shape)

    # ---- PDF
    from .report_pdf import build_pdf
    retest = pd.read_csv(work / "report/retest.csv") if (work / "report/retest.csv").exists() else None
    build_pdf(out / "kemp_eval_report.pdf", dict(R=R, N=N, NH=NH, key=key, M=M, lab=lab, retest=retest, C=C, CH=CH, ipae=ipae, top=a.top, figs=out / "figures",
                                                  glossary=(TOK, TYPE, TOOL), nlab=len(lab)))
    print("PDF ->", out / "kemp_eval_report.pdf")


if __name__ == "__main__":
    main()
