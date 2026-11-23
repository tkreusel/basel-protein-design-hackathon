"""Shortlist of ~28 sensible, well-correlating metrics (one per redundancy cluster) for the user to choose from.

  python -m maxi_eval.metric_proposal            -> work/report/deliverables/kemp_eval_metric_proposal.csv
  python -m maxi_eval.metric_proposal --pick 3,5,8 --out clean.csv    -> clean per-protein CSV of the chosen proposal ids
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

# (id, recommended, tool, name, candidate column, short explanation, caveat)
PICKS = [
    ("tmol swap|int.lk_ball_iso|tsa-prod", 1, "tmol", "Solvation-type pocket energy, TSA − product", "tmol (Rosetta energy function in PyTorch) ligand–pocket 'lk_ball_iso' energy of the TSA minus that of the product placed in the same pocket (REU).", "Small magnitudes (sd ≈ 0.6 REU); needs the mean of ≥2 relax runs (run1/run2 shown)."),
    ("tmol swap|int.fa_lk|tsa-sub", 1, "tmol", "Desolvation penalty, TSA − substrate", "tmol Lazaridis–Karplus desolvation term of ligand–pocket pairs, TSA minus substrate, same pocket (REU).", "Sign is counter-intuitive (a larger penalty goes with higher activity); likely compensates lk_ball_iso."),
    ("tmol swap|int.fa_elec|sub-prod", 0, "tmol", "Electrostatics, substrate − product", "tmol Coulomb ligand–pocket energy of the substrate minus the product in the same pocket (REU).", "Correlated with the solvation terms (same 'pocket favours the charge-developing species' theme)."),
    ("tmol swap|int.total|tsa", 1, "tmol", "TSA–pocket interaction energy", "Total tmol ligand–pocket interaction energy of the TSA (REU); more negative = stronger binding of the transition-state mimic.", "Absolute value, no arm comparison needed; weaker than the differences."),
    ("tmol swap|E_total|tsa-prod", 0, "tmol", "Total complex energy, TSA − product", "tmol total energy of the TSA complex minus the product complex, same pocket (REU).", "Includes ligand-independent protein noise; sd ≈ 5.8 REU."),
    ("tmol swap|buns.n_unsat_ligand|tsa-prod", 1, "tmol", "Buried unsatisfied ligand polars, TSA − product", "Difference in the number of buried ligand polar atoms lacking an H-bond partner.", "Coarse counts (steps of 0.5); strongest inside HG3 (−0.86)."),
    ("pr.tsa_pk-sub.int_polar", 0, "tmol", "TSA–polar-residue interaction, TSA − substrate", "Per-residue tmol energies: summed ligand interaction with polar pocket residues (Ser/Thr/Asn/Gln/Tyr…), TSA minus substrate.", "Run-to-run reproducibility only moderate (−0.76 / −0.52)."),
    ("tmol swap|geo.contact_min|tsa-prod", 0, "tmol", "Closest ligand–protein contact, TSA − product", "Shortest ligand–protein heavy-atom distance (Å), TSA minus product: tightness of fit.", "Single-distance measure, sensitive to the exact pose."),
    ("raw.int.rep_over_atr", 0, "tmol", "Packing quality of the TSA (repulsion/attraction)", "Ratio of repulsive to attractive Lennard-Jones ligand–pocket energy on the unrelaxed Boltz structure (lower = better packing).", "Deterministic (no relax noise)."),
    ("rlxcst.field.CH.env_nobase.rdd4", 0, "tmol", "Electric field on the transferring C–H", "Field from the protein environment (without the base) along the C–H bond, point charges with ε = 4r, MV/cm.", "Sign is opposite to the usual hypothesis (weaker field ↔ more active); crude unscreened point charges."),
    ("sub.field.CH.env.rdd4", 0, "tmol", "Electric field on C–H, substrate complex", "Same field for the substrate pose including all environment (MV/cm).", "Right sign but weaker across all 19 (0.57)."),
    ("tmol swap|geo.n_hbond_res|tsa", 0, "tmol", "", "", ""),  # placeholder replaced below if absent
    ("bz.sub-prod.cat_to_pocket_pae", 1, "Boltz", "Pocket-trust: catalytic→pocket PAE, substrate − product", "Boltz predicted aligned error between catalytic residues and pocket residues (Å), substrate complex minus product complex.", "Confidence, not energetics; one Boltz seed per protein."),
    ("ipSAE2/Boltz|boltz.ligand_iptm|sub-prod", 1, "ipSAE2/Boltz", "Ligand ipTM, substrate − product", "Boltz interface confidence for the ligand, substrate complex minus product complex.", "Weaker in HG3 (0.42); one seed."),
    ("ipSAE2/Boltz|ipsae.p2l|sub-prod", 0, "ipSAE2/Boltz", "ipSAE protein→ligand, substrate − product", "ipSAE2 interface score (protein→ligand, d0 = 2 Å), substrate minus product.", "Correlated with ligand ipTM."),
    ("ipSAE2/Boltz|site.pae_cat_to_lig|sub-prod", 0, "ipSAE2/Boltz", "Catalytic→ligand PAE, substrate − product", "Boltz PAE between catalytic residues and the ligand (Å), substrate minus product.", "Confidence measure."),
    ("ipSAE2/Boltz|ipsae.p2l|sub", 0, "ipSAE2/Boltz", "ipSAE protein→substrate (single arm)", "ipSAE2 score of the substrate complex alone.", "Weak in HG3 (0.20); a simple single-arm option."),
    ("pl2.sub_holo.nac30", 1, "PLACER", "Substrate near-attack fraction", "Fraction of the 20 PLACER models with base O···H3 ≤ 3.0 Å and angle ≥ 130° (substrate-bound ensemble).", "20 samples → ±0.1 binomial noise."),
    ("placer_sub.holo.ang_C_O_H_mean", 1, "PLACER", "Substrate C–O···H angle", "Mean C–O···H angle (°) of the base carboxylate towards the substrate's H3 across the PLACER ensemble.", "Narrow range (141–155° in HG3)."),
    ("placer_sub.holo.lig_rmsd_to_input", 0, "PLACER", "Substrate pose stability", "RMSD (Å) of the PLACER-predicted substrate pose to the Boltz input pose (low = pose agrees).", "Indirect measure of pose quality."),
    ("pl2.tsa_holo.nac30", 0, "PLACER", "TSA near-attack fraction", "Same near-attack fraction for the TSA-bound PLACER ensemble.", "Weaker than for the substrate (0.51)."),
    ("tsa_pk.geo.base_sasa", 1, "geometry", "Catalytic-base burial (SASA)", "Solvent-accessible surface (Å²) of the catalytic carboxylate oxygen in the TSA complex; low = buried.", "Strong only inside HG3 (−0.80); weak across families (−0.30)."),
    ("sub.geo.d_O_H", 1, "geometry", "Base···H3 distance, substrate", "Distance (Å) from the catalytic base O to the transferable H3 in the relaxed substrate complex.", "HG3-style site; −0.80 in HG3, −0.58 overall."),
    ("tmol swap|geo.lig_rel_sasa|prod", 0, "geometry", "Product exposure (relative SASA)", "Relative solvent-accessible surface of the product in the pocket (low = buried).", "Mostly a between-family effect (HG3-only −0.29)."),
    ("rlx.geo.n_donors_3p5", 0, "geometry", "Oxyanion donors (Gln50-type)", "Number of H-bond donors within 3.5 Å of the TSA oxyanion O.", "Weak (0.45); mechanistically sensible."),
    ("pr.tsa_pk.n_hbond_res", 0, "geometry", "Residues H-bonding the TSA", "Number of pocket residues that make an H-bond to the TSA (tmol per-residue energies).", "Moderate (0.58)."),
    ("cmp.pocket_frac_pos", 0, "composition", "Positive residues in pocket", "Fraction of Arg/Lys/His among pocket residues.", "Sequence-only; negative ρ."),
    ("cmp.n_mut_vs_HG3", 0, "baseline", "Mutations relative to HG3", "Number of mutations versus HG3 — a baseline for 'how evolved', not a mechanism.", "Defined for HG3-lineage only; use as a baseline to beat."),
    ("placer.apo.frac_pass_core", 0, "preorganisation", "Pre-formed apo fraction (PLACER)", "Fraction of the apo ensemble already in the catalytic geometry (the original preorganisation metric).", "Weak (0.35) and HG3-specific windows."),
]
PICKS = [p for p in PICKS if p[3]]


def load():
    work = Path("work")
    M = pd.read_csv(work / "report/everything_tmolmean/candidate_matrix.csv", index_col=0)
    R = pd.read_csv("work/report/deliverables/kemp_eval_correlations_all.csv").set_index("candidate")
    ex = pd.read_csv(work / "features_extra_mean.csv", index_col=0)
    from .analyze_full import load_master
    ent = load_master(Path("work/report_tmolmean"), work, "config/labels_manual_dataset.csv")
    return M, R, ex, ent


def values(c, M, ex, ent):
    for s in (M, ex, ent):
        if c in s.columns:
            return s[c].astype(float)
    return pd.Series(dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pick", default=None, help="comma-separated proposal numbers -> clean per-protein CSV")
    ap.add_argument("--out", default="work/report/deliverables/kemp_eval_metric_proposal.csv")
    a = ap.parse_args()
    M, R, ex, ent = load()
    lab = ent[ent.log10_kcat_KM.notna()]; y = lab.log10_kcat_KM
    hg = lab.family == "HG3"
    rows = []
    for i, (cid, rec, tool, name, expl, cav) in enumerate(PICKS, 1):
        v = values(cid, M, ex, ent).reindex(lab.index)
        ok = v.notna()
        r19, p19 = spearmanr(v[ok], y[ok]) if ok.sum() > 5 else (np.nan, np.nan)
        okh = ok & hg
        rh = spearmanr(v[okh], y[okh])[0] if okh.sum() > 5 and v[okh].nunique() > 1 else np.nan
        rr = R.loc[cid] if (cid in R.index and tool in ('tmol', 'geometry')) else None   # run-1/run-2 only differ for tmol-relax-derived values
        rows.append(dict(id=i, recommended=bool(rec), tool=tool, metric=name, what_it_measures=expl,
                         higher_value_means=("higher activity" if r19 > 0 else "lower activity"), rho_all19=r19, p_all19=p19, n19=int(ok.sum()),
                         rho_HG3_only=rh, n_HG3=int(okh.sum()), rho_tmol_run1=(rr.rho_run1 if rr is not None else np.nan), rho_tmol_run2=(rr.rho_run2 if rr is not None else np.nan),
                         rho_leave_one_out_min=(rr.rho_loo_min if rr is not None else np.nan), rho_leave_one_out_max=(rr.rho_loo_max if rr is not None else np.nan),
                         value_min=v.min(), value_max=v.max(), caveat=cav, column_id=cid))
    T = pd.DataFrame(rows)
    if a.pick:
        ids = [int(x) for x in a.pick.split(",")]
        clean = lab[["family", "kcat_KM", "log10_kcat_KM"]].copy()
        for i in ids:
            r = T[T.id == i].iloc[0]; clean[f"{i:02d}_{r.metric}"] = values(r.column_id, M, ex, ent).reindex(clean.index)
        clean.insert(0, "protein", clean.index); clean.sort_values("kcat_KM").to_csv(a.out, index=False, float_format="%.5g"); print("clean ->", a.out); return
    T.to_csv(a.out, index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    print(T[["id", "recommended", "tool", "metric", "rho_all19", "rho_HG3_only", "rho_tmol_run1", "rho_tmol_run2", "n19", "n_HG3"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
