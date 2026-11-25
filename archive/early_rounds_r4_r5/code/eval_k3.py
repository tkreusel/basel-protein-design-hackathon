"""Round-3 evaluation (k3*/k3u_* arms) with crystal-calibrated criteria -> out/20261121_eval_k3.csv
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I eval_k3.py ; TMol dG from out/20261121_eval_rfd4_tmol_k3.csv (EVAL_GLOB='k3*' EVAL_TAG=_k3 evaltmol job) if present.
Criteria (crystal reference in ../crystals/crystal_metrics.csv: 29SB, 8RD5, 4BS0, 7K4Q, 5RGF): ligand SASA<=0.05; base SASA_rel<=0.10; side-chain polar atoms within 4.5 A of base == 0;
Lys/Arg within 6 A of base == 0; pocket volume ratio 0.4-1.4 (crystal range) [snug: <=0.75]; >=1 nitro O H-bonded to a donor >=6 A from the base; no donor-near-base nitro bonds;
no clashes (ligand 2.4 A, protein-protein 0.5 %); motif atoms kept (RMSD<=0.5 A on ref-superposition)."""
import glob, os, numpy as np, pandas as pd, biotite.structure as s, biotite.structure.io.pdbx as pdbx
from rfproteina.metrics.custom_metrics.kemp_base import KempBaseMetric
from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
from rfproteina.metrics.custom_metrics.kemp_chem import KempChemMetric
H = os.path.expanduser("~/HACKATHON_ETH"); R = H + "/rfd4_runs"
mb, mp = KempBaseMetric(), KempPocketMetric(target_ratio=0.6, tol=0.15)
mc_g = KempChemMetric(ref_pdb=H + "/kemp_theozyme_GDM_KE_1483.pdb")
mc_h = KempChemMetric(ref_pdb=H + "/HG3_17_kemp_motif.pdb", base_res="ASP", base_match=("CB", "CG", "OD1", "OD2"), base_o=("OD1", "OD2"), hb_res="GLN", hb_match=("CG", "CD", "OE1", "NE2"), hb_atom="NE2")
rows = []
for f in sorted(glob.glob(f"{R}/out/20261121_k3*/inference_samples/*/*.cif")):
    arm = f.split("/out/20261121_")[1].split("/")[0]; hg = "hg" in arm
    a = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False); a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]
    c = (mc_h if hg else mc_g).compute(a)
    r = {**mp.compute(a), **mb.compute(a), "motif_base_rmsd": c["base_rmsd"], "motif_partner_rmsd": c["hb_rmsd"], "base_hb_partners_any": c["base_hbond_partners"]}
    sse = np.array(s.annotate_sse(a[~a.hetero])); r.update(arm=arm, design=os.path.basename(f), length=len(sse), helix=(sse == "a").mean(), strand=(sse == "b").mean())
    rows.append(r)
d = pd.DataFrame(rows); t = f"{R}/out/20261121_eval_rfd4_tmol_k3.csv"
d = d.merge(pd.read_csv(t), on=["arm", "design"], how="left") if os.path.exists(t) else d.assign(tmol_ddg=np.nan)
F = {"lig_sasa": d.ligand_sasa_frac <= 0.05, "base_buried": d.base_sasa_rel <= 0.10, "base_no_sc_polar": d.base_sc_polar_4p5 == 0, "base_no_cation": d.base_cation_6 == 0,
     "pocket": d.volume_ratio.between(0.4, 1.4), "nitro_hb": (d.nitro_hb_far >= 1) & (d.nitro_hb_near_base == 0), "clash": (d.lig_clash_atoms == 0) & (d.pp_clash_frac <= 0.005),
     "motif": (d.motif_base_rmsd <= 0.5) & (d.motif_partner_rmsd <= 0.5)}
for k, v in F.items(): d["ok_" + k] = v
d["pass_all"] = np.logical_and.reduce(list(F.values())); d.to_csv(f"{R}/out/20261121_eval_k3.csv", index=False)
pd.set_option("display.width", 250)
cols = ["length", "helix", "strand", "ligand_sasa_frac", "base_sasa_rel", "base_apolar_frac", "base_sc_polar_4p5", "base_bb_polar_4p5", "base_cation_6", "nitro_hb_far", "nitro_hb_any", "nitro_hb_near_base",
        "volume_ratio", "lig_clash_atoms", "pp_clash_frac", "motif_base_rmsd", "motif_partner_rmsd", "tmol_ddg"]
print(d.groupby("arm")[cols].mean().round(2).T.to_string()); print(d.groupby("arm").agg(n=("design", "size"), **{"ok_" + k: ("ok_" + k, "sum") for k in F}, pass_all=("pass_all", "sum")).to_string())
