"""Aggregate round-4 metrics (runs/20261121_r4_*/metrics/*.csv) -> runs/20261121_r4_SUMMARY/{design_metrics_all.csv, arm_summary.csv, crystal_vs_designs.csv}
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I eval_r4.py
Pass criteria (crystal-calibrated, ../crystals/crystal_metrics.csv): ligand SASA<=0.05; base SASA_rel<=0.10; side-chain polar atoms <=4.5 A of base ==0; Lys/Arg <=6 A ==0;
pocket ratio 0.4-1.4 (snug flag <=0.7); >=1 nitro O H-bonded to a donor >=6 A from the base and none near the base; no clashes (ligand 2.4 A, protein-protein 0.5 %);
motif rigid: base and partner RMSD <=0.3 A (ligand-superposed); TMol dG<0."""
import glob, os, numpy as np, pandas as pd
H = os.path.expanduser("~/HACKATHON_ETH"); S = H + "/runs/20261121_r4_SUMMARY"; os.makedirs(S, exist_ok=True)
d = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(H + "/runs/20261121_r4_*/metrics/*.csv"))], ignore_index=True)
d["lineage"] = d.design.str.extract(r"_(\d+)_beam")[0]
d["motif"] = d.campaign.str.replace("20261121_r4_", "")
F = {"lig_sasa": d.ligand_sasa_frac <= 0.05, "base_buried": d.base_sasa_rel <= 0.10, "base_no_sc_polar": d.base_sc_polar_4p5 <= np.where(d.motif.str.contains("with-trp"), 1, 0),  # with-trp: the fixed motif Trp44 NE1 is always within 4.5 A of the base "base_no_cation": d.base_cation_6 == 0,
     "pocket": d.volume_ratio.between(0.4, 1.4), "nitro_hb": (d.nitro_hb_far >= 1) & (d.nitro_hb_near_base == 0), "clash": (d.lig_clash_atoms == 0) & (d.pp_clash_frac <= 0.005),
     "motif_rigid": (d.motif_base_rmsd <= 0.3) & (d.motif_partner_rmsd <= 0.3), "tmol": d.tmol_ddg < 0}
for k, v in F.items(): d["ok_" + k] = v
d["pass_all"] = np.logical_and.reduce(list(F.values())); d["snug"] = d.volume_ratio <= 0.7
d.to_csv(S + "/design_metrics_all.csv", index=False)
cols = ["length", "helix", "strand", "ligand_sasa_frac", "base_sasa_rel", "base_apolar_frac", "base_sc_polar_4p5", "base_bb_polar_4p5", "base_cation_6", "nitro_hb_far", "nitro_hb_near_base",
        "volume_ratio", "lig_clash_atoms", "motif_base_rmsd", "motif_partner_rmsd", "tmol_ddg"]
g = d.groupby(["motif", "arm"]); sm = g[cols].mean().round(2).join(g.agg(n=("design", "size"), **{"ok_" + k: ("ok_" + k, "sum") for k in F}, pass_all=("pass_all", "sum")))
sm["n_lineages_pass"] = d[d.pass_all].groupby(["motif", "arm"]).lineage.nunique(); sm.to_csv(S + "/arm_summary.csv"); pd.set_option("display.width", 280); pd.set_option("display.max_columns", 50); print(sm.T.to_string())
c = pd.read_csv(H + "/crystals/crystal_metrics.csv").groupby("pdb").mean(numeric_only=True)
cc = ["volume_ratio", "ligand_sasa_frac", "base_sasa_rel", "base_apolar_frac", "base_sc_polar_4p5", "base_bb_polar_4p5", "base_cation_6", "nitro_hb_far", "nitro_hb_near_base"]
tab = pd.concat([c[cc].T.add_prefix("xtal_"), d[d.pass_all | True].groupby("arm")[cc].mean().T.add_prefix("design_")], axis=1); tab.round(2).to_csv(S + "/crystal_vs_designs.csv")
