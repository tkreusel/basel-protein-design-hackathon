"""Aggregate ALL round-4/5 design metrics (all seeds, incl. channel arm) -> runs/20261121_r5_SUMMARY/{design_metrics_all.csv, arm_summary.csv, channel_arm_vs_others.csv}
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I eval_r5.py. Core criteria as select_passers.py; channel_ok = bottleneck >= 1.4 A; pass_channel = core + channel_ok."""
import glob, os, numpy as np, pandas as pd
H = os.path.expanduser("~/HACKATHON_ETH"); S = H + "/runs/20261121_r5_SUMMARY"; os.makedirs(S, exist_ok=True)
d = pd.concat([pd.read_csv(f).assign(file=os.path.basename(f)) for f in sorted(glob.glob(H + "/runs/20261121_r4_*/metrics/*.csv"))], ignore_index=True)
cc = pd.read_csv(H + "/runs/20261121_r4_SUMMARY/channel_check.csv").set_index("id").bottleneck_radius        # r4 designs (0.5 A grid)
if "bottleneck_radius" not in d: d["bottleneck_radius"] = np.nan
d["bottleneck_radius"] = d.bottleneck_radius.fillna(d.design.map(cc))
d["seed"] = d.get("seed", pd.Series(7, index=d.index)).fillna(7).astype(int); d["motif"] = d.campaign.str.replace("20261121_r4_", "")
d["rep"] = d.design.str.extract(r"_(\d+)_beam")[0]; d["lineage"] = d.motif + "|" + d.arm + "|s" + d.seed.astype(str) + "|r" + d.rep
F = {"lig_sasa": d.ligand_sasa_frac <= 0.10, "base_buried": d.base_sasa_rel <= 0.10, "no_sc_polar": d.base_sc_polar_4p5 <= np.where(d.motif.str.contains("with-trp"), 1, 0), "no_cation": d.base_cation_6 == 0,
     "pocket": d.volume_ratio.between(0.4, 1.4), "nitro": (d.nitro_hb_far >= 1) & (d.nitro_hb_near_base == 0), "clash": (d.lig_clash_atoms == 0) & (d.pp_clash_frac <= 0.005),
     "motif": (d.motif_base_rmsd <= 0.3) & (d.motif_partner_rmsd <= 0.3), "tmol": d.tmol_ddg < 0}
for k, v in F.items(): d["ok_" + k] = v
d["pass_core"] = np.logical_and.reduce(list(F.values())); d["channel_ok"] = d.bottleneck_radius >= 1.4; d["pass_channel"] = d.pass_core & d.channel_ok
d.to_csv(S + "/design_metrics_all.csv", index=False)
g = d.groupby(["motif", "arm"])
cols = ["length", "helix", "strand", "bottleneck_radius", "ligand_sasa_frac", "base_sasa_rel", "base_sc_polar_4p5", "base_cation_6", "nitro_hb_far", "volume_ratio", "motif_base_rmsd", "motif_partner_rmsd", "tmol_ddg"]
sm = g[cols].mean().round(2).join(g.agg(n=("design", "size"), lineages=("lineage", "nunique"), seeds=("seed", "nunique"), pass_core=("pass_core", "sum"), channel_ok=("channel_ok", "sum"), pass_channel=("pass_channel", "sum")))
sm["lineages_pass_core"] = d[d.pass_core].groupby(["motif", "arm"]).lineage.nunique(); sm.to_csv(S + "/arm_summary.csv")
pd.set_option("display.width", 300); pd.set_option("display.max_columns", 60); print(sm.to_string())
