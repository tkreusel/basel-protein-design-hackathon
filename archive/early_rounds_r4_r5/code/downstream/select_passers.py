"""Select designs for LigandMPNN + Boltz-2: all steered-arm designs passing the core criteria (r4 + r5 metrics) -> runs/20261121_r5_downstream/selection.csv + designs_in/*.cif
Core criteria (same as eval_r4 but ligand SASA <=0.10 so the channel arm is not penalised; motif rigidity <=0.3 A; TMol dG<0; with-Trp allows 1 polar atom = fixed Trp44 NE1).
Channel arm: flag channel_ok = bottleneck >= 1.4 A (reported, not required for selection).
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I select_passers.py"""
import glob, os, shutil, numpy as np, pandas as pd
H = os.path.expanduser("~/HACKATHON_ETH"); D = H + "/runs/20261121_r5_downstream"; os.makedirs(D + "/designs_in", exist_ok=True)
d = pd.concat([pd.read_csv(f).assign(src=f) for f in sorted(glob.glob(H + "/runs/20261121_r4_*/metrics/*.csv"))], ignore_index=True)
if "seed" not in d: d["seed"] = 7
d["seed"] = d.seed.fillna(7).astype(int); d["motif"] = d.campaign.str.replace("20261121_r4_", "")
if "bottleneck_radius" not in d: d["bottleneck_radius"] = np.nan
F = {"lig_sasa": d.ligand_sasa_frac <= 0.10, "base_buried": d.base_sasa_rel <= 0.10, "no_sc_polar": d.base_sc_polar_4p5 <= np.where(d.motif.str.contains("with-trp"), 1, 0),
     "no_cation": d.base_cation_6 == 0, "pocket": d.volume_ratio.between(0.4, 1.4), "nitro": (d.nitro_hb_far >= 1) & (d.nitro_hb_near_base == 0),
     "clash": (d.lig_clash_atoms == 0) & (d.pp_clash_frac <= 0.005), "motif": (d.motif_base_rmsd <= 0.3) & (d.motif_partner_rmsd <= 0.3), "tmol": d.tmol_ddg < 0}
d["pass_core"] = np.logical_and.reduce(list(F.values())); d["channel_ok"] = d.bottleneck_radius >= 1.4
sel = d[d.pass_core & (d.arm != "unsteered")].copy()
def path(r):
    base = f"{H}/runs/{r.campaign}/structures/{r.arm}"
    return glob.glob(f"{base}/**/designs/{r.design}.cif", recursive=True)[0]
sel["cif"] = sel.apply(path, axis=1)
for _, r in sel.iterrows(): shutil.copy(r.cif, f"{D}/designs_in/{r.design}.cif")
sel.drop(columns=["src"]).to_csv(D + "/selection.csv", index=False); d.drop(columns=["src"]).to_csv(D + "/all_designs_core_flags.csv", index=False)
print(sel.groupby(["motif", "arm"]).agg(n=("design", "size"), channel_ok=("channel_ok", "sum")).to_string()); print("total", len(sel))
