"""Evaluate the theozyme (GDM_KE_1483) runs: pocket + chemistry metrics + clashes + filter flags -> out/20261121_eval_k2.csv
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I eval_k2.py   (TMol dG from out/20261121_eval_rfd4_tmol.csv if present)
Filters: ligand clash 0 (<2.4 A), pp clash <=0.5 %, base_rmsd<=0.5 & hb_rmsd<=0.5 (motif kept), volume ratio 0.8-1.3, ligand SASA<=0.15,
base_sasa_rel<=0.15, base_hbond_partners<=1, tmol_ddg<0, nitro_unsat_buried==0."""
import glob, os, numpy as np, pandas as pd
import biotite.structure.io.pdbx as pdbx
from scipy.spatial import cKDTree
from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
from rfproteina.metrics.custom_metrics.kemp_chem import KempChemMetric
H = os.path.expanduser("~/HACKATHON_ETH"); R = H + "/rfd4_runs"
mp, mc = KempPocketMetric(), KempChemMetric(ref_pdb=H + "/kemp_theozyme_GDM_KE_1483.pdb"); rows = []
for f in sorted(glob.glob(f"{R}/out/20261121_k2*/inference_samples/*/*.cif")):
    arm = f.split("/out/20261121_")[1].split("/")[0]
    a = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False); a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]
    p, l = a[~a.hetero], a[a.hetero]; pt = cKDTree(p.coord)
    lc = int(sum(len(x) for x in pt.query_ball_point(l.coord, 2.4))); pr = pt.query_pairs(2.2, output_type="ndarray")
    pc = float((np.abs(p.res_id[pr[:, 0]] - p.res_id[pr[:, 1]]) > 1).sum()) / len(p) if len(pr) else 0.0
    r = {**mp.compute(a), **mc.compute(a)}; r.update(arm=arm, design=os.path.basename(f), length=len(set(p.res_id)), lig_clash_atoms=lc, pp_clash_frac=pc); rows.append(r)
d = pd.DataFrame(rows)
t = f"{R}/out/20261121_eval_rfd4_tmol.csv"
if os.path.exists(t): d = d.merge(pd.read_csv(t).rename(columns={"tmol_ddg": "tmol_ddg"}), on=["arm", "design"], how="left")
else: d["tmol_ddg"] = np.nan
F = {"clash": (d.lig_clash_atoms == 0) & (d.pp_clash_frac <= 0.005), "motif": (d.base_rmsd <= 0.5) & (d.hb_rmsd <= 0.5),
     "pocket": d.volume_ratio.between(0.8, 1.3) & (d.ligand_sasa_frac <= 0.15), "base_buried": d.base_sasa_rel <= 0.15,
     "base_hb": d.base_hbond_partners <= 1, "nitro_ok": d.nitro_unsat_buried == 0, "tmol": (d.tmol_ddg < 0) | d.tmol_ddg.isna()}
for k, v in F.items(): d["ok_" + k] = v
d["pass_all"] = np.logical_and.reduce(list(F.values())); d.to_csv(f"{R}/out/20261121_eval_k2.csv", index=False)
pd.set_option("display.width", 250)
cols = ["length", "volume_ratio", "ligand_sasa_frac", "loop_frac_shell", "base_sasa_rel", "base_apolar_frac", "base_n_polar_env", "base_hbond_partners", "lig_hbonds", "nitro_hbonds", "nitro_unsat_buried", "hb_partner_dist", "base_rmsd", "hb_rmsd", "tmol_ddg"]
print(d.groupby("arm")[cols].mean().round(2).T.to_string()); print(d.groupby("arm").agg(n=("design", "size"), **{"ok_" + k: ("ok_" + k, "sum") for k in F}, pass_all=("pass_all", "sum")).to_string())
