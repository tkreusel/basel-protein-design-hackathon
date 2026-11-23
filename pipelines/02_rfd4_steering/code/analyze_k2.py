"""Detailed look at the theozyme runs: per-design filter failures, base polar neighbours, oxyanion (N1) environment, nitro O environment, fold composition.
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I analyze_k2.py  (reads out/20261121_eval_k2.csv and the CIFs)"""
import glob, os, numpy as np, pandas as pd, biotite.structure as s, biotite.structure.io.pdbx as pdbx
from scipy.spatial import cKDTree
H = os.path.expanduser("~/HACKATHON_ETH"); R = H + "/rfd4_runs"
from rfproteina.metrics.custom_metrics.kemp_chem import KempChemMetric
mc = KempChemMetric(ref_pdb=H + "/kemp_theozyme_GDM_KE_1483.pdb"); ref = mc.ref
rl = ref[ref.hetero]; rows = []; pol_rows = []
for f in sorted(glob.glob(f"{R}/out/20261121_k2*/inference_samples/*/*.cif")):
    arm = f.split("/out/20261121_")[1].split("/")[0]; des = os.path.basename(f)
    a = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False); a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]
    lig, p = a[a.hetero], a[~a.hetero]
    names = [n for n in rl.atom_name if n in set(lig.atom_name)]
    A = np.array([lig.coord[lig.atom_name == n][0] for n in names]); B = np.array([rl.coord[rl.atom_name == n][0] for n in names])
    _, tr = s.superimpose(A, B)
    rb, _ = mc._locate(p, tr, "GLU", ("CG", "CD", "OE1", "OE2")); rh, _ = mc._locate(p, tr, "SER", ("CA", "CB", "OG"))
    base = p[(p.res_id == rb) & (p.res_name == "GLU")]; bo = base[np.isin(base.atom_name, ("OE1", "OE2"))]
    ser = p[(p.res_id == rh) & (p.res_name == "SER") & (p.atom_name == "OG")]
    other = p[(p.res_id != rb) & (p.res_id != rh)]; pol = other[np.isin(other.element, ("N", "O"))]; pt = cKDTree(pol.coord)
    # polar protein atoms near the base O (<=4.5 A), excluding motif Ser
    near = sorted({(pol.res_name[i] + str(pol.res_id[i]) + ":" + pol.atom_name[i]) for c in bo.coord for i in pt.query_ball_point(c, 4.5)})
    n35 = sorted({(pol.res_name[i] + ":" + pol.atom_name[i]) for c in bo.coord for i in pt.query_ball_point(c, 3.5)})
    N1 = lig.coord[lig.atom_name == "N1"][0]; nit = lig[np.isin(lig.atom_name, ("O21", "O11"))]
    donors_N1 = [(pol.res_name[i] + ":" + pol.atom_name[i]) for i in pt.query_ball_point(N1, 3.5)]
    full = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False); sa = s.sasa(a, vdw_radii="ProtOr", point_number=200)
    nit_sa = sa[a.hetero][np.isin(lig.atom_name, ("O21", "O11"))]
    nd = [float(pt.query(c)[0]) for c in nit.coord]
    # fold composition
    sse = np.array(s.annotate_sse(p)); nstr = int(((sse[1:] == "b") & (sse[:-1] != "b")).sum() + (sse[0] == "b")); nhel = int(((sse[1:] == "a") & (sse[:-1] != "a")).sum() + (sse[0] == "a"))
    rows.append(dict(arm=arm, design=des, length=len(sse), helix=(sse == "a").mean(), strand=(sse == "b").mean(), coil=(sse == "c").mean(), n_strand_segments=nstr, n_helix_segments=nhel,
                     serOG_N1=float(np.linalg.norm(ser.coord[0] - N1)), n_donor_acceptors_N1_3p5=len(donors_N1), N1_partners=";".join(donors_N1),
                     base_polar_3p5=len(n35), base_polar_4p5=len(near), base_polar_list_3p5=";".join(n35),
                     nitro_nearest_polar=";".join(f"{d:.1f}" for d in nd), nitro_sasa=";".join(f"{x:.1f}" for x in nit_sa)))
d = pd.DataFrame(rows).merge(pd.read_csv(f"{R}/out/20261121_eval_k2.csv"), on=["arm", "design"], suffixes=("", "_e"))
d.to_csv(f"{R}/out/20261121_analysis_k2.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(d.groupby("arm")[["length", "helix", "strand", "coil", "n_strand_segments", "n_helix_segments"]].mean().round(2).to_string())
print(d[["arm", "design", "length", "pocket_volume", "ligand_volume", "volume_ratio", "ligand_sasa_frac", "n_shell", "loop_frac_shell", "tmol_ddg"]].round(2).to_string())
print(d[["arm", "design", "serOG_N1", "n_donor_acceptors_N1_3p5", "N1_partners", "base_polar_3p5", "base_polar_list_3p5", "nitro_nearest_polar", "nitro_sasa"]].round(2).to_string())
print(d.groupby("arm")[["serOG_N1", "n_donor_acceptors_N1_3p5", "base_polar_3p5", "base_polar_4p5"]].mean().round(2))
F = [c for c in d.columns if c.startswith("ok_")]
print(d.groupby("arm")[F].sum())
