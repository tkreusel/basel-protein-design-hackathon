"""Score Maxi's overnight top-30 (10 Boltz TSA samples each) with the same stand-in as scoring/rank_v2.py, for a like-for-like
comparison with our designs. iPDE and ipSAE2 are taken from his table (no confidence json beside his CIFs)."""
import csv, glob, sys, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scoring")); sys.path.insert(0, str(ROOT / "src"))
import host_score as H
from host_geometry import geometry, burial_count
from kemp import geometry as g
from biotite.structure.io.pdbx import CIFFile, get_structure
T = Path("/mnt/work/projects/Maxi/kemp_results_20261122/final_strict_all_gates_ipsae0.8/top30")
tab = {int(r["rank30"]): r for r in csv.DictReader(open(T / "top30_table.csv"))}
fa = {}
cur = None
for l in open(T / "top30.fasta"):
    if l.startswith(">"): cur = int(l[1:3]); fa[cur] = ["", l[1:].split()[0]]
    else: fa[cur][0] += l.strip()
def gscore(x, k, sf=0.6):
    z = (x - H.REF[k]["mean"]) / H.REF[k]["sigma"]
    if k.endswith("dist") and z < 0: z = z / sf
    return float(np.exp(-0.5 * z * z))
out = []
for rk in range(1, 31):
    d = next(T.glob(f"{rk:02d}_*"))
    seq = fa[rk][0]
    cifs = sorted((d / "boltz_TSA_10_samples").glob("*.cif"))
    S = [get_structure(CIFFile.read(c), model=1) for c in cifs]
    s0 = S[0]; prot = s0[(~s0.hetero) & (s0.chain_id == "A")]; lig = s0[s0.res_name == "6NT"]
    n3 = lig.coord[lig.atom_name == "N3"][0]; n1 = lig.coord[lig.atom_name == "N1"][0]
    base = min(((np.linalg.norm(s0.coord[(s0.res_id == r) & np.isin(s0.atom_name, ("OD1", "OD2")) & (s0.chain_id == "A")] - n3, axis=1).min(), r)
                for r in set(prot.res_id[prot.res_name == "ASP"])))[1]
    gl = [(np.linalg.norm(s0.coord[(s0.res_id == r) & (s0.atom_name == "NE2") & (s0.chain_id == "A")][0] - n1), r) for r in set(prot.res_id[prot.res_name == "GLN"])]
    donor = min(gl)[1]
    amap = g.load_atom_map(ROOT / "configs/atom_map.yaml")
    G, ligs, cas, bur = [], [], [], []
    for s in S:
        lc, lr = g.iter_ligand_copies(s, amap)[0]
        m = geometry(s, chain="A", base_resid=int(base), donor_resid=int(donor), lig_chain=lc, lig_resid=lr)
        contact = gscore(m["b_dist"], "b_dist") * gscore(m["d_dist"], "d_dist")
        ba = min(gscore(abs(m["b_ang_NH"]), "b_ang_NH"), gscore(abs(m["b_oop"]), "b_oop"))
        da = min(gscore(abs(m["d_ang_ax"]), "d_ang_ax"), gscore(abs(m["d_oop"]), "d_oop"))
        G.append((contact, ba * da, contact * ba * da, m["b_dist"], m["d_dist"]))
        bur.append(burial_count(s, chain="A", base_resid=int(base))["burial_CG"])
        ligs.append(s[s.res_name == "6NT"].coord); cas.append(s[(~s.hetero) & (s.chain_id == "A") & (s.atom_name == "CA")].coord)
    pose = []
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            R, pc, qc = H.kabsch(cas[j], cas[i]); pose.append(float(np.sqrt((((ligs[j] - pc) @ R + qc - ligs[i]) ** 2).sum(1).mean())))
    cons = float(np.exp(-np.median(pose)))             # simplified (pose term only)
    contact, appr, geom, bd, dd = (float(np.mean([x[i] for x in G])) for i in range(5))
    ipde = float(tab[rk]["boltz_ipde.tsa"]); s_ipde = float(np.clip((4 - ipde) / 3.5, 0, 1)); s_bur = min(1.0, float(np.mean(bur)) / 45.9)
    tot = geom ** 2 * s_ipde * s_bur * cons
    out.append(dict(rank30=rk, name=fa[rk][1], total_mult=round(tot, 4), geometry=round(geom, 3), contact=round(contact, 3), approach=round(appr, 3),
                    s_ipde=round(s_ipde, 3), s_burial=round(s_bur, 3), cons_pose_only=round(cons, 3), b_dist=round(bd, 2), d_dist=round(dd, 2),
                    ipde=ipde, ipsae=tab[rk]["ipsae.tsa"], burial=round(float(np.mean(bur)), 1), sequence=seq, base=int(base), donor=int(donor), length=len(seq)))
out.sort(key=lambda r: -r["total_mult"])
with open(ROOT / "submission/maxi_top30_scored.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for r in out: print(f"{r['rank30']:>2} {r['name'][-46:]:46} tot={r['total_mult']:.3f} geom={r['geometry']:.2f} (c {r['contact']:.2f} a {r['approach']:.2f}) ipde={r['ipde']} cons={r['cons_pose_only']} b={r['b_dist']} d={r['d_dist']} len={r['length']}")
