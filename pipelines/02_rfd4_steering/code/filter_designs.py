"""Early filter for RFD3/RFD4 Kemp designs (no refolding yet). Writes filtered/20261121_design_metrics.csv + pass flags.
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I filter_designs.py
Inputs: rfd4_runs/out/20261121_{unsteered,pocket,pocket_loops,tmol}/inference_samples/*/*.cif ; rfd3/out/20261121_*/raw.tar.gz|cif.gz;
        rfd4_runs/out/20261121_eval_rfd4_tmol.csv (TMol ddG, RFD4 only); motif reference HG3_17_kemp_motif.pdb.
Filters (thresholds in THR): length 150-260; ligand-protein clashes 0 (heavy atoms < 2.4 A (the crystal motif's closest contact is shown in the log)); protein-protein clash fraction <= 0.5 %
(non-adjacent residues, heavy atoms < 2.2 A); motif tip RMSD <= 0.5 A after ligand superposition (best same-name residue per motif residue);
pocket volume ratio 0.8-1.3; ligand SASA fraction <= 0.15; shell residues >= 15; TMol ddG < 0 where available."""
import glob, gzip, io, os, tarfile, numpy as np, pandas as pd
import biotite.structure as s, biotite.structure.io.pdbx as pdbx, biotite.structure.io.pdb as pdb
from scipy.spatial import cKDTree
from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
H = os.path.expanduser("~/HACKATHON_ETH"); m = KempPocketMetric()
THR = dict(len=(150, 260), lig_clash_A=2.4, pp_clash_A=2.2, pp_clash_frac=0.005, motif_rmsd=0.5, ratio=(0.8, 1.3), sasa=0.15, n_shell=15)
ref = pdb.PDBFile.read(f"{H}/HG3_17_kemp_motif.pdb").get_structure(model=1)
TIPS = {"TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"), "GLN": ("CG", "CD", "OE1", "NE2"), "ASP": ("CB", "CG", "OD1", "OD2")}
rlig = ref[ref.hetero]; rtips = {r: ref[(ref.res_name == r) & ~ref.hetero & np.isin(ref.atom_name, TIPS[r])] for r in TIPS}
def motif_rmsd(a, lig):
    # superpose on ligand (matched by atom name), then best same-name residue per motif residue
    n = [x for x in rlig.atom_name if x in set(lig.atom_name)]
    A = np.array([lig.coord[lig.atom_name == x][0] for x in n]); B = np.array([rlig.coord[rlig.atom_name == x][0] for x in n])
    _, tr = s.superimpose(A, B); out = {}
    for r, tip in rtips.items():
        best = np.inf
        for rid in np.unique(a.res_id[(a.res_name == r) & ~a.hetero]):
            res = a[(a.res_id == rid) & (a.res_name == r) & ~a.hetero]
            try: X = np.array([res.coord[res.atom_name == x][0] for x in tip.atom_name])
            except IndexError: continue
            Y = tr.apply(tip.coord)
            best = min(best, float(np.sqrt(((X - Y) ** 2).sum(1).mean())))
        out[r] = best
    return max(out.values()), out
def clashes(a):
    p, l = a[~a.hetero], a[a.hetero]; pt = cKDTree(p.coord)
    lig_cl = int(sum(len(x) for x in pt.query_ball_point(l.coord, THR["lig_clash_A"])))
    pairs = pt.query_pairs(THR["pp_clash_A"], output_type="ndarray"); ri = p.res_id
    ok = np.abs(ri[pairs[:, 0]] - ri[pairs[:, 1]]) > 1 if len(pairs) else np.array([], bool)
    return lig_cl, float(ok.sum()) / max(len(p), 1)
def load(raw, gz=False):
    return pdbx.get_structure(pdbx.CIFFile.read(io.StringIO(raw) if not gz else io.StringIO(gzip.decompress(raw).decode())), model=1, use_author_fields=False)
rows = []
def score(a, **meta):
    a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]; lig = a[a.hetero]
    r = m.compute(a); lc, pc = clashes(a); mr, per = motif_rmsd(a, lig)
    rows.append(dict(**meta, length=len(set(a.res_id[~a.hetero])), lig_clash_atoms=lc, pp_clash_frac=pc, motif_rmsd=mr, **{k: r[k] for k in ("volume_ratio", "ligand_sasa_frac", "n_shell", "loop_frac_shell", "pocket_volume")}))
for f in sorted(glob.glob(f"{H}/rfd4_runs/out/20261121_*/inference_samples/*/*.cif")):
    arm = f.split("/out/20261121_")[1].split("/")[0]
    if arm != "smoke": score(load(open(f).read()), tool="rfd4", arm=arm, design=os.path.basename(f))
for d in sorted(glob.glob(f"{H}/rfd3/out/20261121_*_*")):
    base = os.path.basename(d)[9:]
    if base.startswith("test") or "_" not in base or base.endswith(".csv"): continue
    key, arm = base.split("_", 1); items = []
    if os.path.exists(f"{d}/raw.tar.gz"):
        with tarfile.open(f"{d}/raw.tar.gz") as t: items = [(x.name, t.extractfile(x).read()) for x in t if x.name.endswith(".cif.gz")]
    else: items = [(f, open(f, "rb").read()) for f in glob.glob(f"{d}/**/*.cif.gz", recursive=True)]
    for name, raw in items:
        if f"_{key}_150_260_model" in name: score(load(raw, True), tool="rfd3", arm=f"{key}_{arm}", design=os.path.basename(name))
df = pd.DataFrame(rows)
t = pd.read_csv(f"{H}/rfd4_runs/out/20261121_eval_rfd4_tmol.csv").rename(columns={"arm": "arm_"}); t["tool"] = "rfd4"
df = df.merge(t[["tool", "arm_", "design", "tmol_ddg"]].rename(columns={"arm_": "arm"}), on=["tool", "arm", "design"], how="left")
f = {"len": df.length.between(*THR["len"]), "lig_clash": df.lig_clash_atoms == 0, "pp_clash": df.pp_clash_frac <= THR["pp_clash_frac"],
     "motif": df.motif_rmsd <= THR["motif_rmsd"], "pocket": df.volume_ratio.between(*THR["ratio"]), "buried": df.ligand_sasa_frac <= THR["sasa"],
     "shell": df.n_shell >= THR["n_shell"], "tmol": (df.tmol_ddg < 0) | df.tmol_ddg.isna()}
for k, v in f.items(): df["ok_" + k] = v
df["pass_all"] = np.logical_and.reduce(list(f.values()))
os.makedirs(f"{H}/filtered", exist_ok=True); df.to_csv(f"{H}/filtered/20261121_design_metrics.csv", index=False)
pd.set_option("display.width", 220)
print(df.groupby(["tool", "arm"]).agg(n=("design", "size"), **{("pass_" + k): ("ok_" + k, "sum") for k in f}, pass_all=("pass_all", "sum")).to_string())
print(df[df.pass_all][["tool", "arm", "design", "length", "volume_ratio", "loop_frac_shell", "motif_rmsd", "tmol_ddg"]].round(2).to_string())
