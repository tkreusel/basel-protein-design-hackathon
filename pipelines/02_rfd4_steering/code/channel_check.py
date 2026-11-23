"""Is the ligand pocket connected to bulk solvent? For probe radii r, grid (0.5 A) points whose distance to every protein heavy atom is >= 1.7 + r are 'free'
(ligand removed, since it is replaced by substrate); start = free points within 3.0 A of any ligand atom; connected (26-neighbour) to the padded box boundary => channel of radius r.
Reports the bottleneck radius = largest r with a path (0 = none even at r=0.5; the cavity itself may be smaller than the probe, flagged by 'start_free').
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I channel_check.py -> runs/20261121_r4_SUMMARY/channel_check.csv"""
import glob, os, numpy as np, pandas as pd, biotite.structure.io.pdbx as pdbx, biotite.structure.io.pdb as pdb
from scipy import ndimage
from scipy.spatial import cKDTree
H = os.path.expanduser("~/HACKATHON_ETH"); R = [0.5, 0.7, 0.9, 1.1, 1.4]
def analyse(a):
    a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]; lig, p = a[a.hetero], a[~a.hetero]
    lo, hi = p.coord.min(0) - 8, p.coord.max(0) + 8; ax = [np.arange(lo[i], hi[i], 0.5) for i in range(3)]
    G = np.stack(np.meshgrid(*ax, indexing="ij"), -1); sh = G.shape[:3]; pts = G.reshape(-1, 3)
    # only near the ligand matters for the start; the rest for connectivity
    dp = cKDTree(p.coord).query(pts, distance_upper_bound=6)[0].reshape(sh); dl = cKDTree(lig.coord).query(pts)[0].reshape(sh); out = {}; best = 0.0; start_free = {}
    for r in R:
        free = dp >= 1.7 + r; lab, n = ndimage.label(free, structure=np.ones((3, 3, 3)))
        st = free & (dl <= 3.0); start_free[r] = int(st.sum())
        ids = set(np.unique(lab[st])) - {0}
        edge = set(np.unique(np.concatenate([lab[0].ravel(), lab[-1].ravel(), lab[:, 0].ravel(), lab[:, -1].ravel(), lab[:, :, 0].ravel(), lab[:, :, -1].ravel()]))) - {0}
        ok = bool(ids & edge); out[f"channel_r{r}"] = ok
        if ok: best = r
    out["bottleneck_radius"] = best; out["start_free_r0.5"] = start_free[0.5]; out["start_free_r1.4"] = start_free[1.4]
    return out
rows = []
for f in sorted(glob.glob(H + "/crystals/[0-9]*.pdb")):
    pid = os.path.basename(f)[:4]; st = pdb.PDBFile.read(f)
    for m in (range(1, 4) if pid == "29SB" else [1]):
        a = st.get_structure(model=m, altloc="first"); a = a[a.element != "H"]; prot = a[~a.hetero & ~np.isin(a.res_name, ("HOH", "WAT"))]; ligs = a[a.res_name == "6NT"]
        for ch, rid in sorted({(c, r) for c, r in zip(ligs.chain_id, ligs.res_id)}):
            l = ligs[(ligs.chain_id == ch) & (ligs.res_id == rid)]
            if len(l) < 12: continue
            l.hetero[:] = True
            import biotite.structure as s
            rows.append({"set": "crystal", "id": f"{pid}_m{m}_{ch}{rid}", **analyse(s.concatenate([prot, l]))})
for f in sorted(glob.glob(H + "/runs/20261121_r4_*/structures/*/designs/*.cif")):
    a = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False); c = f.split("/runs/20261121_r4_")[1].split("/")
    rows.append({"set": c[0] + "/" + c[2], "id": os.path.basename(f)[:-4], **analyse(a)})
d = pd.DataFrame(rows); d.to_csv(H + "/runs/20261121_r4_SUMMARY/channel_check.csv", index=False)
pd.set_option("display.width", 220)
print(d.groupby("set")[["bottleneck_radius", "channel_r0.7", "channel_r1.1", "channel_r1.4", "start_free_r0.5", "start_free_r1.4"]].agg(["mean"]).round(2).to_string())
dm = pd.read_csv(H + "/runs/20261121_r4_SUMMARY/design_metrics_all.csv"); pa = dm[dm.pass_all].design.tolist()
print("passing designs:"); print(d[d.id.isin(pa)][["set", "id", "bottleneck_radius", "start_free_r0.5"]].to_string())
