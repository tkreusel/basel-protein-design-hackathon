"""Re-rank candidates under the rules in SUBMISSION.md (2026-11-22 09:46): geometry dominates and the terms MULTIPLY.

SUBMISSION.md: (1) catalytic contact distances carry the most weight, (2) approach geometry (angle and in-plane term,
the WORSE of the two decides), (3) iPDE, (4) base enclosure, (5) pose consistency across repeated predictions; "geometry
outweighs everything else combined; strength in one term does not compensate for failure in another: the terms multiply";
each distance has an optimum, too short is penalised harder than too long.

Our stand-in (the hosts' exact formulas are not public):
  contact   = s(base O..N3) * s(donor..N1)   gaussian around the crystal optimum (configs/host_reference.yaml: 2.55 / 3.05 A),
              sigma_long = 0.25 A (floor), sigma_short = 0.6 * sigma_long  (harsher when too short)
  approach  = min(s_ang_NH, s_oop) for the base  *  min(s_ang_ax, s_oop) for the donor   (worse of angle / in-plane decides)
  geometry  = contact * approach                   (per Boltz seed, then mean over seeds)
  total     = geometry^2 * s_ipde * s_burial * s_consistency        (geometry weighted as "more than everything else combined")
Usage: python -I rank_v2.py <out csv>   (reads all multi-seed runs listed in RUNS)"""
import csv, glob, sys, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scoring")); sys.path.insert(0, str(ROOT / "src"))
import host_score as H

R = ROOT / "runs"
RUNS = [(R / "2026-11-22_host-style-multiseed_pool-v1_5seeds", ROOT / "scoring/pool_v1/candidates.csv"),
        (R / "2026-11-22_host-style-multiseed_overnight-top24_5seeds", ROOT / "scoring/pool_overnight_top24/candidates.csv")]
pd_run = sorted(R.glob("2026-11-22_partial-diffusion_job62*"))
if pd_run:
    RUNS.append((pd_run[0] / "multiseed", pd_run[0] / "multiseed_candidates.csv"))
REF = H.REF
def sig(k):    return REF[k]["sigma"]
def gscore(x, k, short_factor=0.6):
    z = (x - REF[k]["mean"]) / sig(k)
    if k.endswith("dist") and z < 0: z = z / short_factor            # too short: harsher
    return float(np.exp(-0.5 * z * z))
rows = []
for run, cands in RUNS:
    H.CANDS.clear(); H.SEED_DIRS.clear()
    for r in csv.DictReader(open(cands)): H.CANDS[r["cand_id"]] = r
    sd = sorted(p for p in Path(run).glob("seed_*") if p.is_dir()); H.SEED_DIRS.extend(str(p) for p in sd)
    rk = {r["cand_id"]: r for r in csv.DictReader(open(Path(run) / "scores/ranking.csv"))}
    for cid, c in H.CANDS.items():
        g = []
        for s in H.SEED_DIRS:
            try: o = H.one_seed(cid, s)
            except Exception: o = None
            if not o: continue
            contact = gscore(o["b_dist"], "b_dist") * gscore(o["d_dist"], "d_dist")
            base_ap = min(gscore(abs(o["b_ang_NH"]), "b_ang_NH"), gscore(abs(o["b_oop"]), "b_oop"))
            don_ap = min(gscore(abs(o["d_ang_ax"]), "d_ang_ax"), gscore(abs(o["d_oop"]), "d_oop"))
            g.append((contact, base_ap * don_ap, contact * base_ap * don_ap))
        r = rk.get(cid)
        if not g or not r or not r.get("s_ipde"): continue
        contact, approach, geom = (float(np.mean([x[i] for x in g])) for i in range(3))
        cons = float(r["s_consistency"]) if r.get("s_consistency") else 0.0
        tot = geom ** 2 * float(r["s_ipde"]) * float(r["s_burial"]) * cons
        rows.append(dict(cand_id=cid, source=c.get("source", ""), n_seeds=len(g), contact=round(contact, 3), approach=round(approach, 3),
                         geometry=round(geom, 3), s_ipde=round(float(r["s_ipde"]), 3), s_burial=round(float(r["s_burial"]), 3),
                         s_consistency=round(cons, 3), total_mult=round(tot, 4), composite_add=round(float(r["composite_abs"]), 3),
                         b_dist=round(float(r["b_dist"]), 2), d_dist=round(float(r["d_dist"]), 2), sequence=c["sequence"],
                         base_pos=c["base_pos"], donor_pos=c["donor_pos"]))
rows.sort(key=lambda r: -r["total_mult"])
with open(sys.argv[1], "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print(f"{len(rows)} candidates ranked -> {sys.argv[1]}")
for i, r in enumerate(rows[:25], 1):
    print(f"{i:>2} {r['cand_id'][:44]:44} {r['source'][:6]:6} tot={r['total_mult']:.3f} geom={r['geometry']:.2f} (contact {r['contact']:.2f}, appr {r['approach']:.2f}) ipde={r['s_ipde']:.2f} bur={r['s_burial']:.2f} cons={r['s_consistency']:.2f} b={r['b_dist']} d={r['d_dist']} add={r['composite_add']}")
