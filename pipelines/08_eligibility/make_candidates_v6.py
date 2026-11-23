"""candidates_v6.csv = merged top-80 pool designs that are not yet in candidates_v5 (matched by sequence); absolute structure paths."""
import csv, glob, sys
from pathlib import Path
R = Path("/mnt/work/projects/tobi_v1/runs"); T = Path("/mnt/work/projects/tobi_v1")
M30 = "/mnt/work/projects/Maxi/kemp_results_20261122/final_strict_all_gates_ipsae0.8/top30"
have = {r["sequence"] for r in csv.DictReader(open(T / "eligibility/candidates_v5.csv"))}
pool = [r for r in csv.DictReader(open(T / "submission/merged80_pool.csv")) if r["sequence"] not in have]
SEARCH = [R / "2026-11-22_host-style-multiseed_pool-v1_5seeds/seed_0", R / "2026-11-22_host-style-multiseed_overnight-top24_5seeds/seed_0",
          R / "2026-11-21_unsteered-control_job60_250designs_seed7/boltz_out", *sorted(R.glob("2026-11-22_partial-diffusion_job62*/boltz_out"))]
out = []
for r in pool:
    cid = r["orig"]
    if r["src"] == "maxi_overnight":
        n = cid[:2]; hit = glob.glob(f"{M30}/{n}_*/boltz_TSA_6NT_best.cif"); cid = cid.split("_", 1)[1]
        full = cid
    else:
        hit = [p for s in SEARCH for p in glob.glob(f"{s}/**/predictions/{cid}/{cid}_model_0.cif", recursive=True)]
    if not hit: print("MISSING", cid); continue
    out.append((cid, r["sequence"], hit[0]))
with open(T / "eligibility/candidates_v6.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["cand_id", "sequence", "structure"]); w.writerows(out)
print(len(out), "of", len(pool))
