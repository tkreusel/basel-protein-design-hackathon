"""Build eligibility/candidates.csv (cand_id, sequence, structure) from the finalist list + extra designs + controls.
Controls must FAIL (published HG3.17 sequence/structure): they validate the gate itself."""
import csv, glob, sys
from pathlib import Path
R = Path("/mnt/work/projects/tobi_v1/runs"); T = Path("/mnt/work/projects/tobi_v1")
IDS = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
SEARCH = [R / "2026-11-22_host-style-multiseed_pool-v1_5seeds/seed_0", R / "2026-11-22_host-style-multiseed_overnight-top24_5seeds/seed_0",
          R / "2026-11-21_unsteered-control_job60_250designs_seed7/boltz_out", *sorted(R.glob("2026-11-22_partial-diffusion_job62*/boltz_out"))]
SEQ = {}
for f in [R / "2026-11-21_unsteered-control_job60_250designs_seed7/sequences.csv", *sorted(R.glob("2026-11-22_partial-diffusion_job62*/sequences.csv"))]:
    SEQ.update({r["seq_id"]: r["sequence"] for r in csv.DictReader(open(f))})
SEQ.update({r["cand_id"]: r["sequence"] for r in csv.DictReader(open(T / "scoring/pool_v1/candidates.csv"))})
rows = []
for cid in IDS:
    hit = [p for s in SEARCH for p in glob.glob(f"{s}/**/predictions/{cid}/{cid}_model_0.cif", recursive=True)]
    if not hit or cid not in SEQ: print("MISSING", cid); continue
    rows.append((cid, SEQ[cid], hit[0]))
ds = {r[1]: r[3].replace(" ", "") for r in csv.reader(open("/mnt/work/manually_created_dataset.csv"), delimiter="\t")}
rows.append(("CONTROL_HG3.17_published", ds["HG3.17"], str(T / "benchmark/pdb/7K4Z.cif")))
with open(Path(__file__).parent / "candidates.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["cand_id", "sequence", "structure"]); w.writerows(rows)
print(len(rows), "rows")
