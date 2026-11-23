"""Stage extra candidates (pool designs never run through the Maxi suite) for substrate/product refolds: copies their 6NT refold into boltz_root/<id>/
and writes Boltz inputs into boltz_arms_in2/ (kept separate from the first 12 finalists). Usage: python -I stage_extra.py"""
import csv, glob, shutil
from pathlib import Path
import yaml
T = Path("/mnt/work/projects/tobi_v1"); R = T / "runs"
have = {p.name for p in Path("boltz_root").iterdir() if p.is_dir()}
pool = list(csv.DictReader(open(T / "submission/merged80_pool.csv")))
SEARCH = [R / "2026-11-22_host-style-multiseed_pool-v1_5seeds/seed_0", R / "2026-11-22_host-style-multiseed_overnight-top24_5seeds/seed_0",
          R / "2026-11-21_unsteered-control_job60_250designs_seed7/boltz_out", *sorted(R.glob("2026-11-22_partial-diffusion_job62*/boltz_out"))]
ids = [(r["orig"], r["sequence"]) for r in pool if r["src"] != "maxi_overnight" and r["judge"] and float(r["judge"]) >= 0.0003]
SEQ = dict(ids)
for extra in ("kemp4_with_trp_226__rfd4", "kemp4_with_trp_126__mpnn1"):   # gate survivors of Maxi's pipeline
    for s in SEARCH:
        for p in glob.glob(f"{s}/**/predictions/{extra}", recursive=True):
            pass
arms = Path("boltz_arms_in2"); arms.mkdir(exist_ok=True)
n = 0
for cid, seq in SEQ.items():
    if cid in have: continue
    hit = [p for s in SEARCH for p in glob.glob(f"{s}/**/predictions/{cid}", recursive=True)]
    if not hit: print("MISSING", cid); continue
    d = Path("boltz_root") / cid; d.mkdir(exist_ok=True)
    for f in Path(hit[0]).iterdir(): shutil.copy2(f, d / f.name)
    for arm in ("H5J", "6VP"):
        (arms / f"{cid}__{arm}.yaml").write_text(yaml.safe_dump({"version": 1, "sequences": [
            {"protein": {"id": "A", "sequence": seq, "msa": "empty"}}, {"ligand": {"id": "B", "ccd": arm}}]}, sort_keys=False))
    n += 1
print(n, "designs staged")
