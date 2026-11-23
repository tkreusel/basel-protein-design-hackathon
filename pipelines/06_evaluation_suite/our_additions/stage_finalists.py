"""Copy the Boltz-2 TSA (6NT) refold files of the finalists into boltz_root/<id>/ (the layout the Maxi suite
scans with rglob("*_model_*.cif")) and write Boltz inputs for the substrate (H5J) and product (6VP) arms.
Usage: python -I stage_finalists.py   (finalists listed below, provisional, chosen 2026-11-22 07:00 UTC)"""
import csv, glob, shutil, sys
from pathlib import Path
import yaml

R = Path("/mnt/work/projects/tobi_v1/runs")
FINALISTS = {   # id -> why
    "hgnt__steer_base_tmol__no_trp_0_beam3_5__rfd4": "Maxi steered, 5-seed rank 1",
    "hgnt__steer_base_tmol__no_trp_0_beam3_6__rfd4": "Maxi steered, 5-seed rank 2",
    "hgnt__steer_base_tmol__no_trp_0_beam3_5__mpnn3": "Maxi steered, 5-seed rank 3",
    "full-ce941bbd5c43": "ours, 5-seed rank 4",
    "full-af14b3c83c24_s4": "ours, 5-seed rank 5",
    "kemp4_with_trp_128__rfd4": "unsteered overnight, single-seed rank 1",
    "kemp4_with_trp_126__mpnn2": "unsteered overnight, single-seed top",
    "kemp4_with_trp_248__mpnn2": "unsteered overnight, single-seed top",
    "kemp4_with_trp_230__mpnn1": "unsteered overnight, single-seed top",
    "p08_f0.5_2__mpnn3": "partial diffusion child, best gain over parent",
    "p05_f0.3_0__mpnn2": "partial diffusion child",
    "p11_f0.5_4__mpnn1": "partial diffusion child",
}
SEARCH = [R / "2026-11-22_host-style-multiseed_pool-v1_5seeds/seed_0",
          R / "2026-11-21_unsteered-control_job60_250designs_seed7/boltz_out",
          *sorted(R.glob("2026-11-22_partial-diffusion_job62*/boltz_out"))]
SEQ = {}
for f in (R / "2026-11-21_unsteered-control_job60_250designs_seed7/sequences.csv",
          *sorted(R.glob("2026-11-22_partial-diffusion_job62*/sequences.csv"))):
    SEQ.update({r["seq_id"]: r["sequence"] for r in csv.DictReader(open(f))})
SEQ.update({r["cand_id"]: r["sequence"] for r in csv.DictReader(open("/mnt/work/projects/tobi_v1/scoring/pool_v1/candidates.csv"))})
out = Path("boltz_root"); arms = Path("boltz_arms_in"); arms.mkdir(exist_ok=True)
rows = []
for cid, why in FINALISTS.items():
    hit = [p for s in SEARCH for p in glob.glob(f"{s}/**/predictions/{cid}", recursive=True)]
    if not hit:
        print(f"MISSING {cid}"); continue
    d = out / cid; d.mkdir(parents=True, exist_ok=True)
    for f in Path(hit[0]).iterdir():
        shutil.copy2(f, d / f.name)
    for arm, ccd in (("H5J", "H5J"), ("6VP", "6VP")):
        (arms / f"{cid}__{arm}.yaml").write_text(yaml.safe_dump({"version": 1, "sequences": [
            {"protein": {"id": "A", "sequence": SEQ[cid], "msa": "empty"}}, {"ligand": {"id": "B", "ccd": ccd}}]}, sort_keys=False))
    rows.append((cid, why, hit[0]))
csv.writer(open("finalists.csv", "w")).writerows([("id", "why", "source")] + rows)
print(len(rows), "finalists staged ->", out, "; arm inputs ->", arms)
