"""Write Boltz-2 inputs for candidates.csv, in priority-ordered chunks.
Usage: python -I boltz_seed_inputs.py <candidates.csv> <out dir> [chunk=25]"""
import csv, sys, yaml
from pathlib import Path
cands, out = Path(sys.argv[1]), Path(sys.argv[2]); chunk = int(sys.argv[3]) if len(sys.argv) > 3 else 25
rows = sorted(csv.DictReader(open(cands)), key=lambda r: (float(r["priority"]), r["cand_id"]))
for k in range(0, len(rows), chunk):
    d = out / "boltz_in" / f"chunk_{k // chunk:03d}"; d.mkdir(parents=True, exist_ok=True)
    for r in rows[k:k + chunk]:
        (d / f"{r['cand_id']}.yaml").write_text(yaml.safe_dump({"version": 1, "sequences": [
            {"protein": {"id": "A", "sequence": r["sequence"], "msa": "empty"}},
            {"ligand": {"id": "B", "ccd": "6NT"}}]}, sort_keys=False))
print(f"{len(rows)} candidates in {-(-len(rows) // chunk)} chunks of <= {chunk}")
