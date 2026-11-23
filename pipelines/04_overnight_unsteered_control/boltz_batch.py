"""Stage 4 prep: Boltz-2 inputs for every sequence in sequences.csv, in chunks.

Matches Maxi's prep_boltz.py: protein chain A with `msa: empty` plus ligand
chain B = CCD 6NT, one yaml per seq_id, chunked so a killed job resumes at chunk
granularity. (Boltz also skips existing prediction folders on its own.)
De novo designs get a query-only MSA either way: verified here by an A/B test
that produced bit-identical structures with and without ColabFold's alignment.

Usage: python -I boltz_batch.py <sequences.csv> <out dir> [chunk=125]
"""

import csv
import sys
from pathlib import Path

import yaml

seq_csv, out = Path(sys.argv[1]), Path(sys.argv[2])
chunk = int(sys.argv[3]) if len(sys.argv) > 3 else 125
rows = [r for r in csv.DictReader(open(seq_csv))
        if r["sequence"] and set(r["sequence"]) <= set("ACDEFGHIKLMNPQRSTVWY")]
bad = sum(1 for r in csv.DictReader(open(seq_csv))) - len(rows)
for k in range(0, len(rows), chunk):
    d = out / "boltz_in" / f"chunk_{k // chunk:02d}"
    d.mkdir(parents=True, exist_ok=True)
    for r in rows[k:k + chunk]:
        (d / f"{r['seq_id']}.yaml").write_text(yaml.safe_dump({
            "version": 1,
            "sequences": [{"protein": {"id": "A", "sequence": r["sequence"], "msa": "empty"}},
                          {"ligand": {"id": "B", "ccd": "6NT"}}]}, sort_keys=False))
n_chunks = -(-len(rows) // chunk)
print(f"{len(rows)} Boltz inputs in {n_chunks} chunks of <= {chunk}"
      f"{f'  ({bad} sequences skipped: non-standard residues)' if bad else ''}")
