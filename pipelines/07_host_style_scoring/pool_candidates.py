"""Pool candidate sequences from every source into one candidates.csv.

Columns: cand_id, sequence, base_pos, donor_pos, wall_pos, source, priority

Only candidates the host-style scorer can judge are included: an Asp base and a
Gln donor at known positions (the HG3.17 family). Others are listed in
`excluded.csv` with the reason, never silently dropped.

Sources
  ours_shortlist   shortlist_for_placer/ (31 refold-validated designs)
  maxi_steered     Maxi's steered designs from the handoff tarball (HG3.17
                   motifs only; GDM Glu/Ser designs cannot be scored here)
  overnight        the unsteered control run (sequences.csv + selection.csv)

priority: 0 = ours_shortlist, 1 = maxi_steered, 2 = overnight (sorted by
stage-1 ranking). Lower runs first in the multi-seed Boltz job, so if time runs
out the highest-priority candidates are the ones that finished.

Usage: python -I pool_candidates.py <out dir> [--overnight <run dir>]
"""

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HAND = ROOT / "imports/steering_handoff_20261121/20261121_rfd4_steering_handoff/runs/20261121_r5_downstream"

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("--overnight")
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
rows, excl = [], []


def add(cid, seq, base, donor, wall, source, prio):
    seq = seq.strip()
    if not seq or not set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"):
        excl.append((cid, source, "nonstandard residue")); return
    if len(seq) > 260:
        excl.append((cid, source, f"length {len(seq)} > 260")); return
    if seq[base - 1] != "D" or seq[donor - 1] != "Q":
        excl.append((cid, source, f"base/donor not D/Q ({seq[base-1]}/{seq[donor-1]})")); return
    rows.append(dict(cand_id=cid, sequence=seq, base_pos=base, donor_pos=donor,
                     wall_pos=wall, source=source, priority=prio))


# 1 ours: shortlist
seqs, cur = {}, None
for line in (ROOT / "shortlist_for_placer/sequences.fasta").read_text().splitlines():
    if line.startswith(">"):
        cur = line[1:].split()[0]; seqs[cur] = ""
    else:
        seqs[cur] += line.strip()
for r in csv.DictReader(open(ROOT / "shortlist_for_placer/manifest.csv")):
    add(r["design_id"], seqs[r["design_id"]], int(r["base_pos"]), int(r["donor_pos"]),
        int(r["wall_pos"]), "ours_shortlist", 0)

# 2 Maxi's steered designs (HG3.17 motifs only)
for r in csv.DictReader(open(HAND / "sequences.csv")):
    if not r["design"].startswith("hg"):
        excl.append((r["seq_id"], "maxi_steered", "GDM Glu/Ser theozyme (no Asp/Gln)")); continue
    ids = [int(x[1:]) for x in r["fixed"].split()]
    add(r["seq_id"], r["sequence"], ids[0], ids[1], ids[2] if len(ids) > 2 else -1, "maxi_steered", 1)

# 3 overnight unsteered control
if a.overnight:
    run = Path(a.overnight)
    rank = {r["design"]: int(r["select_rank"]) for r in csv.DictReader(open(run / "selection/selection.csv"))}
    for r in csv.DictReader(open(run / "sequences.csv")):
        ids = [int(x[1:]) for x in r["fixed"].split()]
        add(r["seq_id"], r["sequence"], ids[0], ids[1], ids[2], "overnight", 2 + rank[r["design"]] / 1000)

rows.sort(key=lambda r: (r["priority"], r["cand_id"]))
with open(out / "candidates.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
with open(out / "excluded.csv", "w") as fh:
    csv.writer(fh).writerows([("cand_id", "source", "reason")] + excl)
from collections import Counter
print(f"{len(rows)} candidates:", dict(Counter(r['source'] for r in rows)))
print(f"{len(excl)} excluded:", dict(Counter((e[1], e[2].split(' (')[0]) for e in excl)))
print(f"wrote {out / 'candidates.csv'}")
