"""Merged pool: Maxi overnight top30 (maxi_top30_scored.csv) + all of our scored candidates (ranking_v2.csv), de-duplicated by sequence,
top 80 by judging stand-in score. Writes merged80_pool.csv (+ ids file for eligibility)."""
import csv
from pathlib import Path
H = Path(__file__).resolve().parent
rows = {}
for r in csv.DictReader(open(H / "maxi_top30_scored.csv")):
    rows[r["sequence"]] = dict(src="maxi_overnight", orig=r["name"], judge=float(r["total_mult"]), sequence=r["sequence"])
for r in csv.DictReader(open(H / "ranking_v2.csv")):
    if r["sequence"] in rows and rows[r["sequence"]]["src"] == "maxi_overnight": continue
    j = float(r["total_mult"])
    if r["sequence"] not in rows or rows[r["sequence"]]["judge"] < j:
        rows[r["sequence"]] = dict(src=r["source"], orig=r["cand_id"], judge=j, sequence=r["sequence"])
top = sorted(rows.values(), key=lambda x: -x["judge"])[:80]
with open(H / "merged80_pool.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(top[0])); w.writeheader(); w.writerows(top)
from collections import Counter
print(len(rows), "unique;", Counter(t["src"] for t in top), "judge cutoff", top[-1]["judge"])
