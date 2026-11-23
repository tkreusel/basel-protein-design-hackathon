"""Merge parent-candidate scores from different passes into one ranking.

Multi-seed pool (job 61) and the single-seed overnight designs are put on the
same footing by composite3 = mean(s_ipde, s_geometry, s_burial) (the consistency
component only exists for multi-seed). The output's `composite_abs` column IS
composite3 (select_parents sorts by it); the real value is kept as
`composite_abs_full`.

Usage: python -I merge_rankings.py <out.csv> <ranking.csv> [<ranking.csv> ...]"""
import csv, sys
rows = []
for f in sys.argv[2:]:
    for r in csv.DictReader(open(f)):
        if not r.get("composite_abs") or not all(r.get(k) for k in ("s_ipde", "s_geometry", "s_burial")):
            continue
        r["composite_abs_full"] = r["composite_abs"]
        r["composite_abs"] = sum(float(r[k]) for k in ("s_ipde", "s_geometry", "s_burial")) / 3
        r["ranking_file"] = f
        rows.append(r)
cols = list(dict.fromkeys(k for r in rows for k in r))
with open(sys.argv[1], "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
print(f"{len(rows)} scored candidates merged -> {sys.argv[1]}")
