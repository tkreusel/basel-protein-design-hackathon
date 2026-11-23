"""Did partial diffusion beat the parents? Compare on the SAME footing.

Children get one Boltz seed first, so they are compared with their parent on the
three single-seed components only: composite3 = mean(s_ipde, s_geometry,
s_burial). (Parent values are their multi-seed means.) A child counts as an
improvement when composite3 > parent's composite3. The best improvements (<= --per
per parent, <= --top overall) go to the multi-seed pass, where consistency is
measured and the final composite_abs is computed.

Usage: python -I compare_children.py <run dir> [--top 36 --per 4]"""
import argparse, csv, re
from pathlib import Path

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--top", type=int, default=36)
ap.add_argument("--per", type=int, default=4); a = ap.parse_args()
run = Path(a.run)
par = {r["parent"]: r for r in csv.DictReader(open(run / "parents/parents.csv"))}
cand = {r["cand_id"]: r for r in csv.DictReader(open(run / "child_candidates.csv"))}
c3 = lambda r: sum(float(r[k]) for k in ("s_ipde", "s_geometry", "s_burial")) / 3  # noqa: E731
rows = []
for r in csv.DictReader(open(run / "scores_seed0/ranking.csv")):
    if not r.get("composite_abs"):
        continue
    pid = cand[r["cand_id"]]["source"].split(":")[1]
    p = par[pid]
    rows.append(dict(cand_id=r["cand_id"], parent=pid, composite3=round(c3(r), 4),
                     parent_composite3=round(c3(p), 4), gain=round(c3(r) - c3(p), 4),
                     s_ipde=r["s_ipde"], s_geometry=r["s_geometry"], s_burial=r["s_burial"],
                     parent_s_ipde=p["s_ipde"], parent_s_geometry=p["s_geometry"], parent_s_burial=p["s_burial"],
                     flags=r["flags"]))
rows.sort(key=lambda r: -r["gain"])
with open(run / "children_vs_parent.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
n_better = sum(r["gain"] > 0 for r in rows)
print(f"{len(rows)} scored child sequences; {n_better} beat their parent on composite3")
cnt, top = {}, []
for r in rows:
    if r["gain"] <= 0 or r["flags"] or cnt.get(r["parent"], 0) >= a.per or len(top) >= a.top:
        continue
    cnt[r["parent"]] = cnt.get(r["parent"], 0) + 1; top.append(r)
with open(run / "multiseed_candidates.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["cand_id", "sequence", "base_pos", "donor_pos", "wall_pos", "source", "priority"])
    for i, r in enumerate(top):
        c = cand[r["cand_id"]]
        w.writerow([c["cand_id"], c["sequence"], c["base_pos"], c["donor_pos"], c["wall_pos"], c["source"], i])
print(f"{len(top)} child sequences -> multi-seed pass (<= {a.per}/parent)")
for r in top[:10]:
    print(f"  {r['cand_id'][:40]:40} {r['parent']} comp3 {r['composite3']:.3f} vs parent {r['parent_composite3']:.3f} (+{r['gain']:.3f})")
