"""Filter partial-diffusion children BEFORE they cost GPU time (rejection sampling).

Why this is the steering: the active-site tips are pinned (C_CRD), so what
partial diffusion can change is the scaffold around them. The hosts' Burial score
(CA atoms within 12 A of the base) and the pocket criteria are backbone
properties, so they can be checked on the RFD4 output directly. A child survives
only if, relative to ITS PARENT, it

  hard   ligand clash-free; both motif tips within 0.3 A of the theozyme
  keeps  at least as many of Maxi's applied core criteria as the parent
  keeps  base burial (CA within 12 A of CG) >= parent - 2

Survivors are ranked (all criteria passed, then criteria passed, then burial,
then fewer polar/cationic neighbours at the base) and the top M per parent are
written in the format overnight/mpnn_batch.py reads. Everything else is listed
in children_all.csv with the reason, never dropped silently.

Child id = <parent>_f<level>_<k>   (level = partial_trajectory_fraction_steps)

Usage: python -I select_children.py <parents.csv> <out dir> <raw dir>=<level> [...] [--m 13]
"""

import argparse
import csv
import glob
import re
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "overnight")); sys.path.insert(0, str(ROOT / "src"))
from host_geometry import burial_count  # noqa: E402
from select_designs import CRITERIA, _one  # noqa: E402


def measure_one(args):
    path, child, parent, level = args
    r = _one(path)
    r.update(child=child, parent=parent, level=level)
    if r["status"] == "ok":
        s = get_structure(CIFFile.read(path), model=1)
        chain = s.chain_id[~s.hetero][0]
        r.update(burial_CG=burial_count(s, chain=chain, base_resid=int(r["base_resid"]))["burial_CG"])
    return r


ap = argparse.ArgumentParser()
ap.add_argument("parents"); ap.add_argument("out"); ap.add_argument("raw", nargs="+")
ap.add_argument("--m", type=int, default=13); ap.add_argument("--procs", type=int, default=4)
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

jobs = []
for spec in a.raw:
    d, level = spec.rsplit("=", 1)
    for f in sorted(glob.glob(f"{d}/inference_samples/*/*_generated.cif")):
        name = Path(f).parent.name                      # p01_3
        if not re.fullmatch(r"p\d+_\d+", name):
            raise SystemExit(f"unexpected RFD4 output name {name!r} (expected <parent>_<replicate>)")
        pid, k = name.rsplit("_", 1)
        jobs.append((f, f"{pid}_f{level}_{k}", pid, level))
print(f"{len(jobs)} children from {len(a.raw)} raw dir(s)")

parent_rows = {r["parent"]: r for r in csv.DictReader(open(a.parents))}
pm = {}
with Pool(a.procs) as p:
    for pid, pr in parent_rows.items():
        pm[pid] = measure_one((str(Path(a.parents).parent / "parents" / f"{pid}.cif"), pid, pid, "parent"))
    rows = p.map(measure_one, jobs)

kept, allrows = [], []
for pid in parent_rows:
    pr = pm[pid]
    print(f"{pid}: parent n_pass={pr.get('n_pass')} burial_CG={pr.get('burial_CG')} status={pr['status']}")
for r in rows:
    pr = pm[r["parent"]]
    why = []
    if r["status"] != "ok":
        why.append(r["status"][:60])
    elif pr["status"] != "ok":
        why.append("parent not measurable")
    else:
        if not r["ok_clash"]: why.append("ligand-clash")
        if not r["ok_motif"]: why.append("motif-tips-moved")
        if r["n_pass"] < pr["n_pass"]: why.append("fewer-criteria-than-parent")
        if r["burial_CG"] < pr["burial_CG"] - 2: why.append("less-buried-than-parent")
    r["reject"] = "; ".join(why)
    allrows.append(r)
    if not why:
        kept.append(r)

by_parent = {}
for r in kept:
    by_parent.setdefault(r["parent"], []).append(r)
sel = []
for pid, lst in by_parent.items():
    lst.sort(key=lambda r: (-r["n_pass"], -r["burial_CG"], r["base_sc_polar_4p5"] + r["base_cation_6"]))
    for i, r in enumerate(lst[: a.m], 1):
        r["child_rank"] = i; sel.append(r)
sel.sort(key=lambda r: (r["child_rank"], r["parent"]))       # interleave parents: rank-1 of each first
for i, r in enumerate(sel, 1):
    r["select_rank"] = i; r["design"] = r["child"]

cols = ["select_rank", "design", "cif", "parent", "level", "child_rank", "length", "n_pass",
        "burial_CG", "base_resid", "partner_resid", "wall_resid"]
with open(out / "selection.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(sel)
allcols = sorted({k for r in allrows for k in r})
with open(out / "children_all.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=allcols); w.writeheader(); w.writerows(allrows)
with open(out / "parents_measured.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=sorted({k for r in pm.values() for k in r})); w.writeheader(); w.writerows(pm.values())

print(f"\nchildren measured {len(allrows)}, survive filter {len(kept)}, selected {len(sel)} "
      f"(<= {a.m} per parent) -> {out / 'selection.csv'}")
from collections import Counter
print("rejection reasons (a child can have several):",
      dict(Counter(x.strip() for r in allrows if r["reject"] for x in r["reject"].split(";"))))
for pid in parent_rows:
    n_all = sum(1 for r in allrows if r["parent"] == pid)
    print(f"  {pid}: {len(by_parent.get(pid, []))}/{n_all} survive; "
          f"{sum(1 for r in sel if r['parent'] == pid)} selected")
