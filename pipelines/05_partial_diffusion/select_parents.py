"""Pick the parents for partial diffusion: the best-scoring designs we have.

Input: the host-style ranking (scoring/host_score.py output `ranking.csv`) of the
pooled candidates, and the pool's candidates.csv. A "parent" is a BACKBONE: all
sequences refolded from one RFD4 structure share it, so the best-scoring
sequence represents the backbone.

RFD4 partial diffusion needs the annotated `*_generated.cif` that RFD4 itself
wrote (it carries the conditioning columns: fixed motif atoms, burial, H-bond
roles, unindexed groups). Boltz refolds do NOT work as parents. Where each source
keeps that file:

  ours_shortlist   runs/*/manifest.jsonl  design_id -> paths.rfd4_cif
  maxi_steered     handoff runs/20261121_r4_*/structures/**/<design>.cif
  overnight        <overnight run>/raw_gen/inference_samples/<design>/<design>_generated.cif

A parent whose file is missing or has no condition columns is listed in
`parents_excluded.csv` with the reason, never dropped silently.

Writes <out>/parents.csv, <out>/parents/<pid>.cif (copies, so the run folder is
self-contained) and <out>/partial_spec.yaml (one RFD4 example per parent).

Usage: python -I select_parents.py <ranking.csv> <out dir> [--n 12] [--overnight <run dir>]
"""

import argparse
import csv
import json
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HAND = ROOT / "imports/steering_handoff_20261121/20261121_rfd4_steering_handoff"

ap = argparse.ArgumentParser()
ap.add_argument("ranking"); ap.add_argument("out")
ap.add_argument("--n", type=int, default=12)
ap.add_argument("--overnight", help="overnight run dir (for its raw_gen structures)")
a = ap.parse_args()
out = Path(a.out); (out / "parents").mkdir(parents=True, exist_ok=True)

# ---- where each backbone's annotated RFD4 output lives ---------------------
# design ids are hash(cell, index, seed) and COLLIDE across runs (same id in the
# baseline, full-motif and tip-atom runs), so each run has its own map, chosen by
# the shortlist's name prefix ("full-" = job 21, "tip-" = job 23).
ours = {"full": {}, "tip": {}}
for pref, tag in (("full", "BURIED_job21"), ("tip", "TIPATOMS-BURIED_job23")):
    for m in (ROOT / "runs").glob(f"*{tag}*/manifest.jsonl"):
        for line in m.read_text().splitlines():
            j = json.loads(line)
            p = (j.get("paths") or {}).get("rfd4_cif")
            if p:
                ours[pref][j["design_id"]] = ROOT / p
maxi = {}
for p in (HAND / "runs").glob("20261121_r4_*/structures/**/*.cif"):
    maxi.setdefault(p.stem, p)


def resolve(cand_id: str, source: str):
    """-> (backbone key, cif path or None)"""
    if source == "ours_shortlist":
        pref, key = re.sub(r"_s\d+$", "", cand_id).split("-", 1)
        return f"{pref}-{key}", ours[pref].get(key)
    if source == "maxi_steered":
        key = re.sub(r"__(mpnn\d+|rfd4)$", "", cand_id)
        return key, maxi.get(key)
    if source == "overnight" and a.overnight:
        key = re.sub(r"__(mpnn\d+|rfd4)$", "", cand_id)
        return key, Path(a.overnight) / "raw_gen/inference_samples" / key / f"{key}_generated.cif"
    return cand_id, None


def annotated(p: Path) -> bool:
    with open(p) as fh:
        return "condition_coordinate_1_atom" in fh.read(400000)


rows = [r for r in csv.DictReader(open(a.ranking)) if r.get("composite_abs")]
rows.sort(key=lambda r: -float(r["composite_abs"]))
seen, parents, excl = set(), [], []
for r in rows:
    bad = [f for f in r["flags"].split(";") if f and f != "seed-missing"]
    key, cif = resolve(r["cand_id"], r["source"])
    if key in seen:
        continue                       # a better sequence of this backbone is already in
    seen.add(key)
    if bad:
        excl.append((r["cand_id"], f"flags: {';'.join(bad)}")); continue
    if cif is None or not cif.exists():
        excl.append((r["cand_id"], "annotated RFD4 structure not found")); continue
    if not annotated(cif):
        excl.append((r["cand_id"], "no C_CRD (fixed-atom) column in the RFD4 output: motif would not stay pinned")); continue
    parents.append((r, key, cif))
    if len(parents) >= a.n:
        break

with open(out / "parents.csv", "w") as fh:
    w = csv.writer(fh)
    w.writerow(["parent", "backbone", "best_cand", "source", "composite_abs", "s_ipde", "s_geometry",
                "s_burial", "s_consistency", "n_seeds", "burial_CG", "ipde", "cif"])
    for i, (r, key, cif) in enumerate(parents, 1):
        pid = f"p{i:02d}"
        shutil.copy2(cif, out / "parents" / f"{pid}.cif")
        w.writerow([pid, key, r["cand_id"], r["source"]] +
                   [r.get(k, "") for k in ("composite_abs", "s_ipde", "s_geometry", "s_burial",
                                            "s_consistency", "n_seeds", "burial_CG", "ipde")] + [str(cif)])
with open(out / "parents_excluded.csv", "w") as fh:
    csv.writer(fh).writerows([("cand_id", "reason")] + excl)
with open(out / "partial_spec.yaml", "w") as fh:
    for i in range(1, len(parents) + 1):
        fh.write(f"p{i:02d}:\n  input: {out / 'parents' / f'p{i:02d}.cif'}\n  input_as_annotated_cif: true\n")
from collections import Counter
print(f"{len(parents)} parents:", dict(Counter(r["source"] for r, _, _ in parents)),
      f"| excluded {len(excl)}: {dict(Counter(e[1] for e in excl))}")
for i, (r, key, _) in enumerate(parents, 1):
    print(f"  p{i:02d} {r['cand_id'][:44]:44} {r['source']:15} comp {float(r['composite_abs']):.3f}")
