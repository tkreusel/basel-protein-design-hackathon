"""Score Boltz-2 refolds and lay out the run folder (attempts/ + success/).

**Refold success** = the active site is STILL THERE after an independent
prediction from sequence alone:

  1. catalytic base, oxyanion donor and aromatic wall all within the same
     crystal-calibrated distances used at stage 1, measured on the REFOLDED
     structure
  2. the ligand is still buried in a pocket

This is the first non-circular test in the pipeline. RFD4 built the stage-1
geometry, so scoring its own output on that geometry proves little. Boltz-2
sees only the sequence and the ligand.

Boltz confidence (pLDDT, ipTM, ligand ipTM) is recorded but NOT gated on:
plan 6.6 notes absolute confidence is not calibrated on de novo folds.

Usage: python -I finalize_refold.py <project root> <boltz out> <designs.json> <run dir>
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402

PROT_CHAIN = "A"


def find_predictions(boltz_out: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for pred in boltz_out.rglob("predictions/*"):
        if not pred.is_dir():
            continue
        cifs = sorted(pred.glob(f"{pred.name}_model_*.cif"))
        confs = sorted(pred.glob(f"confidence_{pred.name}_model_*.json"))
        if cifs:
            out[pred.name] = {"cif": cifs[0],
                              "conf": confs[0] if confs else None}
    return out


def score(cif: Path, amap: g.AtomMap, cat: dict) -> dict:
    struct = get_structure(CIFFile.read(cif), model=1)
    copies = g.iter_ligand_copies(struct, amap)
    if not copies:
        raise g.MetricUnavailable(f"{cif.name}: ligand absent from prediction")
    lc, lr = copies[0]
    kw = dict(lig_chain=lc, lig_resid=lr, chain=PROT_CHAIN)
    s = {
        "base_to_abstraction": g.base_to_abstraction(
            struct, amap, base_resid=cat["base"], **kw),
        "donor_to_oxyanion": g.donor_to_oxyanion(
            struct, amap, donor_resid=cat["donor"], **kw),
        "wall_stacking": g.wall_stacking(
            struct, amap, wall_resid=cat["wall"], **kw),
        "min_contact": g.min_protein_ligand_contact(
            struct, amap, lig_chain=lc, lig_resid=lr),
    }
    s.update(g.ligand_burial(struct, amap, lig_chain=lc, lig_resid=lr))
    return s


def gate(s: dict, thr: dict) -> tuple[bool, list[str]]:
    geo = thr["geometry"]
    bur = geo["ligand_burial"]
    reasons: list[str] = []
    for key, label in (("base_to_abstraction", "base->abstraction"),
                       ("donor_to_oxyanion", "donor->oxyanion"),
                       ("wall_stacking", "wall stacking")):
        v, spec = s.get(key), geo[key]
        if v is None:
            reasons.append(f"{label}: not measured")
            continue
        if v > spec["pass_max"]:
            reasons.append(f"{label} {v:.2f} > pass_max {spec['pass_max']}")
        lo = spec.get("hard_fail_min")
        if lo is not None and v < lo:
            reasons.append(f"{label} {v:.2f} < hard_fail_min {lo} (too close)")
    if s.get("rel_sasa") is not None and s["rel_sasa"] > bur["rel_sasa_pass_max"]:
        reasons.append(f"ligand rel_sasa {s['rel_sasa']:.3f} > "
                       f"{bur['rel_sasa_pass_max']} (lost the pocket)")
    if s.get("enclosure") is not None and s["enclosure"] < bur["enclosure_pass_min"]:
        reasons.append(f"enclosure {s['enclosure']:.3f} < {bur['enclosure_pass_min']}")
    if s.get("n_contacts") is not None and s["n_contacts"] < bur["n_contacts_pass_min"]:
        reasons.append(f"only {s['n_contacts']:.0f} contacts < "
                       f"{bur['n_contacts_pass_min']}")
    return (not reasons), reasons


def label(reasons: list[str]) -> str:
    if not reasons:
        return "SUCCESS"
    t = reasons[0]
    if "rel_sasa" in t or "enclosure" in t or "contacts" in t:
        return "FAIL-lost-pocket"
    close = "hard_fail_min" in t or "too close" in t
    if "base->" in t:
        return "FAIL-base-too-close" if close else "FAIL-base-drifted"
    if "donor->" in t:
        return "FAIL-donor-too-close" if close else "FAIL-donor-drifted"
    if "wall" in t:
        return "FAIL-wall-drifted"
    return "FAIL-gate"


def main(root: Path, boltz_out: Path, designs_json: Path, run_dir: Path) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    thr = g.load_thresholds(root / "configs" / "thresholds.yaml")
    designs = json.loads(designs_json.read_text())
    preds = find_predictions(boltz_out)

    run_dir.mkdir(parents=True, exist_ok=True)
    attempts, success = run_dir / "attempts", run_dir / "success"
    for d in (attempts, success):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    rows = []
    for name, meta in sorted(designs.items(),
                             key=lambda kv: kv[1].get("stage1_rank", 999)):
        p = preds.get(name)
        row = dict(design_id=meta["design_id"], stage1_rank=meta["stage1_rank"],
                   length=meta["length"], verdict="MISSING", reason="",
                   plddt=None, iptm=None, ligand_iptm=None)
        row.update({f"s1_{k.split('.')[-1]}": v
                    for k, v in meta["stage1_scores"].items()})
        if not p:
            rows.append(row)
            print(f"  [MISSING] {meta['design_id']} no prediction")
            continue
        if p["conf"] and p["conf"].exists():
            c = json.loads(p["conf"].read_text())
            row["plddt"] = c.get("complex_plddt")
            row["iptm"] = c.get("iptm")
            row["ligand_iptm"] = c.get("ligand_iptm")
        try:
            s = score(p["cif"], amap, meta["catalytic"])
            ok, reasons = gate(s, thr)
            row.update({f"rf_{k}": v for k, v in s.items()})
            row["verdict"] = label(reasons)
            row["reason"] = reasons[0] if reasons else ""
        except Exception as exc:  # noqa: BLE001 - record, never drop
            row["verdict"] = "ERROR"
            row["reason"] = f"{type(exc).__name__}: {exc}"
            ok = False
        rows.append(row)

        shutil.copy2(p["cif"],
                     attempts / f"rank{meta['stage1_rank']:02d}_"
                                f"{meta['design_id']}_{row['verdict']}.cif")
        d_ = row.get("rf_base_to_abstraction")
        print(f"  [{row['verdict']:>20}] {meta['design_id']} "
              f"pLDDT={row['plddt'] if row['plddt'] is None else round(row['plddt'],3)} "
              + (f"base={d_:.2f} donor={row['rf_donor_to_oxyanion']:.2f} "
                 f"wall={row['rf_wall_stacking']:.2f} "
                 f"sasa={row['rf_rel_sasa']:.3f}" if isinstance(d_, float) else ""))

    # success/ ranked by refold burial then donor distance
    good = [r for r in rows if r["verdict"] == "SUCCESS"]
    good.sort(key=lambda r: (r.get("rf_rel_sasa", 9.9),
                             abs(r.get("rf_donor_to_oxyanion", 9.9) - 3.05)))
    for rank, r in enumerate(good, start=1):
        src = attempts / f"rank{r['stage1_rank']:02d}_{r['design_id']}_SUCCESS.cif"
        if src.exists():
            shutil.copy2(src, success / f"refold{rank:02d}_{r['design_id']}.cif")

    cols = sorted({k for r in rows for k in r})
    cols = ["design_id", "stage1_rank", "verdict", "reason", "length"] + \
           [c for c in cols if c not in
            ("design_id", "stage1_rank", "verdict", "reason", "length")]
    with open(run_dir / "summary.csv", "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join("" if r.get(c) is None else str(r.get(c, ""))
                              for c in cols) + "\n")

    n = len(rows)
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    (run_dir / "SUCCESS.md").write_text(f"""# What counts as SUCCESS in this run

Stage: **Boltz-2 refold survival** of designs that already passed stage 1.

Result: **{len(good)} of {n} refolded designs kept their active site.**

Verdicts: {", ".join(f"{k} {v}" for k, v in sorted(counts.items()))}

## The question this answers

Does the **sequence alone** encode the active site? RFD4 built the stage-1
geometry, so judging RFD4's own output on it is close to circular. Boltz-2
is given only the amino-acid sequence and the 6NT ligand — no structure, no
MSA (single-sequence mode), no knowledge of the motif. If the catalytic
residues still land on the ligand in its prediction, the sequence genuinely
specifies the site.

## Criteria — identical to stage 1, measured on the REFOLDED structure

1. **Active site survived**: base -> 6NT N3, oxyanion donor -> 6NT N1, and
   aromatic wall -> benzo centroid all within the crystal-calibrated
   distances in `configs/thresholds.yaml`.
2. **Ligand still buried**: relative SASA, enclosure and contact count within
   the same bounds.

Using the stage-1 thresholds unchanged is deliberate: "does it hold up" means
"does it still look like a Kemp eliminase active site", and the crystals
define that. This is a strictly harder test than stage 1, so a lower pass rate
is expected and is not a bug.

## What is recorded but NOT gated

Boltz confidence (`plddt`, `iptm`, `ligand_iptm`) is in `summary.csv` but does
not decide pass/fail. Plan section 6.6: absolute confidence is not calibrated
on de novo folds, so it should inform ranking by *gaps*, not absolute cut-offs.

## Still not established

Catalytic activity. No physics scoring (tmol), no active-site ensemble
(PLACER), and no eligibility check (sequence identity, TM-score).

## Folders

- `attempts/` — every refolded design, verdict in the filename, named by its
  stage-1 rank so you can see whether stage-1 ranking predicted survival.
- `success/` — survivors only, re-ranked by refold burial then donor distance.
- `summary.csv` — stage-1 scores (`s1_*`) beside refold scores (`rf_*`), so
  you can see exactly what moved.
""")

    print(f"\n{'=' * 60}")
    print(f"run folder: {run_dir}")
    print(f"  refolded: {n}   survived: {len(good)}")
    for k, v in sorted(counts.items()):
        print(f"    {k}: {v}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]))
