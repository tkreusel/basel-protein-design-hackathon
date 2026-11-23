"""Turn a raw model run into a human-readable run folder.

This enforces the workspace convention (see runs/README.md). Every run folder
gets the SAME shape, so a person can open any one of them and know where to
look without reading code:

    <run folder>/
        SUCCESS.md       what success means for this stage, in words + numbers
        summary.csv      one row per attempt, with verdict and failure reason
        attempts/        EVERY attempt, filename carries the verdict
        success/         only the successes, ranked best-first
        manifest.jsonl   machine-readable DesignRecord per attempt
        designs/         per-attempt record.json

Success for stage 1 (RFD4 backbone design) is defined as BOTH:
  1. the active site survived  -- catalytic base, oxyanion donor and aromatic
     wall all within the distances calibrated from the HG3 crystal series
  2. the ligand is in the right place -- buried in a pocket, not stuck on the
     protein surface

Criterion 2 exists because 20/32 designs in the first real batch satisfied
every distance in criterion 1 while leaving up to 54% of the ligand
solvent-exposed. Distances alone do not mean there is a pocket.

Usage: python -I finalize_run.py <project root> <rfd4 out_dir> <run folder>
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402
from kemp.record import read_manifest  # noqa: E402

GATE = "stage1_motif"


def success_doc(root: Path, n_total: int, n_success: int) -> str:
    thr = g.load_thresholds(root / "configs" / "thresholds.yaml")
    geo = thr["geometry"]
    bur = geo["ligand_burial"]
    b, d, w = (geo["base_to_abstraction"], geo["donor_to_oxyanion"],
               geo["wall_stacking"])
    return f"""# What counts as SUCCESS in this run

Stage: **RFD4-Proteina backbone + sequence co-design** around the HG3.17 Kemp
eliminase active site.

Result: **{n_success} of {n_total} attempts succeeded.**

A design succeeds only if **both** conditions hold. Both are checked by
`src/extract_designs.py::apply_motif_gate`; nothing is judged by eye.

## 1. The active site survived

The three catalytic groups must still sit on the ligand at catalytic distance.
Thresholds are calibrated from the HG3 crystal series, not guessed
(`configs/thresholds.yaml`, `docs/MILESTONE1.md`).

| measurement | must be | ideal | HG3.17 crystal (7K4Z) |
| --- | --- | --- | --- |
| catalytic base Asp OD -> 6NT N3 | <= {b['pass_max']} A | {b['ideal'][0]}-{b['ideal'][1]} A | 2.53 A |
| oxyanion donor Gln NE2 -> 6NT N1 | <= {d['pass_max']} A | {d['ideal'][0]}-{d['ideal'][1]} A | 3.10 A |
| aromatic wall centroid -> benzo centroid | <= {w['pass_max']} A | {w['ideal'][0]}-{w['ideal'][1]} A | 5.53 A |
| closest protein-ligand heavy atom pair | >= {thr['stage1_backbone']['min_protein_ligand_contact']} A (no clash) | | 2.53 A |

## 2. The ligand is in the right place (buried in a pocket)

Distances alone are not enough. In the first real batch, 20/32 designs passed
every distance above with the ligand sitting on the OUTSIDE of the protein, the
catalytic residues merely reaching onto it. Those are not enzymes.

| measurement | must be | ideal | crystals (7K4Z/4BS0/7K4X) |
| --- | --- | --- | --- |
| ligand relative SASA (0=buried, 1=exposed) | <= {bur['rel_sasa_pass_max']} | {bur['rel_sasa_ideal'][0]}-{bur['rel_sasa_ideal'][1]} | 0.002-0.011 |
| enclosure (fraction of directions blocked by protein) | >= {bur['enclosure_pass_min']} | {bur['enclosure_ideal'][0]}-{bur['enclosure_ideal'][1]} | 0.992-1.000 |
| protein heavy atoms within 5 A | >= {bur['n_contacts_pass_min']} | | 65-68 |

Cut-offs are deliberately looser than the crystals: a de novo pocket need not
be as snug as an evolved one, but it must be a pocket.

## 3. Length

Competition limit: <= {thr['stage6_compliance']['length_max']} residues.

## What success here does NOT mean

- **Not** that the sequence folds to this structure. That needs the Boltz-2
  refold (stage 5); until then the structure is RFD4's own output and judging
  its geometry is close to circular.
- **Not** that the design is eligible to submit. Sequence identity and
  TM-score to published structures are unchecked.
- **Not** that it is catalytically active. No physics scoring (tmol) or
  active-site ensemble (PLACER) has run.

## Folders

- `attempts/` - every attempt. Filename says the verdict, so a directory
  listing alone tells you what happened.
- `success/` - copies of the successes only, ranked best-first by ligand
  burial. `rank01_` is the best.
- `summary.csv` - one row per attempt with all numbers and the failure reason.
"""


def main(root: Path, rfd4_out: Path, run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)

    # The manifest is append-only (so a killed stage can resume), which means
    # re-running finalize on the same folder would DOUBLE it -- 32 attempts
    # became 64 the first time this happened. Finalizing is a fresh layout of
    # one completed run, not a resume, so clear the derived files first.
    for stale in (run_dir / "manifest.jsonl",):
        stale.unlink(missing_ok=True)
    for stale_dir in (run_dir / "designs",):
        if stale_dir.exists():
            shutil.rmtree(stale_dir)

    # Re-run extraction into this run folder so manifest/designs live with it.
    subprocess.run(
        [sys.executable, "-I", str(Path(__file__).parent / "extract_designs.py"),
         str(root), str(rfd4_out), str(run_dir)],
        check=True,
    )

    recs = list(read_manifest(run_dir / "manifest.jsonl"))
    attempts = run_dir / "attempts"
    success = run_dir / "success"
    for d in (attempts, success):
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)

    def short_reason(r) -> str:
        reasons = r.gates.get(f"{GATE}_reasons") or []
        if not reasons:
            return "unknown"
        t = reasons[0]
        if "rel_sasa" in t or "enclosure" in t or "contacts" in t:
            return "ligand-on-surface"
        # Distinguish too-far from too-close: a strained short contact and a
        # missing contact are different problems, and a label that says
        # "too-far" for a 2.41 A contact actively misleads a human reading
        # the folder listing.
        too_close = "hard_fail_min" in t or "clash" in t
        if "base->abstraction" in t:
            return "base-too-close" if too_close else "base-too-far"
        if "donor->oxyanion" in t:
            return "donor-too-close" if too_close else "donor-too-far"
        if "wall" in t:
            return "wall-too-close" if too_close else "wall-too-far"
        if "length" in t:
            return "too-long"
        if "clash" in t or "contact" in t:
            return "clash"
        return "gate-fail"

    ok = [r for r in recs if r.status == "ok" and r.gates.get(GATE) == "pass"]
    bad = [r for r in recs if not (r.status == "ok" and r.gates.get(GATE) == "pass")]
    ok.sort(key=lambda r: r.scores.get("geom.rel_sasa", 9.9))

    # --- attempts/: every attempt, verdict in the filename -------------
    for n, r in enumerate(sorted(recs, key=lambda x: x.design_id), start=1):
        src = Path(r.paths.get("rfd4_cif", ""))
        if not src.exists():
            continue
        if r.status != "ok":
            verdict = "ERROR"
        elif r.gates.get(GATE) == "pass":
            verdict = "SUCCESS"
        else:
            verdict = f"FAIL-{short_reason(r)}"
        shutil.copy2(src, attempts / f"attempt{n:02d}_{r.design_id}_{verdict}.cif")

    # --- success/: ranked copies ----------------------------------------
    fasta = []
    for rank, r in enumerate(ok, start=1):
        src = Path(r.paths.get("rfd4_cif", ""))
        if src.exists():
            shutil.copy2(src, success / f"rank{rank:02d}_{r.design_id}.cif")
        pdb = Path(r.paths.get("backbone_pdb", ""))
        if pdb.exists():
            shutil.copy2(pdb, success / f"rank{rank:02d}_{r.design_id}.pdb")
        cat = r.catalytic_residues
        wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat.get("wall")
        s = r.scores
        fasta.append(
            f">rank{rank:02d}_{r.design_id} len={r.length} "
            f"base={cat.get('base')} donor={cat.get('donor')} wall={wall} "
            f"rel_sasa={s.get('geom.rel_sasa'):.3f} "
            f"enclosure={s.get('geom.enclosure'):.3f}\n{r.sequence}"
        )
    (success / "sequences.fasta").write_text("\n".join(fasta) + ("\n" if fasta else ""))

    # --- summary.csv: every attempt ------------------------------------
    cols = ["attempt", "design_id", "verdict", "failure_reason", "length",
            "base_pos", "donor_pos", "wall_pos",
            "base_dist", "donor_dist", "wall_dist",
            "rel_sasa", "enclosure", "n_contacts", "n_shell", "min_contact",
            "rank_if_success", "cif"]
    rank_of = {r.design_id: i + 1 for i, r in enumerate(ok)}
    with open(run_dir / "summary.csv", "w") as fh:
        fh.write(",".join(cols) + "\n")
        for n, r in enumerate(sorted(recs, key=lambda x: x.design_id), start=1):
            cat = r.catalytic_residues or {}
            wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat.get("wall")
            s = r.scores
            passed = r.status == "ok" and r.gates.get(GATE) == "pass"
            row = {
                "attempt": n, "design_id": r.design_id,
                "verdict": "SUCCESS" if passed else ("ERROR" if r.status != "ok" else "FAIL"),
                "failure_reason": "" if passed else short_reason(r),
                "length": r.length,
                "base_pos": cat.get("base"), "donor_pos": cat.get("donor"),
                "wall_pos": wall,
                "base_dist": s.get("geom.base_to_abstraction"),
                "donor_dist": s.get("geom.donor_to_oxyanion"),
                "wall_dist": s.get("geom.wall_stacking"),
                "rel_sasa": s.get("geom.rel_sasa"),
                "enclosure": s.get("geom.enclosure"),
                "n_contacts": s.get("geom.n_contacts"),
                "n_shell": s.get("geom.n_shell"),
                "min_contact": s.get("geom.min_contact"),
                "rank_if_success": rank_of.get(r.design_id, ""),
                "cif": r.paths.get("rfd4_cif", ""),
            }
            fh.write(",".join("" if row[c] is None else str(row[c]) for c in cols) + "\n")

    (run_dir / "SUCCESS.md").write_text(success_doc(root, len(recs), len(ok)))

    print(f"\n{'=' * 60}")
    print(f"run folder: {run_dir}")
    print(f"  attempts/  {len(recs)} attempts ({len(bad)} failed)")
    print(f"  success/   {len(ok)} successes, ranked")
    print(f"  SUCCESS.md what success means here")
    print(f"  summary.csv one row per attempt")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
