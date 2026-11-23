"""Write Boltz-2 inputs for the SUCCESS designs of a finalized RFD4 run.

Only designs that passed the stage-1 gate are refolded. Refolding rejects
wastes GPU and pollutes the comparison -- that mistake has already been made
once (runs/_superseded/..._job14_CANCELLED-ran-on-gate-failed-designs).

Carries each design's catalytic residue positions forward into a sidecar JSON,
because the refold has to be scored at those positions and Boltz output
carries no provenance of its own.

  msa: empty   single-sequence mode (boltz schema.py:1127). De novo designs
               have no database homologues, so an MMseqs2 search would return
               the query row and nothing else, and single-sequence mode avoids
               needing the 28 GB-RAM node.
  ccd: 6NT     Boltz-2 takes the ligand by CCD code, verified in the
               organisers' working example. No SDF needed on this path.

Usage: python -I prepare_refold.py <run dir> <out dir>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp.record import read_manifest  # noqa: E402

GATE = "stage1_motif"
LIGAND_CCD = "6NT"


def main(run_dir: Path, out_dir: Path) -> None:
    recs = list(read_manifest(run_dir / "manifest.jsonl"))
    ok = [r for r in recs
          if r.status == "ok" and r.gates.get(GATE) == "pass" and r.sequence]
    ok.sort(key=lambda r: r.scores.get("geom.rel_sasa", 9.9))

    out_dir.mkdir(parents=True, exist_ok=True)
    sidecar = {}

    for rank, r in enumerate(ok, start=1):
        cat = r.catalytic_residues
        wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat.get("wall")
        name = f"{r.design_id}"
        spec = {
            "version": 1,
            "sequences": [
                {"protein": {"id": "A", "sequence": r.sequence, "msa": "empty"}},
                {"ligand": {"id": "B", "ccd": LIGAND_CCD}},
            ],
        }
        (out_dir / f"{name}.yaml").write_text(yaml.safe_dump(spec, sort_keys=False))
        sidecar[name] = {
            "design_id": r.design_id,
            "stage1_rank": rank,
            "length": r.length,
            # 1-indexed sequence positions. Boltz renumbers the protein 1..N,
            # so these carry over directly to the refolded structure.
            "catalytic": {"base": cat.get("base"), "donor": cat.get("donor"),
                          "wall": wall},
            "stage1_scores": {k: v for k, v in r.scores.items()
                              if k.startswith("geom.")},
            "source_run": run_dir.name,
        }

    # The sidecar must live OUTSIDE out_dir. `boltz predict <dir>` tries to
    # parse every file in the directory it is given and hard-fails on anything
    # that is not .fasta/.yaml:
    #   RuntimeError: Unable to parse filetype .json, please provide a
    #   .fasta or .yaml file.
    (out_dir.parent / "designs.json").write_text(json.dumps(sidecar, indent=2))
    print(f"source run : {run_dir.name}")
    print(f"total       : {len(recs)} designs")
    print(f"stage-1 pass: {len(ok)}  -> refolding these only")
    print(f"wrote       : {len(ok)} Boltz inputs to {out_dir}, designs.json to {out_dir.parent}")
    if not ok:
        raise SystemExit("no passing designs to refold")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
