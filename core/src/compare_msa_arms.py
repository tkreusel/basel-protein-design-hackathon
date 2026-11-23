"""Compare the with-MSA and no-MSA Boltz-2 refolds of one design.

Reports confidence, active-site geometry and the structural difference
between the two predictions. If the ColabFold MSA is query-only (depth 1),
the two arms should be indistinguishable -- this measures whether they are.

Usage: python -I compare_msa_arms.py <project root> <work dir> <design id>
"""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402
from kemp.record import read_manifest  # noqa: E402

import os
# Source run is overridable so this works for any stage-1 run, not just the
# one it was first written against.
RUN = os.environ.get(
    "KEMP_SRC_RUN",
    "runs/2026-11-21_rfd4-backbones-BURIED_job21_32designs_seed201",
) + "/manifest.jsonl"


def find(work: Path, arm: str, did: str):
    hits = glob.glob(f"{work}/out_{arm}/**/predictions/{did}/{did}_model_0.cif",
                     recursive=True)
    conf = glob.glob(
        f"{work}/out_{arm}/**/predictions/{did}/confidence_{did}_model_0.json",
        recursive=True)
    return (Path(hits[0]) if hits else None,
            Path(conf[0]) if conf else None)


def main(root: Path, work: Path, did: str) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    mf = Path(RUN)
    if not mf.is_absolute():
        mf = root / RUN
    rec = next((r for r in read_manifest(mf) if r.design_id == did), None)
    if rec is None:
        raise SystemExit(f"design {did} not in manifest")
    cat = rec.catalytic_residues
    wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat["wall"]

    results = {}
    coords = {}
    for arm in ("with_msa", "no_msa"):
        cif, conf = find(work, arm, did)
        if cif is None:
            print(f"{arm}: NO PREDICTION")
            continue
        s = get_structure(CIFFile.read(cif), model=1)
        copies = g.iter_ligand_copies(s, amap)
        row = {}
        if conf and conf.exists():
            c = json.loads(conf.read_text())
            row |= {"plddt": c.get("complex_plddt"), "iptm": c.get("iptm"),
                    "ligand_iptm": c.get("ligand_iptm"), "ptm": c.get("ptm")}
        if copies:
            lc, lr = copies[0]
            kw = dict(lig_chain=lc, lig_resid=lr, chain="A")
            try:
                row |= {
                    "base": g.base_to_abstraction(s, amap, base_resid=cat["base"], **kw),
                    "donor": g.donor_to_oxyanion(s, amap, donor_resid=cat["donor"], **kw),
                    "wall": g.wall_stacking(s, amap, wall_resid=wall, **kw),
                }
                row |= g.ligand_burial(s, amap, lig_chain=lc, lig_resid=lr)
            except g.MetricUnavailable as exc:
                row["error"] = str(exc)
        results[arm] = row
        prot = s[(~s.hetero) & (s.chain_id == "A")]
        ca = prot[prot.atom_name == "CA"]
        coords[arm] = ca.coord

    print(f"\ndesign {did}  (length {rec.length})")
    keys = ["plddt", "iptm", "ligand_iptm", "base", "donor", "wall",
            "rel_sasa", "enclosure", "n_contacts"]
    print(f"{'metric':13} {'with_msa':>12} {'no_msa':>12} {'delta':>10}")
    print("-" * 50)
    for k in keys:
        a = results.get("with_msa", {}).get(k)
        b = results.get("no_msa", {}).get(k)
        if isinstance(a, float) and isinstance(b, float):
            print(f"{k:13} {a:>12.3f} {b:>12.3f} {a - b:>+10.3f}")
        else:
            print(f"{k:13} {str(a):>12} {str(b):>12} {'-':>10}")

    if len(coords) == 2 and coords["with_msa"].shape == coords["no_msa"].shape:
        rms = g.kabsch_rmsd(coords["no_msa"], coords["with_msa"])
        print(f"\nCA RMSD between the two predictions: {rms:.3f} A")
        print("  (0.00 would mean the MSA changed nothing at all)")

    print("\nVERDICT:", end=" ")
    a = results.get("with_msa", {})
    b = results.get("no_msa", {})
    if "base" in a and "base" in b:
        moved = max(abs(a[k] - b[k]) for k in ("base", "donor", "wall"))
        if moved < 0.1:
            print("MSA made no meaningful difference (site geometry identical).")
        elif a["donor"] < b["donor"] - 0.5:
            print("MSA IMPROVED the donor geometry.")
        elif a["donor"] > b["donor"] + 0.5:
            print("MSA WORSENED the donor geometry.")
        else:
            print(f"MSA changed the site geometry by up to {moved:.2f} A, "
                  "no clear improvement.")
    else:
        print("insufficient data to compare.")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3])
