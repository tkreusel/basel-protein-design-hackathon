"""Score the organisers' supplied Kemp motif with our own calibrated metrics.

Answers one question: is the provided inputs/HG3_17_kemp_motif.pdb a
catalytically competent active site by the thresholds we calibrated from the
crystal series (configs/thresholds.yaml)? If yes, it is a ready-to-scaffold
theozyme and stage 0 does not have to build one from scratch.

Usage: python -I check_motif.py <project root> <motif pdb>
"""

from __future__ import annotations

import sys
from pathlib import Path

from biotite.structure.io.pdb import PDBFile

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402

BASE_RESID, DONOR_RESID, WALL_RESID = 127, 50, 44


def main(root: Path, motif_pdb: Path) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    thr = g.load_thresholds(root / "configs" / "thresholds.yaml")["geometry"]

    struct = PDBFile.read(motif_pdb).get_structure(model=1)
    copies = g.iter_ligand_copies(struct, amap)
    print(f"motif: {motif_pdb}")
    print(f"atoms: {struct.array_length()}   ligand copies: {copies}")

    lc, lr = copies[0]
    kw = dict(lig_chain=lc, lig_resid=lr, chain="A")

    base = g.base_to_abstraction(struct, amap, base_resid=BASE_RESID, **kw)
    donor = g.donor_to_oxyanion(struct, amap, donor_resid=DONOR_RESID, **kw)
    wall = g.wall_stacking(struct, amap, wall_resid=WALL_RESID, **kw)

    rows = [
        ("base -> abstraction (N3)", base, thr["base_to_abstraction"]),
        ("donor -> oxyanion  (N1)", donor, thr["donor_to_oxyanion"]),
        ("wall stacking", wall, thr["wall_stacking"]),
    ]
    print(f"\n{'metric':26} {'value':>7} {'ideal':>12} {'pass_max':>9}  verdict")
    print("-" * 70)
    all_pass = True
    for name, value, spec in rows:
        lo, hi = spec["ideal"]
        pmax = spec["pass_max"]
        ok = value <= pmax
        all_pass &= ok
        inside = "IDEAL" if lo <= value <= hi else ("pass" if ok else "FAIL")
        print(f"{name:26} {value:7.2f} {f'{lo}-{hi}':>12} {pmax:9.2f}  {inside}")

    print(f"\nverdict: {'CATALYTICALLY COMPETENT' if all_pass else 'FAILS A GATE'}")
    print("\nFor comparison, the crystal series (chain A):")
    print("  HG3.17 7K4Z : base 2.51  donor 3.10  wall 5.54   (126,000 M-1s-1)")
    print("  HG3.7  7K4X : base 2.49  donor 2.96  wall 5.01   ( 27,000)")
    print("  HG3    7K4Q : base 2.25  donor 6.05  wall 5.02   (    146)")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
