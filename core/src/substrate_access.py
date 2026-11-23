"""Can a substrate reach the active site? Bottleneck probe radius to the pocket.

The stage-1 gate rewards a BURIED ligand (rel_sasa ~0). A fully enclosed pocket
has no way in for the substrate, so burial and access pull against each other.
This measures access directly.

Method (ligand REMOVED so the pocket is empty, as it would be before binding):
  1. Grid the protein (heavy atoms, hetero removed) at 0.6 A spacing.
  2. clearance(x) = min over atoms of (distance - vdW radius).
  3. For probe radius r, free space = clearance > r. Label connected
     components (6-connectivity).
  4. The EXTERIOR is the component containing a grid corner (8 A padding).
  5. The SITE is the set of free grid points within 2.2 A of the removed
     ligand's heavy atoms. It is accessible at radius r if any site point lies
     in the exterior component.
  6. bottleneck = the largest r (from the tested list) that is accessible.
     0 means not even a 0.8 A probe gets in; 'site_closed_at' notes the pocket
     itself being too small for the probe, which is not a tunnel problem.

Reference sizes: water ~1.4 A. The substrate (5-nitrobenzisoxazole) and 6NT are
planar aromatics ~1.7-1.9 A half-thickness, so a ~1.5 A sphere is a lenient
stand-in and 1.8-2.0 A is closer to a rigid pass-through. Static structures
under-estimate access (side chains breathe), so interpret against the crystals.

Usage: python -I substrate_access.py <project root> <out csv>
"""

from __future__ import annotations

import csv
import glob
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from scipy import ndimage
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from base_environment import RUNS  # noqa: E402

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
RADII = [0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
SPACING, PAD, SITE_R = 0.6, 8.0, 2.2


def access(struct, lig_name="6NT", chain="A") -> dict | None:
    prot = struct[(~struct.hetero) & (struct.chain_id == chain) & (struct.element != "H")]
    lig = struct[(struct.res_name == lig_name) & (struct.element != "H")]
    if lig.array_length() == 0 or prot.array_length() == 0:
        return None
    xyz = prot.coord
    lo, hi = xyz.min(0) - PAD, xyz.max(0) + PAD
    axes = [np.arange(lo[i], hi[i] + SPACING, SPACING) for i in range(3)]
    shape = tuple(len(a) for a in axes)
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)

    clear = np.full(len(grid), np.inf)
    elems = np.array([str(e) for e in prot.element])
    for el, rad in VDW.items():
        m = elems == el
        if m.any():
            d, _ = cKDTree(xyz[m]).query(grid, distance_upper_bound=8.0)
            clear = np.minimum(clear, d - rad)
    clear = clear.reshape(shape)

    # site grid points: near the (removed) ligand
    near = cKDTree(lig.coord).query(grid, distance_upper_bound=SITE_R)[0] < SITE_R
    near = near.reshape(shape)

    best, closed_at = 0.0, None
    res = {}
    for r in RADII:
        free = clear > r
        site = free & near
        if not site.any():
            res[r] = "closed"
            if closed_at is None:
                closed_at = r
            continue
        lab, _ = ndimage.label(free)
        ext = lab[0, 0, 0]
        ok = bool(np.isin(lab[site], ext).any()) if ext else False
        res[r] = ok
        if ok:
            best = r
    return {"bottleneck": best, "site_closed_at": closed_at, "by_radius": res}


def structure_rows(root: Path):
    # crystals: calibration set
    for pdb, chain, lig in (("7K4Z", "A", "6NT"), ("4BS0", "A", "6NT"),
                            ("7K4X", "A", "6NT"), ("8RD5", "A", "6NT"),
                            ("7K4Q", "A", "6NT"), ("6SS1", "A", "6NT")):
        s = get_structure(CIFFile.read(root / "benchmark" / "pdb" / f"{pdb}.cif"), model=1)
        yield "crystal", pdb, s
    # designs: refold structures of every arm
    for arm, (_, raw, _) in RUNS.items():
        meta = json.loads((root / raw / "designs.json").read_text())
        for name in meta:
            h = glob.glob(f"{root}/{raw}/**/predictions/{name}/{name}_model_0.cif",
                          recursive=True)
            if h:
                yield arm, name, get_structure(CIFFile.read(h[0]), model=1)


def main(root: Path, out_csv: Path, only: set[str] | None = None) -> None:
    rows = []
    for arm, name, s in structure_rows(root):
        if only is not None and arm != "crystal" and name not in only:
            continue
        a = access(s)
        if a is None:
            continue
        rows.append({"arm": arm, "design_id": name, "bottleneck": a["bottleneck"],
                     "site_closed_at": a["site_closed_at"] or "",
                     **{f"open_{r}": a["by_radius"][r] for r in RADII}})
        print(f"{arm:13}{name:18} bottleneck={a['bottleneck']:.1f}"
              f"  {''.join('Y' if a['by_radius'][r] is True else ('c' if a['by_radius'][r] == 'closed' else '.') for r in RADII)}",
              flush=True)
    cols = list(rows[0])
    with open(out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in cols) + "\n")
    print(f"\nwrote {out_csv}   (columns Y/./c = open / blocked / pocket too small at "
          f"{RADII})")


if __name__ == "__main__":
    only = None
    if len(sys.argv) > 3:
        only = set(Path(sys.argv[3]).read_text().split())
    main(Path(sys.argv[1]), Path(sys.argv[2]), only)
