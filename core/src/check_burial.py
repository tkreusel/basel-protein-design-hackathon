"""Is the ligand actually buried in a pocket, or stuck on the surface?

The motif gate in extract_designs.py checks catalytic DISTANCES only. A design
can satisfy every distance while the ligand sits on the outside of the protein
with the three catalytic residues draped over it -- geometrically "correct",
catalytically useless. This script measures enclosure directly and compares
every design against the HG3.17 crystal reference.

Metrics:
  rel_sasa      ligand SASA in the complex / ligand SASA alone.
                0.0 = fully buried, 1.0 = fully exposed.
  n_contacts    protein heavy atoms within 5 A of any ligand heavy atom.
  n_shell       protein heavy atoms within 8 A -- how much pocket surrounds it.
  enclosure     fraction of directions from the ligand centroid that hit
                protein within 12 A. A real pocket surrounds the ligand; a
                surface site only blocks one side. This is the metric that
                actually distinguishes "in a pocket" from "against a wall".

Usage: python -I check_burial.py <project root> <cif or dir> [...]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from biotite.structure import sasa
from biotite.structure.io.pdbx import CIFFile, get_structure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402


def fibonacci_directions(n: int = 256) -> np.ndarray:
    """Evenly distributed unit vectors on a sphere."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]
    )


def enclosure_fraction(lig_xyz, prot_xyz, max_dist=12.0, tol=2.8, n_dir=256):
    """Fraction of ray directions from the ligand centroid that hit protein.

    A ray "hits" if any protein heavy atom lies within `tol` of the ray within
    `max_dist`. Buried pockets approach 1.0; a ligand on a flat surface sits
    near 0.5; a fully exposed one near 0.
    """
    centre = lig_xyz.mean(axis=0)
    rel = prot_xyz - centre
    dist = np.linalg.norm(rel, axis=1)
    near = rel[dist <= max_dist]
    if len(near) == 0:
        return 0.0
    dirs = fibonacci_directions(n_dir)
    # projection of each atom onto each direction
    proj = near @ dirs.T                       # (n_atoms, n_dir)
    perp2 = (near ** 2).sum(axis=1)[:, None] - proj ** 2
    hit = (proj > 0) & (perp2 <= tol ** 2)
    return float(hit.any(axis=0).mean())


def analyse(path: Path, amap: g.AtomMap) -> dict | None:
    struct = get_structure(CIFFile.read(path), model=1)
    copies = g.iter_ligand_copies(struct, amap)
    if not copies:
        return None
    lchain, lresid = copies[0]

    lig_mask = (
        (struct.chain_id == lchain)
        & (struct.res_id == lresid)
        & (struct.res_name == amap.ligand_ccd_id)
        & (struct.element != "H")
    )
    prot_mask = (~struct.hetero) & (struct.element != "H")
    if not lig_mask.any() or not prot_mask.any():
        return None

    lig_xyz = struct.coord[lig_mask]
    prot_xyz = struct.coord[prot_mask]

    # SASA of the ligand alone vs in the complex.
    complex_atoms = struct[lig_mask | prot_mask]
    lig_in_complex = np.concatenate(
        [np.zeros(int(prot_mask.sum()), bool), np.ones(int(lig_mask.sum()), bool)]
    ) if False else None  # ordering handled below

    # Rebuild explicitly so index order is known: protein first, then ligand.
    ordered = struct[prot_mask] + struct[lig_mask]
    n_prot = int(prot_mask.sum())
    sasa_complex = sasa(ordered, ignore_ions=False)
    lig_sasa_complex = float(np.nansum(sasa_complex[n_prot:]))

    lig_only = struct[lig_mask]
    lig_sasa_free = float(np.nansum(sasa(lig_only, ignore_ions=False)))
    rel = lig_sasa_complex / lig_sasa_free if lig_sasa_free > 0 else float("nan")

    d = np.linalg.norm(prot_xyz[:, None, :] - lig_xyz[None, :, :], axis=-1).min(axis=1)
    return {
        "rel_sasa": rel,
        "n_contacts": int((d < 5.0).sum()),
        "n_shell": int((d < 8.0).sum()),
        "enclosure": enclosure_fraction(lig_xyz, prot_xyz),
        "n_res": len(set(struct.res_id[prot_mask])),
    }


def main(root: Path, targets: list[Path]) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")

    files: list[Path] = []
    for t in targets:
        if t.is_dir():
            files += sorted(t.rglob("*_generated.cif")) or sorted(t.rglob("*.cif"))
        else:
            files.append(t)

    print("rel_sasa: 0=buried 1=exposed | enclosure: 1=surrounded 0=open\n")
    hdr = f"{'structure':46} {'res':>4} {'rel_sasa':>9} {'contacts':>9} {'shell':>6} {'enclosure':>10}"
    print(hdr)
    print("-" * len(hdr))

    rows = []
    for f in files:
        try:
            r = analyse(f, amap)
        except Exception as exc:  # noqa: BLE001
            print(f"{f.name[:46]:46} ERROR {type(exc).__name__}: {exc}")
            continue
        if r is None:
            print(f"{f.name[:46]:46} (no ligand)")
            continue
        rows.append((f, r))
        print(f"{f.name[:46]:46} {r['n_res']:>4} {r['rel_sasa']:>9.3f} "
              f"{r['n_contacts']:>9} {r['n_shell']:>6} {r['enclosure']:>10.3f}")

    if len(rows) > 1:
        arr = {k: np.array([r[k] for _, r in rows], float)
               for k in ("rel_sasa", "n_contacts", "n_shell", "enclosure")}
        print("-" * len(hdr))
        print(f"{'MEAN':46} {'':>4} {arr['rel_sasa'].mean():>9.3f} "
              f"{arr['n_contacts'].mean():>9.1f} {arr['n_shell'].mean():>6.1f} "
              f"{arr['enclosure'].mean():>10.3f}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), [Path(p) for p in sys.argv[2:]])
