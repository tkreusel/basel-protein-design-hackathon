"""Geometric metrics for Kemp eliminase designs (plan 3, 6.1).

Design rules, from the plan:
  * Pure functions. Every metric takes coordinates (or a structure plus a role
    spec) and returns a float. No I/O, no global state, no model calls.
  * Atom roles are READ FROM configs/atom_map.yaml, never hard-coded here
    (plan 6.1). Constants in this module are config KEYS, not atom names.
  * Fail loud (plan 0.2). A metric that cannot be computed raises
    MetricUnavailable with a reason; callers record that reason in the
    DesignRecord rather than substituting a silent default. There is exactly
    one deliberate exception, documented at `donor_to_oxyanion`: a cell with
    no oxyanion donor is a valid design choice, not an error.

Because no AF3 container exists on this cluster, RFD4's motif-RMSD metric is
unavailable (docs/TOOL_NOTES.md), so `motif_heavy_atom_rmsd` here is the
pipeline's only motif-fidelity check. It is load-bearing.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import yaml

__all__ = [
    "MetricUnavailable",
    "AtomMap",
    "load_atom_map",
    "load_thresholds",
    "centroid",
    "distance",
    "min_pairwise_distance",
    "kabsch_rmsd",
    "base_to_abstraction",
    "donor_to_oxyanion",
    "wall_stacking",
    "motif_heavy_atom_rmsd",
    "min_protein_ligand_contact",
    "iter_ligand_copies",
]


class MetricUnavailable(RuntimeError):
    """A metric could not be computed. Carries why, for the DesignRecord."""


# ---------------------------------------------------------------------------
# config loading
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AtomMap:
    """The atom-role mapping from configs/atom_map.yaml.

    Thin, validated accessor. Construct via `load_atom_map`; the validation in
    `__post_init__` is what stops a malformed config from silently producing
    plausible-looking numbers for the wrong atoms.
    """

    ligand_ccd_id: str
    abstraction_centre: str
    oxyanion_centre: str
    nitro_nitrogen: str
    nitro_oxygens: tuple[str, ...]
    benzo_ring: tuple[str, ...]
    base_tip_atoms: tuple[str, ...]
    base_proximal_atom: str
    base_all_atoms: tuple[str, ...]
    donating_atom_by_resname: dict[str, str]
    wall_ring_atoms: tuple[str, ...]
    raw: dict

    def __post_init__(self) -> None:
        if self.abstraction_centre == self.oxyanion_centre:
            raise ValueError(
                "atom_map: abstraction_centre and oxyanion_centre are the same "
                f"atom ({self.abstraction_centre!r}); the mapping is wrong"
            )
        for field, value in (
            ("abstraction_centre", self.abstraction_centre),
            ("oxyanion_centre", self.oxyanion_centre),
            ("base_proximal_atom", self.base_proximal_atom),
        ):
            if not value:
                raise ValueError(f"atom_map: {field} is empty")
        if self.base_proximal_atom not in self.base_tip_atoms:
            raise ValueError(
                f"atom_map: base proximal_atom {self.base_proximal_atom!r} is "
                f"not among tip_atoms {self.base_tip_atoms}"
            )


def load_atom_map(path: str | Path) -> AtomMap:
    """Load and validate configs/atom_map.yaml."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    try:
        roles = cfg["roles"]
        ref = cfg["reference_scaffold"]
        return AtomMap(
            ligand_ccd_id=cfg["ligand"]["ccd_id"],
            abstraction_centre=roles["abstraction_centre"]["atom_6NT"],
            oxyanion_centre=roles["oxyanion_centre"]["atom_6NT"],
            nitro_nitrogen=roles["nitro_group"]["nitrogen_6NT"],
            nitro_oxygens=tuple(roles["nitro_group"]["oxygens_6NT"]),
            benzo_ring=tuple(roles["benzo_ring"]["atoms_6NT"]),
            base_tip_atoms=tuple(ref["base"]["tip_atoms"]),
            base_proximal_atom=ref["base"]["proximal_atom"],
            base_all_atoms=tuple(ref["base"]["all_atoms"]),
            donating_atom_by_resname=dict(ref["donor"]["donating_atom_by_resname"]),
            wall_ring_atoms=tuple(ref["wall"]["ring_atoms"]),
            raw=cfg,
        )
    except KeyError as exc:
        raise ValueError(f"atom_map is missing required key: {exc}") from exc


def load_thresholds(path: str | Path) -> dict:
    """Load configs/thresholds.yaml verbatim. gate.py interprets it."""
    with open(path) as fh:
        return yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# primitives
# ---------------------------------------------------------------------------


def centroid(coords: np.ndarray) -> np.ndarray:
    """Mean position of an (N, 3) coordinate array."""
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"expected (N, 3) coordinates, got {coords.shape}")
    if coords.shape[0] == 0:
        raise MetricUnavailable("centroid of an empty atom set")
    return coords.mean(axis=0)


def distance(a: np.ndarray, b: np.ndarray) -> float:
    """Euclidean distance between two 3-vectors."""
    return float(np.linalg.norm(np.asarray(a, float) - np.asarray(b, float)))


def min_pairwise_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Smallest distance between any atom of `a` and any atom of `b`."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        raise MetricUnavailable("min_pairwise_distance over an empty atom set")
    return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min())


def kabsch_rmsd(mobile: np.ndarray, target: np.ndarray) -> float:
    """RMSD after optimal superposition (Kabsch). Both (N, 3), 1:1 ordered.

    Used for the motif-fidelity check, so it must superpose rather than compare
    in place: RFD4/RFD3 place the motif in an arbitrary global frame.
    """
    mobile = np.asarray(mobile, float)
    target = np.asarray(target, float)
    if mobile.shape != target.shape:
        raise MetricUnavailable(
            f"kabsch_rmsd shape mismatch: {mobile.shape} vs {target.shape}"
        )
    if mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError(f"expected (N, 3) coordinates, got {mobile.shape}")
    if mobile.shape[0] < 3:
        raise MetricUnavailable(
            f"kabsch_rmsd needs >=3 atoms to superpose, got {mobile.shape[0]}"
        )
    p = mobile - mobile.mean(axis=0)
    q = target - target.mean(axis=0)
    # Optimal rotation via SVD of the covariance, with a reflection guard.
    v, _, wt = np.linalg.svd(p.T @ q)
    if np.linalg.det(v @ wt) < 0.0:
        v[:, -1] *= -1.0
    p_rot = p @ (v @ wt)
    return float(np.sqrt(((p_rot - q) ** 2).sum(axis=1).mean()))


# ---------------------------------------------------------------------------
# atom selection against a biotite AtomArray
# ---------------------------------------------------------------------------


def _select(struct, *, chain=None, resid=None, names=None, resname=None, hetero=None):
    """Boolean mask over a biotite AtomArray. Kept private and explicit."""
    mask = np.ones(struct.array_length(), dtype=bool)
    if chain is not None:
        mask &= struct.chain_id == chain
    if resid is not None:
        mask &= struct.res_id == resid
    if resname is not None:
        mask &= struct.res_name == resname
    if names is not None:
        mask &= np.isin(struct.atom_name, list(names))
    if hetero is not None:
        mask &= struct.hetero == hetero
    return mask


def _coords(struct, mask, what: str) -> np.ndarray:
    if not mask.any():
        raise MetricUnavailable(f"no atoms found for {what}")
    return struct.coord[mask]


def _one_coord(struct, mask, what: str) -> np.ndarray:
    coords = _coords(struct, mask, what)
    if len(coords) > 1:
        # Altlocs or a duplicated atom name: refuse to pick arbitrarily.
        raise MetricUnavailable(
            f"{what} matched {len(coords)} atoms; expected exactly 1 "
            "(altloc or duplicate atom name?)"
        )
    return coords[0]


def _ligand_atom(struct, amap: AtomMap, atom_name: str, chain, resid) -> np.ndarray:
    mask = _select(struct, chain=chain, resid=resid, names=[atom_name])
    mask &= struct.res_name == amap.ligand_ccd_id
    return _one_coord(struct, mask, f"ligand {amap.ligand_ccd_id} atom {atom_name}")


# ---------------------------------------------------------------------------
# the metrics
# ---------------------------------------------------------------------------


def base_to_abstraction(
    struct, amap: AtomMap, *, base_resid: int, lig_chain, lig_resid, chain="A"
) -> float:
    """Catalytic base carboxylate -> ligand abstraction centre, in Angstrom.

    Reports the SHORTEST distance over the base's tip atoms rather than the
    configured proximal atom: in a de novo design the carboxylate may be
    flipped relative to the HG3 reference, and OD1/OD2 are chemically
    equivalent, so insisting on OD2 would penalise an equally good geometry.

    Calibration note: across the benchmark this metric does NOT track activity
    (HG3, the weakest variant, has the shortest distance). Necessary but not
    discriminating -- see docs/MILESTONE1.md before weighting it.
    """
    target = _ligand_atom(struct, amap, amap.abstraction_centre, lig_chain, lig_resid)
    tips = _coords(
        struct,
        _select(
            struct, chain=chain, resid=base_resid, names=amap.base_tip_atoms,
            hetero=False,
        ),
        f"base residue {chain}{base_resid} tip atoms {amap.base_tip_atoms}",
    )
    return float(np.linalg.norm(tips - target, axis=1).min())


def donor_to_oxyanion(
    struct,
    amap: AtomMap,
    *,
    donor_resid: int | None,
    lig_chain,
    lig_resid,
    chain="A",
    absent_sentinel: float = 99.0,
) -> float:
    """Oxyanion donor tip -> ligand oxyanion centre, in Angstrom.

    THE discriminating metric. Benchmark: Gln50 2.96-3.10 A (active) vs Lys50
    6.05-7.42 A (inactive) vs Ala50 no tip.

    Deliberate exception to fail-loud: a theozyme cell with no oxyanion donor
    (plan cell family 'c', water-mediated) is a valid design, and an ablated
    residue genuinely has no donating atom. Both return `absent_sentinel` so
    the design ranks last on this metric instead of crashing the batch. A
    donor that is *specified but missing from the structure* still raises.
    """
    if donor_resid is None:
        return absent_sentinel

    resnames = set(
        struct.res_name[_select(struct, chain=chain, resid=donor_resid, hetero=False)]
    )
    if not resnames:
        raise MetricUnavailable(
            f"donor residue {chain}{donor_resid} is absent from the structure"
        )
    if len(resnames) > 1:
        raise MetricUnavailable(
            f"donor residue {chain}{donor_resid} has multiple residue names: {resnames}"
        )
    resname = str(resnames.pop())

    donating = amap.donating_atom_by_resname.get(resname)
    if not donating:
        # A residue type with no donating atom (e.g. ALA): ablated by design.
        return absent_sentinel

    target = _ligand_atom(struct, amap, amap.oxyanion_centre, lig_chain, lig_resid)
    tip = _one_coord(
        struct,
        _select(struct, chain=chain, resid=donor_resid, names=[donating], hetero=False),
        f"donor {resname}{donor_resid} atom {donating}",
    )
    return distance(tip, target)


def wall_stacking(
    struct,
    amap: AtomMap,
    *,
    wall_resid: int,
    lig_chain,
    lig_resid,
    chain="A",
    ring_atoms: Sequence[str] | None = None,
) -> float:
    """Pocket-wall ring centroid -> ligand benzo-ring centroid, in Angstrom.

    `ring_atoms` overrides the configured Trp ring for a design that uses a
    different aromatic wall residue (Phe, Tyr, His).
    """
    names = tuple(ring_atoms) if ring_atoms is not None else amap.wall_ring_atoms
    wall = _coords(
        struct,
        _select(struct, chain=chain, resid=wall_resid, names=names, hetero=False),
        f"wall residue {chain}{wall_resid} ring atoms",
    )
    lig_mask = _select(struct, chain=lig_chain, resid=lig_resid, names=amap.benzo_ring)
    lig_mask &= struct.res_name == amap.ligand_ccd_id
    lig = _coords(struct, lig_mask, f"ligand benzo ring {amap.benzo_ring}")
    return distance(centroid(wall), centroid(lig))


def motif_heavy_atom_rmsd(design_coords: np.ndarray, motif_coords: np.ndarray) -> float:
    """Catalytic heavy-atom RMSD of a design's motif to the input motif.

    Plan 6.2 acceptance: < 1.5 A. The caller must pass the two coordinate sets
    in the SAME atom order; a mismatch there is a silent correctness bug, so
    build both arrays from one shared ordered atom list.
    """
    return kabsch_rmsd(design_coords, motif_coords)


def min_protein_ligand_contact(struct, amap: AtomMap, *, lig_chain, lig_resid) -> float:
    """Closest protein-heavy-atom to ligand-heavy-atom distance, in Angstrom.

    Plan 4.3 / 6.2 [lit]: reject if any pair is closer than 1.5 A. Hydrogens
    are excluded so the metric means the same thing on structures with and
    without riding hydrogens -- the benchmark crystals have them, designs will
    not.
    """
    lig_mask = _select(struct, chain=lig_chain, resid=lig_resid)
    lig_mask &= struct.res_name == amap.ligand_ccd_id
    lig_mask &= struct.element != "H"
    lig = _coords(struct, lig_mask, f"ligand {amap.ligand_ccd_id}")

    prot_mask = (~struct.hetero) & (struct.element != "H")
    prot = _coords(struct, prot_mask, "protein heavy atoms")
    return min_pairwise_distance(prot, lig)


def ligand_burial(struct, amap: AtomMap, *, lig_chain, lig_resid) -> dict[str, float]:
    """Is the ligand enclosed in a pocket, or stuck on the surface?

    The distance metrics above are necessary but NOT sufficient: a design can
    satisfy every catalytic distance while the ligand sits on the outside of
    the protein with the three catalytic residues draped over it. Measured on
    a real 32-design batch, mean ligand rel_sasa was 0.29 against 0.008 for the
    HG3.17/HG3.7 crystals -- most designs had no pocket at all, and a
    distance-only gate passed all of them.

    Returns:
      rel_sasa    ligand SASA in complex / ligand SASA alone.
                  0 = fully buried, 1 = fully exposed. Crystals: 0.002-0.011.
      enclosure   fraction of ray directions from the ligand centroid that
                  hit protein within 12 A. Distinguishes "inside a pocket"
                  (~1.0) from "against a flat surface" (~0.5).
                  Crystals: 0.99-1.00.
      n_contacts  protein heavy atoms within 5 A. Crystals: 65-68.
      n_shell     protein heavy atoms within 8 A. Crystals: 200-211.
    """
    from biotite.structure import sasa as _sasa

    lig_mask = _select(struct, chain=lig_chain, resid=lig_resid)
    lig_mask &= struct.res_name == amap.ligand_ccd_id
    lig_mask &= struct.element != "H"
    prot_mask = (~struct.hetero) & (struct.element != "H")
    if not lig_mask.any():
        raise MetricUnavailable(f"ligand {amap.ligand_ccd_id} not found for burial")
    if not prot_mask.any():
        raise MetricUnavailable("no protein heavy atoms for burial")

    # Protein first, then ligand, so slice indices are known.
    ordered = struct[prot_mask] + struct[lig_mask]
    n_prot = int(prot_mask.sum())
    lig_in_complex = float(np.nansum(_sasa(ordered, ignore_ions=False)[n_prot:]))
    lig_free = float(np.nansum(_sasa(struct[lig_mask], ignore_ions=False)))
    rel = lig_in_complex / lig_free if lig_free > 0 else float("nan")

    lig_xyz = struct.coord[lig_mask]
    prot_xyz = struct.coord[prot_mask]
    dmin = np.linalg.norm(
        prot_xyz[:, None, :] - lig_xyz[None, :, :], axis=-1
    ).min(axis=1)

    # Ray-cast enclosure from the ligand centroid.
    centre = lig_xyz.mean(axis=0)
    rel_xyz = prot_xyz - centre
    dist = np.linalg.norm(rel_xyz, axis=1)
    near = rel_xyz[dist <= 12.0]
    if len(near) == 0:
        enclosure = 0.0
    else:
        n_dir = 256
        i = np.arange(n_dir) + 0.5
        phi = np.arccos(1 - 2 * i / n_dir)
        theta = np.pi * (1 + 5 ** 0.5) * i
        dirs = np.column_stack(
            [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)]
        )
        proj = near @ dirs.T
        perp2 = (near ** 2).sum(axis=1)[:, None] - proj ** 2
        enclosure = float(((proj > 0) & (perp2 <= 2.8 ** 2)).any(axis=0).mean())

    return {
        "rel_sasa": float(rel),
        "enclosure": enclosure,
        "n_contacts": float((dmin < 5.0).sum()),
        "n_shell": float((dmin < 8.0).sum()),
    }


def iter_ligand_copies(struct, amap: AtomMap) -> list[tuple[str, int]]:
    """Every (chain_id, res_id) of the ligand present in the structure."""
    mask = struct.res_name == amap.ligand_ccd_id
    return sorted(
        {(str(c), int(r)) for c, r in zip(struct.chain_id[mask], struct.res_id[mask])}
    )
