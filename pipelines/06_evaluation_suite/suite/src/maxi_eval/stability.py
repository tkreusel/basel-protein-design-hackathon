"""Foldability / pocket-quality proxies that tmol does not provide natively.

These are APPROXIMATIONS of Rosetta filters (BuriedUnsatHbonds, holes/packstat). They are good for
ranking and for catching gross problems; they are not numerically comparable to the Rosetta filters.
Inputs are the tmol-protonated protein + ligand AtomArrays (tmol_core.pose_to_replicate_arrays).
"""
from __future__ import annotations

import numpy as np

# acceptor-capable atoms (protein) by residue; backbone O / OXT always
ACCEPTORS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"], "ASN": ["OD1"], "GLN": ["OE1"], "SER": ["OG"],
             "THR": ["OG1"], "TYR": ["OH"], "HIS": ["ND1", "NE2"], "MET": [], "*": ["O", "OXT"]}
LIG_ACCEPTORS = {"N1", "N2", "O11", "O21", "O1"}   # ring N / nitro O / isoxazole-or-phenolate O
BURIED_SASA = 1.0   # A^2 : heavy-atom SASA below this counts as buried
HB_H_A = 2.5        # A   : H...A cutoff
HB_ANGLE = 110.0    # deg : D-H...A minimum


def _angle(a, b, c):
    v1, v2 = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9), -1, 1))))


def buried_unsatisfied_polars(prot, lig, pocket_resids: list[int]) -> dict:
    """Approximate BUNS in the pocket shell + on the ligand."""
    from biotite.structure import sasa

    both = prot + lig
    heavy_mask = both.element != "H"
    s = np.full(len(both), np.nan)
    from .geometry import vdw_radii
    hv = both[heavy_mask]
    s[heavy_mask] = sasa(hv, vdw_radii=vdw_radii(hv), ignore_ions=False)
    is_lig = np.concatenate([np.zeros(len(prot), bool), np.ones(len(lig), bool)])
    H = np.where(both.element == "H")[0]
    polar = np.where(np.isin(both.element, ["N", "O"]))[0]
    # donors: polar heavy atoms with an H within 1.2 A
    hx = both.coord[H]
    donor_H: dict[int, list[int]] = {}
    for i in polar:
        d = np.linalg.norm(hx - both.coord[i], axis=1)
        donor_H[i] = [int(H[k]) for k in np.where(d < 1.2)[0]]

    def is_acceptor(i):
        if is_lig[i]:
            return str(both.atom_name[i]) in LIG_ACCEPTORS
        rn, an = str(both.res_name[i]), str(both.atom_name[i])
        return an in ACCEPTORS.get(rn, []) or an in ACCEPTORS["*"]

    acc_idx = [i for i in polar if is_acceptor(i)]
    pocket = set(pocket_resids)
    unsat, n_buried_polar = [], 0
    for i in polar:
        in_scope = is_lig[i] or int(both.res_id[i]) in pocket
        if not in_scope or not np.isfinite(s[i]) or s[i] >= BURIED_SASA:
            continue
        hs = donor_H[i]
        acc = is_acceptor(i)
        if not hs and not acc:
            continue
        n_buried_polar += 1
        ok = False
        for h in hs:                                   # as donor
            for a in acc_idx:
                if a == i:
                    continue
                if np.linalg.norm(both.coord[h] - both.coord[a]) <= HB_H_A and _angle(both.coord[i], both.coord[h], both.coord[a]) >= HB_ANGLE:
                    ok = True
                    break
            if ok:
                break
        if not ok and acc:                             # as acceptor
            for j in polar:
                if j == i:
                    continue
                for h in donor_H[j]:
                    if np.linalg.norm(both.coord[h] - both.coord[i]) <= HB_H_A and _angle(both.coord[j], both.coord[h], both.coord[i]) >= HB_ANGLE:
                        ok = True
                        break
                if ok:
                    break
        if not ok:
            unsat.append(f"{'LIG' if is_lig[i] else str(both.res_name[i]) + str(int(both.res_id[i]))}:{both.atom_name[i]}")
    return {"buns.n_buried_polar": float(n_buried_polar), "buns.n_unsat": float(len(unsat)),
            "buns.n_unsat_ligand": float(sum(u.startswith("LIG") for u in unsat)), "buns.atoms": ";".join(unsat)}


# vdW radii (A) for the void probe
_VDW = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "H": 1.1}


def pocket_void_volume(prot, lig, shell: float = 6.0, probe: float = 1.2, spacing: float = 0.6,
                       min_enclosure: float = 0.85) -> dict:
    """Empty-but-enclosed volume next to the ligand (proxy for holes / poor packing).

    Grid points within `shell` of the ligand where a `probe`-radius sphere fits between all atoms
    (ligand included) and which are enclosed (>= min_enclosure of 64 rays hit an atom within 10 A).
    """
    from scipy.spatial import cKDTree

    atoms = prot + lig
    xyz = atoms.coord
    rad = np.array([_VDW.get(e, 1.7) for e in atoms.element])
    lo, hi = lig.coord.min(0) - shell, lig.coord.max(0) + shell
    g = np.mgrid[lo[0]:hi[0]:spacing, lo[1]:hi[1]:spacing, lo[2]:hi[2]:spacing].reshape(3, -1).T
    near_lig = cKDTree(lig.coord).query(g)[0] <= shell
    g = g[near_lig]
    tree = cKDTree(xyz)
    cand = tree.query_ball_point(g, r=4.0)
    free = np.array([len(c) == 0 or np.all(np.linalg.norm(xyz[c] - p, axis=1) - rad[c] > probe) for p, c in zip(g, cand)])
    g = g[free]
    if len(g) == 0:
        return {"void.volume_A3": 0.0, "void.n_points": 0.0}
    from .geometry import fibonacci_directions
    dirs = fibonacci_directions(64)
    keep = []
    for p in g:
        rel = xyz - p
        dist = np.linalg.norm(rel, axis=1)
        m = dist <= 10.0
        rr = rel[m]
        proj = rr @ dirs.T
        perp2 = (rr ** 2).sum(1)[:, None] - proj ** 2
        keep.append(((proj > 0) & (perp2 <= 1.7 ** 2)).any(0).mean() >= min_enclosure)
    n = int(np.sum(keep))
    return {"void.volume_A3": n * spacing ** 3, "void.n_points": float(n)}
