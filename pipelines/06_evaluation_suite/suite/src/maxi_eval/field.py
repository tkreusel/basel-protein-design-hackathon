"""Electric field projected on the reaction coordinates (Boxer/Fried-style), from tmol partial charges.

Point-charge Coulomb sum over the protein; the ligand (and optionally the catalytic base) is
excluded so the number is the field the *environment* exerts on the reacting bonds.

Sign convention - POSITIVE = favourable to the developing charge separation:
  field_CH : E . u(X -> H3)            positive pushes the transferring proton toward the base
  field_LG : -E . u(C7A -> O1/N1)      positive stabilises negative charge building on the leaving O
                                       (u points away from the ring, so the field must oppose it)
Evaluated at the bond midpoints. Units MV/cm (1 V/A = 100 MV/cm).

LIMITS (stated in the output, not hidden): vacuum point charges with Rosetta/MMFF partial charges,
no solvent screening and no polarisation -> use for RANKING and for ensemble VARIANCE, never as an
absolute number to compare with vibrational-Stark data. APBS/PBSA screened fields are the proper next
step and are not available in the tmol env. `rdd4` (eps = 4r) is an unscreened-vs-screened sensitivity check.
"""
from __future__ import annotations

import numpy as np

from . import geometry as G
from . import ligands as L

K_V_A = 14.3996  # V*Angstrom/e   (e^2/(4 pi eps0) in eV*A)


def field_vector(points: np.ndarray, xyz: np.ndarray, q: np.ndarray, mode: str = "vac") -> np.ndarray:
    """E [V/A] at each point from charges q at xyz. points [P,3]."""
    d = points[:, None, :] - xyz[None, :, :]
    r2 = (d ** 2).sum(-1)
    r2 = np.maximum(r2, 1e-6)
    r = np.sqrt(r2)
    if mode == "vac":
        w = q[None, :] / (r2 * r)           # q / r^3
    elif mode == "rdd4":
        w = q[None, :] / (2.0 * r2 * r2)    # eps=4r -> E = q r_hat / (2 r^3)  ==  q d / (2 r^4)
    else:
        raise ValueError(mode)
    return K_V_A * (w[:, :, None] * d).sum(1)


def reaction_axes(ly: dict, roles: dict) -> dict:
    X, O1, C7A = ly[roles["abstraction"]], ly[roles["oxyanion"]], ly[roles["oxyanion_partner"]]
    H = ly.get(roles.get("abstraction_H") or "", None)
    if H is None:
        H = G.ideal_h(X, ly[roles["ring_attached"]], ly[roles["ring_partner"]])
    return {"CH": ((X + H) / 2, G._u(H - X), +1.0), "LG": ((C7A + O1) / 2, G._u(O1 - C7A), -1.0)}


def project(pose_xyz, q, block_of_atom, ly, roles, lig_block: int, base_block: int | None, modes=("vac", "rdd4")) -> dict:
    """pose_xyz/q/block_of_atom from tmol_core.pose_charges(); ly = ligand name->xyz from the same pose."""
    axes = reaction_axes(ly, roles)
    out = {}
    sets = {"env": block_of_atom != lig_block}
    if base_block is not None:
        sets["env_nobase"] = (block_of_atom != lig_block) & (block_of_atom != base_block)
        sets["base_only"] = block_of_atom == base_block
    for sname, sel in sets.items():
        for mode in modes:
            for ax, (pt, u, sign) in axes.items():
                E = field_vector(pt[None], pose_xyz[sel], q[sel], mode)[0]
                out[f"field.{ax}.{sname}.{mode}"] = float(sign * (E @ u) * 100.0)  # MV/cm
    return out
