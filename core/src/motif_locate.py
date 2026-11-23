"""Locate unindexed motif residues in an RFD4 output without its provenance.

RFD4 recentres its output, so conditioned atoms do not keep their input
coordinates, and the provenance columns (`annotation_src-component_1_atom`) are
sometimes missing (6 of 32 tip-atom designs). The ligand is always conditioned,
so it anchors the frame:

  1. Kabsch-superpose the input ligand onto the output ligand (same atoms).
  2. Carry each input motif residue's tip atoms through that transform.
  3. For each role pick the output residue of the expected identity whose tip
     atoms are closest to those positions; report the tip RMSD.

A role is only accepted when the tip RMSD is below `tol` (default 0.6 A), so a
failure to honour the conditioning is reported, not silently mapped to the
nearest wrong residue.
"""

from __future__ import annotations

import numpy as np


def kabsch_transform(P: np.ndarray, Q: np.ndarray):
    pc, qc = P.mean(0), Q.mean(0)
    v, _, wt = np.linalg.svd((P - pc).T @ (Q - qc))
    if np.linalg.det(v @ wt) < 0:
        v[:, -1] *= -1
    R = v @ wt
    return R, pc, qc


def locate(design, motif_in, roles: dict, lig_name="6NT", tol=0.6) -> dict:
    """roles: {role: dict(resid=<input res_id>, resname=<expected output name>,
                           tips=[atom names])}
    Returns {role: dict(resid=<output res_id>, rmsd=float)} or raises ValueError.
    """
    dl = design[(design.res_name == lig_name)]
    il = motif_in[(motif_in.res_name == lig_name)]
    names = [n for n in il.atom_name if n in set(dl.atom_name)]
    if len(names) < 6:
        raise ValueError("ligand atoms do not match between input and output")
    P = np.array([il.coord[il.atom_name == n][0] for n in names])
    Q = np.array([dl.coord[dl.atom_name == n][0] for n in names])
    R, pc, qc = kabsch_transform(P, Q)
    lig_rmsd = float(np.sqrt((((P - pc) @ R + qc - Q) ** 2).sum(1).mean()))
    if lig_rmsd > 0.3:
        raise ValueError(f"ligand frame fit poor (rmsd {lig_rmsd:.2f} A)")

    out = {}
    prot = design[~design.hetero]
    for role, spec in roles.items():
        tips = spec["tips"]
        ref = []
        for n in tips:
            m = (motif_in.res_id == spec["resid"]) & (motif_in.atom_name == n) & (~motif_in.hetero)
            if not m.any():
                raise ValueError(f"input motif lacks {role} atom {n}")
            ref.append((motif_in.coord[m][0] - pc) @ R + qc)
        ref = np.array(ref)
        best = (1e9, None)
        for rid in sorted(set(prot.res_id[prot.res_name == spec["resname"]])):
            cur = []
            for n in tips:
                m = (prot.res_id == rid) & (prot.atom_name == n)
                if not m.any():
                    break
                cur.append(prot.coord[m][0])
            if len(cur) != len(tips):
                continue
            r = float(np.sqrt(((np.array(cur) - ref) ** 2).sum(1).mean()))
            if r < best[0]:
                best = (r, int(rid))
        if best[1] is None or best[0] > tol:
            raise ValueError(f"{role}: no {spec['resname']} within {tol} A of its "
                             f"conditioned tips (best {best[0]:.2f})")
        out[role] = dict(resid=best[1], rmsd=best[0])
    return out
