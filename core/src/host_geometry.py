"""Host-style scoring components: approach geometry and base burial.

The hosts' slide (2026-11-21) scores designs on four things:
  iPDE         Boltz-2 protein-ligand PDE (we use Boltz's own complex_ipde)
  Geometry     distance / angle / dihedral measurements against the ideal
               approach geometry for the catalytic base and the H-bond donor
  Burial       number of CA atoms within 12 A of the catalytic base
  Consistency  reproducibility of those metrics across Boltz-2 seeds

This module implements Geometry and Burial. Their exact formulas are not public,
so the "ideal" values are MEASURED from the HG3-series crystal structures with
6NT bound (see scoring/calibrate_host_reference.py), not assumed.

Ligand atoms (6NT, from configs/atom_map.yaml):
  N3  abstraction centre (carries the acidic H)  bonded C3A, N2
  N1  oxyanion centre                            bonded C7A, N2
Base = catalytic Asp (CG, OD1, OD2). Donor = Gln (CD, OE1, NE2).
The N-H axis at N3 is the EXTERIOR bisector of C3A-N3-N2 (heavy-atom proxy for
the H3 direction, since predicted structures carry no hydrogens).
"""

from __future__ import annotations

import numpy as np

RING = ["C3A", "C4", "C5", "C6", "C7", "C7A", "N1", "N2", "N3"]

BASE_METRICS = ["b_dist", "b_ang_O", "b_ang_NH", "b_oop", "b_dih"]
DONOR_METRICS = ["d_dist", "d_ang_N", "d_ang_ax", "d_oop", "d_dih"]
ANGULAR = {"b_ang_O", "b_ang_NH", "b_oop", "b_dih", "d_ang_N", "d_ang_ax", "d_oop", "d_dih"}
DIHEDRAL = {"b_dih", "d_dih"}


def _u(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def angle(a, b, c):
    """Angle at b, degrees."""
    return float(np.degrees(np.arccos(np.clip(np.dot(_u(a - b), _u(c - b)), -1, 1))))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = _u(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def _plane_normal(pts):
    c = pts.mean(0)
    return _u(np.linalg.svd(pts - c)[2][-1])


def _coord(s, chain, rid, name, hetero=False):
    m = (s.chain_id == chain) & (s.res_id == rid) & (s.atom_name == name) & (s.hetero == hetero)
    return s.coord[m][0] if m.any() else None


def _lig(s, lig_chain, lig_resid):
    m = (s.chain_id == lig_chain) & (s.res_id == lig_resid) & (s.res_name == "6NT")
    return {n: c for n, c in zip(s.atom_name[m], s.coord[m])}


def geometry(s, *, chain, base_resid, donor_resid, lig_chain, lig_resid) -> dict:
    """The ten approach-geometry measurements for one ligand copy."""
    L = _lig(s, lig_chain, lig_resid)
    need = ["N3", "C3A", "N2", "C4", "N1", "C7A", "C7"]
    if any(n not in L for n in need):
        raise ValueError("ligand atoms missing")
    normal = _plane_normal(np.array([L[n] for n in RING if n in L]))

    def oop(v):  # angle between vector and the ring plane, 0..90
        return float(np.degrees(np.arcsin(abs(np.dot(_u(v), normal)))))

    out = {}
    # --- base: proximal carboxylate O ------------------------------------
    cg = _coord(s, chain, base_resid, "CG")
    os_ = [_coord(s, chain, base_resid, n) for n in ("OD1", "OD2")]
    if cg is None or any(o is None for o in os_):
        raise ValueError("base atoms missing")
    ob = min(os_, key=lambda o: np.linalg.norm(o - L["N3"]))
    nh_axis = _u(-(_u(L["C3A"] - L["N3"]) + _u(L["N2"] - L["N3"])))
    v = ob - L["N3"]
    out["b_dist"] = float(np.linalg.norm(v))
    out["b_ang_O"] = angle(cg, ob, L["N3"])
    out["b_ang_NH"] = float(np.degrees(np.arccos(np.clip(np.dot(_u(v), nh_axis), -1, 1))))
    out["b_oop"] = oop(v)
    out["b_dih"] = dihedral(ob, L["N3"], L["C3A"], L["C4"])
    # --- donor: Gln NE2 -> N1 ---------------------------------------------
    ne2, cd = (_coord(s, chain, donor_resid, n) for n in ("NE2", "CD"))
    if ne2 is None or cd is None:
        raise ValueError("donor atoms missing")
    ax1 = _u(-(_u(L["C7A"] - L["N1"]) + _u(L["N2"] - L["N1"])))
    v = ne2 - L["N1"]
    out["d_dist"] = float(np.linalg.norm(v))
    out["d_ang_N"] = angle(cd, ne2, L["N1"])
    out["d_ang_ax"] = float(np.degrees(np.arccos(np.clip(np.dot(_u(v), ax1), -1, 1))))
    out["d_oop"] = oop(v)
    out["d_dih"] = dihedral(ne2, L["N1"], L["C7A"], L["C7"])
    return out


def burial_count(s, *, chain, base_resid, radius=12.0) -> dict:
    """CA atoms within `radius` of the catalytic base (own CA excluded).

    Primary: measured from the carboxylate carbon CG. Also reported from the
    midpoint of OD1/OD2, since the slide does not say which base atom is used.
    """
    ca = s[(s.chain_id == chain) & (s.atom_name == "CA") & (~s.hetero)]
    ca = ca[ca.res_id != base_resid]
    cg = _coord(s, chain, base_resid, "CG")
    o1, o2 = _coord(s, chain, base_resid, "OD1"), _coord(s, chain, base_resid, "OD2")
    if cg is None or o1 is None or o2 is None:
        raise ValueError("base atoms missing")
    mid = (o1 + o2) / 2
    return {"burial_CG": int((np.linalg.norm(ca.coord - cg, axis=1) <= radius).sum()),
            "burial_O": int((np.linalg.norm(ca.coord - mid, axis=1) <= radius).sum())}


def circ_diff(a, b):
    d = (a - b + 180.0) % 360.0 - 180.0
    return d


def circ_mean(deg):
    r = np.radians(deg)
    return float(np.degrees(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def circ_sd(deg):
    r = np.radians(deg)
    R = np.hypot(np.sin(r).mean(), np.cos(r).mean())
    return float(np.degrees(np.sqrt(-2 * np.log(max(R, 1e-9)))))
