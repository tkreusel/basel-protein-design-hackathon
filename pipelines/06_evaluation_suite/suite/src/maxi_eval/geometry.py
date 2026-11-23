"""Active-site geometry of a Kemp eliminase complex. Pure numpy/biotite, no tmol, no GPU.

All atom roles come from ligands.ROLES (never hard-coded names), so the same code scores the
TSA (6NT), the substrate and the product, and any user-supplied TS template that follows the
canonical naming.

Metric families (all keys are prefixed `geo.`):

  base      d_O_H, d_O_X, ang_O_H_X, ang_C_O_H, tors_syn, oop_H, tors_C_O_H_X, tors_O_H_X_C3A, syn_flag
            (X = abstraction centre: N3 of 6NT / C3 of substrate+TS)
  oxyanion  n_donors, best-donor identity/distances/angle, donor_is_backbone
  wall      pi-stack distance + interplanar angle to the nearest aromatic residue
  burial    lig_rel_sasa, enclosure, n_contacts, n_shell, contact_min (same definitions as tobi_v1)
  water     proxy for lone-pair water sites on the catalytic carboxylate (see water_sites())
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from . import ligands as L

CARBOXYLATES = {"ASP": ("CG", "OD1", "OD2"), "GLU": ("CD", "OE1", "OE2")}
# protein atoms that can donate an H-bond: res -> atoms
DONORS = {
    "SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"], "LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"],
    "ASN": ["ND2"], "GLN": ["NE2"], "HIS": ["ND1", "NE2"], "TRP": ["NE1"], "CYS": ["SG"],
    "ASP": [], "GLU": [], "*": ["N"],  # backbone amide (not PRO)
}
AROMATIC = {"TRP": ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
            "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
            "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}


# ------------------------------------------------------------------ vector helpers
def _u(v):
    n = np.linalg.norm(v)
    return v / n if n > 0 else v


def dist(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)))


def angle(a, b, c) -> float:
    """Angle a-b-c in degrees."""
    v1, v2 = _u(np.asarray(a) - b), _u(np.asarray(c) - b)
    return float(np.degrees(np.arccos(np.clip(v1 @ v2, -1, 1))))


def dihedral(a, b, c, d) -> float:
    """Signed dihedral a-b-c-d in degrees (IUPAC)."""
    b0, b1, b2 = -(np.asarray(b) - a), np.asarray(c) - b, np.asarray(d) - c
    b1 = _u(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    x, y = v @ w, np.cross(b1, v) @ w
    return float(np.degrees(np.arctan2(y, x)))


def wrap180(x: float) -> float:
    return (x + 180.0) % 360.0 - 180.0


def circ_std_deg(angles_deg) -> float:
    """Circular standard deviation in degrees (sqrt(-2 ln R))."""
    a = np.radians(np.asarray(angles_deg, float))
    R = np.hypot(np.cos(a).mean(), np.sin(a).mean())
    return float(np.degrees(np.sqrt(-2 * np.log(max(R, 1e-12)))))


def plane_normal(pts) -> np.ndarray:
    pts = np.asarray(pts)
    _, _, vt = np.linalg.svd(pts - pts.mean(0))
    return vt[-1]


# ------------------------------------------------------------------ lookup tables
@dataclass
class Atoms:
    """name-indexed lookup into a biotite AtomArray (single chain / single residue set)."""
    arr: object
    _idx: dict

    @classmethod
    def of(cls, arr):
        idx = {}
        for i, (rid, an) in enumerate(zip(arr.res_id, arr.atom_name)):
            idx[(int(rid), str(an))] = i
        return cls(arr, idx)

    def xyz(self, rid, name):
        i = self._idx.get((int(rid), name))
        return None if i is None else self.arr.coord[i]

    def resname(self, rid):
        for (r, _), i in self._idx.items():
            if r == rid:
                return str(self.arr.res_name[i])
        return None


def lig_xyz(rep) -> dict[str, np.ndarray]:
    return {str(n): c for n, c in zip(rep.ligand.atom_name, rep.ligand.coord)}


def ideal_h(x_pos, n1, n2, length=1.0):
    """Hydrogen on a planar (sp2) centre X bonded to n1, n2: along the external bisector."""
    b = _u(_u(x_pos - n1) + _u(x_pos - n2))
    return x_pos + length * b


def abstraction_h(rep, roles, ly):
    """(H position, 'model'|'ideal'). Uses the ligand's own H3 if present, else builds it."""
    h = roles.get("abstraction_H")
    if h and h in ly:
        return ly[h], "model"
    X = ly[roles["abstraction"]]
    return ideal_h(X, ly[roles["ring_attached"]], ly[roles["ring_partner"]]), "ideal"


def find_donors(rep, roles, cutoff: float = 3.6) -> list[dict]:
    """Protein heavy atoms that can donate an H-bond and lie within `cutoff` A of the leaving-group atom, from the CURRENT coordinates."""
    ly = lig_xyz(rep)
    O1 = ly[roles["oxyanion"]]
    arr = rep.protein
    P = Atoms.of(arr)
    donors = []
    for (rid, an), i in P._idx.items():
        rn = str(arr.res_name[i])
        if an in DONORS.get(rn, []) or (an == "N" and rn != "PRO"):
            d = dist(arr.coord[i], O1)
            if d <= cutoff:
                donors.append(dict(resid=int(rid), resname=rn, atom=an, d=d, backbone=an == "N"))
    donors.sort(key=lambda x: x["d"])
    return donors


# ------------------------------------------------------------------ site detection
def detect_site(rep, overrides: dict | None = None) -> dict:
    """Identify base / donors / wall / pocket from the TSA-bound geometry (or override)."""
    ov = overrides or {}
    roles = L.ROLES[rep.lig_resname]
    ly = lig_xyz(rep)
    P = Atoms.of(rep.protein)
    arr = rep.protein
    X = ly[roles["abstraction"]]
    O1 = ly[roles["oxyanion"]]
    site: dict = {"roles": roles}

    # --- base: carboxylate O closest to the abstraction centre
    cands = []
    for (rid, an), i in P._idx.items():
        rn = str(arr.res_name[i])
        if rn in CARBOXYLATES and an in CARBOXYLATES[rn][1:]:
            if ov.get("base_resid") and rid != ov["base_resid"]:
                continue
            cands.append((dist(arr.coord[i], X), rid, rn, an))
    if cands:
        d, rid, rn, an = min(cands)
        c_name, o1n, o2n = CARBOXYLATES[rn]
        other = o2n if an == o1n else o1n
        site["base"] = dict(resid=int(rid), resname=rn, C=c_name, O=an, O_other=other, d_to_X=d)
    else:
        site["base"] = None

    site["donors"] = find_donors(rep, roles)

    # --- wall: aromatic ring whose centroid is closest to the ligand benzo centroid
    bz = np.mean([ly[n] for n in roles["benzo"]], axis=0)
    best = None
    for rid in sorted({r for r, _ in P._idx}):
        rn = P.resname(rid)
        if rn in AROMATIC:
            pts = [P.xyz(rid, a) for a in AROMATIC[rn] if P.xyz(rid, a) is not None]
            if len(pts) >= 5:
                cen = np.mean(pts, axis=0)
                d = dist(cen, bz)
                if d < 7.0 and (best is None or d < best["d"]):
                    best = dict(resid=int(rid), resname=rn, ring=AROMATIC[rn], d=d)
    if ov.get("wall_resid"):
        rid = int(ov["wall_resid"])
        rn = P.resname(rid)
        best = dict(resid=rid, resname=rn, ring=AROMATIC.get(rn, []), d=np.nan)
    site["wall"] = best

    # --- pocket residues: any heavy atom within 6 A of the ligand
    lc = rep.ligand.coord
    dd = np.linalg.norm(arr.coord[:, None] - lc[None], axis=-1).min(1)
    site["pocket_resids"] = sorted({int(r) for r in arr.res_id[dd < 6.0]})
    site["pocket_resids_8A"] = sorted({int(r) for r in arr.res_id[dd < 8.0]})
    return site


# ------------------------------------------------------------------ geometry metrics
def base_metrics(rep, site) -> dict:
    out: dict = {}
    b = site.get("base")
    roles = site["roles"]
    ly = lig_xyz(rep)
    P = Atoms.of(rep.protein)
    if b is None:
        return {"geo.base_found": 0.0}
    X = ly[roles["abstraction"]]
    H, hsrc = abstraction_h(rep, roles, ly)
    rid = b["resid"]
    C, O, Oo = P.xyz(rid, b["C"]), P.xyz(rid, b["O"]), P.xyz(rid, b["O_other"])
    n = _u(np.cross(Oo - C, O - C))          # carboxylate plane normal
    HO = H - O
    out.update({
        "geo.base_found": 1.0,
        "geo.base_resid": float(rid),
        "geo.d_O_H": dist(O, H),
        "geo.d_O_X": dist(O, X),
        "geo.d_Oother_H": dist(Oo, H),
        "geo.ang_O_H_X": angle(O, H, X),                  # ~180 = linear proton transfer
        "geo.ang_C_O_H": angle(C, O, H),                  # in-plane lone pair ~ 110-130
        "geo.tors_syn": dihedral(Oo, C, O, H),            # ~0 syn lone pair, ~180 anti
        "geo.oop_H": float(np.degrees(np.arcsin(np.clip(_u(HO) @ n, -1, 1)))),  # 0 = in carboxylate plane
        "geo.tors_C_O_H_X": dihedral(C, O, H, X),         # rotation about the O...H axis
        "geo.tors_O_H_X_C3A": dihedral(O, H, X, ly[roles["ring_attached"]]),  # approach vs ring plane
        "geo.h_source_model": 1.0 if hsrc == "model" else 0.0,
    })
    out["geo.syn_flag"] = 1.0 if abs(out["geo.tors_syn"]) < 60 else 0.0
    return out


def oxyanion_metrics(rep, site, hbond_max=3.6) -> dict:
    """Donors to the leaving-group atom. Uses real H atoms when present, else heavy-atom only."""
    roles = site["roles"]
    ly = lig_xyz(rep)
    A = ly[roles["oxyanion"]]
    arr = rep.protein
    P = Atoms.of(arr)
    donors = find_donors(rep, roles)   # recomputed for THIS structure (relaxed pose / PLACER model), not the input
    site = {**site, "donors": donors}
    out = {"geo.n_donors_3p5": float(sum(1 for d in donors if d["d"] <= 3.5))}
    if not donors:
        out.update({"geo.donor_best_d": float("nan"), "geo.donor_best_is_bb": float("nan"),
                    "geo.donor_best_angle": float("nan")})
        return out
    # side-chain donors rank above backbone ones for the "designed" oxyanion interaction
    sc = [d for d in site["donors"] if not d["backbone"]]
    out["geo.n_sidechain_donors_3p5"] = float(sum(1 for d in sc if d["d"] <= 3.5))
    best = (sc or site["donors"])[0]
    out["geo.donor_best_d"] = best["d"]
    out["geo.donor_best_is_bb"] = 1.0 if best["backbone"] else 0.0
    out["geo.donor_best_resid"] = float(best["resid"])
    # D-H...A angle if H present on that donor (tmol-protonated models)
    Hs = [c for (r, an), i in P._idx.items() if r == best["resid"] and an.startswith("H")
          for c in [arr.coord[i]] if dist(c, P.xyz(best["resid"], best["atom"])) < 1.2]
    if Hs:
        h = min(Hs, key=lambda c: dist(c, A))
        out["geo.donor_best_H_d"] = dist(h, A)
        out["geo.donor_best_angle"] = angle(P.xyz(best["resid"], best["atom"]), h, A)
    else:
        out["geo.donor_best_H_d"] = float("nan")
        out["geo.donor_best_angle"] = float("nan")
    # nearest donor regardless of class
    out["geo.donor_any_d"] = site["donors"][0]["d"]
    return out


def wall_metrics(rep, site) -> dict:
    w = site.get("wall")
    if not w:
        return {"geo.wall_found": 0.0}
    roles = site["roles"]
    ly = lig_xyz(rep)
    P = Atoms.of(rep.protein)
    pts = np.array([P.xyz(w["resid"], a) for a in w["ring"] if P.xyz(w["resid"], a) is not None])
    bz = np.array([ly[n] for n in roles["benzo"]])
    if pts.ndim != 2 or len(pts) < 3 or bz.ndim != 2 or len(bz) < 3:   # ring atoms missing (e.g. PLACER pocket crop) -> no plane
        return {"geo.wall_found": 0.0}
    n1, n2 = plane_normal(pts), plane_normal(bz)
    ang = float(np.degrees(np.arccos(abs(np.clip(n1 @ n2, -1, 1)))))
    return {"geo.wall_found": 1.0, "geo.wall_resid": float(w["resid"]),
            "geo.wall_d_centroid": dist(pts.mean(0), bz.mean(0)), "geo.wall_plane_angle": ang}


# ------------------------------------------------------------------ burial (tobi_v1 definitions)
def fibonacci_directions(n: int = 256) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def enclosure_fraction(lig_xyz_, prot_xyz, max_dist=12.0, tol=2.8, n_dir=256) -> float:
    centre = lig_xyz_.mean(0)
    rel = prot_xyz - centre
    near = rel[np.linalg.norm(rel, axis=1) <= max_dist]
    if len(near) == 0:
        return 0.0
    dirs = fibonacci_directions(n_dir)
    proj = near @ dirs.T
    perp2 = (near ** 2).sum(1)[:, None] - proj ** 2
    return float(((proj > 0) & (perp2 <= tol ** 2)).any(0).mean())


def vdw_radii(arr) -> np.ndarray:
    """ProtOr radii where biotite knows the residue (protein, 6NT), element radii otherwise (H5J/6VP: element radii)."""
    from biotite.structure.info import vdw_radius_protor, vdw_radius_single

    out = np.empty(len(arr))
    for i, (rn, an, el) in enumerate(zip(arr.res_name, arr.atom_name, arr.element)):
        r = None
        try:
            r = vdw_radius_protor(str(rn), str(an))
        except Exception:  # noqa: BLE001
            r = None
        out[i] = r if r is not None else (vdw_radius_single(str(el)) or 1.8)
    return out


def burial_metrics(rep) -> dict:
    from biotite.structure import sasa

    prot = rep.protein[rep.protein.element != "H"]
    lig = rep.ligand[rep.ligand.element != "H"]
    both = prot + lig
    s_c = sasa(both, vdw_radii=vdw_radii(both), ignore_ions=False)
    lig_c = float(np.nansum(s_c[len(prot):]))
    lig_f = float(np.nansum(sasa(lig, vdw_radii=vdw_radii(lig), ignore_ions=False)))
    d = np.linalg.norm(prot.coord[:, None] - lig.coord[None], axis=-1)
    return {"geo.lig_rel_sasa": lig_c / lig_f if lig_f > 0 else float("nan"),
            "geo.enclosure": enclosure_fraction(lig.coord, prot.coord),
            "geo.n_contacts": float((d.min(1) < 5.0).sum()),
            "geo.n_shell": float((d.min(1) < 8.0).sum()),
            "geo.contact_min": float(d.min())}


def base_sasa_metrics(rep, site) -> dict:
    """SASA of the catalytic carboxylate (C,O,O') in the complex vs the isolated residue."""
    from biotite.structure import sasa

    b = site.get("base")
    if b is None:
        return {}
    prot = rep.protein[rep.protein.element != "H"]
    both = prot + rep.ligand[rep.ligand.element != "H"]
    s = sasa(both, vdw_radii=vdw_radii(both), ignore_ions=False)
    sel = (both.res_id == b["resid"]) & np.isin(both.atom_name, [b["C"], b["O"], b["O_other"]]) & (both.chain_id == "A")
    res = prot[prot.res_id == b["resid"]]
    s0 = sasa(res, vdw_radii=vdw_radii(res), ignore_ions=False)
    sel0 = np.isin(res.atom_name, [b["C"], b["O"], b["O_other"]])
    a, a0 = float(np.nansum(s[sel])), float(np.nansum(s0[sel0]))
    return {"geo.base_sasa": a, "geo.base_sasa_free_frac": a / a0 if a0 > 0 else float("nan")}


# ------------------------------------------------------------------ lone-pair water-site proxy
def water_sites(rep, site, d_site=2.8, clash=2.6, partner=1.3) -> dict:
    """Proxy for waters H-bonded to the catalytic carboxylate (NOT a water placement method).

    For each carboxylate O two in-plane lone-pair directions (C-O...X = 120 deg, syn and anti)
    define ideal H-bond water positions X at 2.8 A. A site is
      occupied : a protein/ligand heavy atom already sits at the site (<= `partner` A) -> satisfied
      open     : nothing within `clash` A -> a water can sit there and compete for the lone pair
      blocked  : otherwise.
    Validated against crystal waters in tests/validate_water_proxy.py. Real placement
    (3D-RISM/GIST/Dowser) is not available in tmol -- see PIPELINE.md.
    """
    b = site.get("base")
    if b is None:
        return {}
    P = Atoms.of(rep.protein)
    C = P.xyz(b["resid"], b["C"])
    pos = {"O": P.xyz(b["resid"], b["O"]), "Oo": P.xyz(b["resid"], b["O_other"])}
    allxyz = np.vstack([rep.protein.coord, rep.ligand.coord])
    skip = (rep.protein.res_id == b["resid"])
    out = {"geo.wsite_open_total": 0.0, "geo.wsite_occupied_total": 0.0}
    for key, O in pos.items():
        Oo = pos["Oo"] if key == "O" else pos["O"]
        u = _u(O - C)
        p = _u((Oo - C) - ((Oo - C) @ u) * u)      # in-plane, toward the other O
        n_open = n_occ = 0
        for sgn in (+1, -1):
            d = 0.5 * u + sgn * (math.sqrt(3) / 2) * p
            X = O + d_site * d
            dd = np.linalg.norm(allxyz - X, axis=1)
            dd = np.where(np.concatenate([skip, np.zeros(len(rep.ligand), bool)]), 99.0, dd)
            if dd.min() <= partner:
                n_occ += 1
            elif dd.min() > clash:
                n_open += 1
        out[f"geo.wsite_open_{'prox' if key == 'O' else 'dist'}"] = float(n_open)
        out["geo.wsite_open_total"] += n_open
        out["geo.wsite_occupied_total"] += n_occ
    # crystal waters (validation only)
    if len(rep.waters):
        for key, O in pos.items():
            out[f"geo.xtal_waters_{'prox' if key == 'O' else 'dist'}"] = float(
                (np.linalg.norm(rep.waters - O, axis=1) <= 3.3).sum())
    return out


# ------------------------------------------------------------------ one call for everything
def all_geometry(rep, site=None, overrides=None) -> tuple[dict, dict]:
    site = site or detect_site(rep, overrides)
    m: dict = {}
    # base/proton-transfer geometry only makes sense for ligands with an abstractable proton (6NT: N3-H3, H5J: C3-H2);
    # the product 6VP has none, so it gets oxyanion / wall / burial / water metrics only.
    if L.ROLES[rep.lig_resname]["abstraction_H"] is not None:
        m.update(base_metrics(rep, site))
    m.update(oxyanion_metrics(rep, site))
    m.update(wall_metrics(rep, site))
    m.update(burial_metrics(rep))
    m.update(base_sasa_metrics(rep, site))
    m.update(water_sites(rep, site))
    return m, site


def chi_angles(rep, resid: int) -> dict[str, float]:
    """chi1..chi4 (degrees) of one residue (from heavy atoms)."""
    CHI = {"ASP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
           "GLU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
           "GLN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
           "ASN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
           "TRP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
           "TYR": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
           "PHE": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
           "HIS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
           "LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")],
           "SER": [("N", "CA", "CB", "OG")], "THR": [("N", "CA", "CB", "OG1")],
           "ARG": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")]}
    P = Atoms.of(rep.protein)
    rn = P.resname(resid)
    out = {}
    for k, atoms in enumerate(CHI.get(rn, []), start=1):
        xs = [P.xyz(resid, a) for a in atoms]
        if all(x is not None for x in xs):
            out[f"chi{k}"] = dihedral(*xs)
    return out
