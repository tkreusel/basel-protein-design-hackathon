"""Apply the steering handoff's base-environment criteria to any structure.

Reimplemented from the DEFINITIONS in the handoff's kemp_base.py docstring and
constants (the handoff code itself is not executed: it is third-party data
that imports a framework we do not load). Validated first by reproducing the
handoff's own published arm means on its own design files.

Per structure (base = the catalytic Asp, ligand = 6NT):
  base_sasa_rel       SASA of base carboxylate O in context / isolated residue
  lig_sasa_frac       ligand SASA in complex / alone
  base_sc_polar_4p5   side-chain N/O atoms of OTHER residues within 4.5 A of a
                      base O          (atoms, not residues; includes Trp NE1)
  base_cation_6       Lys NZ / Arg NE,NH1,NH2 within 6 A of a base O
  nitro_hb_far/near   nitro O (O21,O11) with a protein H-bond donor within 3.3 A
                      whose distance to every base O is >= 6 A (far), or only
                      donors closer than 6 A to the base (near)
  lig_clash           protein heavy atom / ligand heavy atom pairs < 2.4 A

Usage as a library: measure(atom_array, lig_chain=None) -> dict
"""

from __future__ import annotations

import numpy as np
from biotite.structure import sasa
from scipy.spatial import cKDTree

ACID_O = ("OD1", "OD2")
NITRO = ("O21", "O11")
DONOR = {("ASN", "ND2"), ("GLN", "NE2"), ("HIS", "ND1"), ("HIS", "NE2"),
         ("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
         ("TRP", "NE1"), ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"),
         ("CYS", "SG")}
CATION = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")}
BB = ("N", "O", "OXT")


def _is_donor(res, atom):
    return (res, atom) in DONOR or (atom == "N" and res != "PRO")


def find_base(prot, lig):
    """The Asp whose carboxylate O is closest to the ligand's N3 (or, failing
    that, to any ligand N)."""
    n3 = lig.coord[lig.atom_name == "N3"]
    ref = n3 if len(n3) else lig.coord
    best, rid = 1e9, None
    for r in set(prot.res_id[prot.res_name == "ASP"]):
        m = (prot.res_id == r) & np.isin(prot.atom_name, ACID_O)
        if m.sum() < 2:
            continue
        d = np.linalg.norm(prot.coord[m][:, None, :] - ref[None, :, :], axis=-1).min()
        if d < best:
            best, rid = d, int(r)
    return rid, best


def measure(struct, prot_chain=None, lig_name="6NT") -> dict | None:
    s = struct[struct.element != "H"]
    lig = s[s.res_name == lig_name]
    if lig.array_length() == 0:
        return None
    if prot_chain is None:
        pc = [c for c in dict.fromkeys(s.chain_id[~s.hetero]) if c not in set(lig.chain_id)]
        prot_chain = pc[0] if pc else None
    prot = s[(~s.hetero) & (s.chain_id == prot_chain)]
    rid, dbase = find_base(prot, lig)
    if rid is None:
        return None
    out = {"base_resid": rid, "base_to_N3": float(dbase)}

    base = prot[prot.res_id == rid]
    bo = prot.coord[(prot.res_id == rid) & np.isin(prot.atom_name, ACID_O)]

    # --- SASA (ProtOr radii, as in the handoff) --------------------------
    full = prot + lig
    sa = sasa(full, vdw_radii="ProtOr", point_number=200)
    nprot = prot.array_length()
    bm = (prot.res_id == rid) & np.isin(prot.atom_name, ACID_O)
    iso = sasa(base, vdw_radii="ProtOr", point_number=200)
    out["base_sasa_rel"] = float(np.nansum(sa[:nprot][bm]) /
                                 max(np.nansum(iso[np.isin(base.atom_name, ACID_O)]), 1e-6))
    ls = np.nansum(sa[nprot:])
    out["lig_sasa_frac"] = float(ls / max(np.nansum(sasa(lig, vdw_radii="ProtOr",
                                                          point_number=200)), 1e-6))

    # --- neighbours of the base O ---------------------------------------
    other = prot[prot.res_id != rid]
    tree = cKDTree(bo)
    d_oth = tree.query(other.coord)[0]
    isbb = np.isin(other.atom_name, BB)
    polar = np.isin(other.element, ("N", "O"))
    out["base_sc_polar_4p5"] = float(((d_oth <= 4.5) & polar & ~isbb).sum())
    out["base_cation_6"] = float(sum(
        (rn, an) in CATION and d <= 6.0
        for rn, an, d in zip(other.res_name, other.atom_name, d_oth)))

    # --- nitro H-bonds ----------------------------------------------------
    donor_mask = np.array([_is_donor(r, a) for r, a in zip(prot.res_name, prot.atom_name)])
    dcoord, dres = prot.coord[donor_mask], prot.res_id[donor_mask]
    far = near = anyc = 0
    for o in lig.coord[np.isin(lig.atom_name, NITRO)]:
        d = np.linalg.norm(dcoord - o, axis=1)
        hit = d <= 3.3
        if not hit.any():
            continue
        anyc += 1
        dbo = np.array([np.linalg.norm(bo - c, axis=1).min() for c in dcoord[hit]])
        if (dbo >= 6.0).any():
            far += 1
        else:
            near += 1
    out.update(nitro_hb_far=float(far), nitro_hb_near_base=float(near),
               nitro_hb_any=float(anyc))

    pt = cKDTree(prot.coord)
    out["lig_clash"] = float(sum(len(x) for x in pt.query_ball_point(lig.coord, 2.4)))
    return out
