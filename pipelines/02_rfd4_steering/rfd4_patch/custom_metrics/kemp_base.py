"""Base environment + nitro H-bond + clash descriptors for Kemp eliminase designs and crystal calibration structures.

Generic (no theozyme reference needed): the catalytic base is the ASP/GLU whose carboxylate O is closest to a ligand ring N (N1/N2/N3 of 6NT).
Outputs (floats; prefixed by the config key):
  base_lig_dist        closest carboxylate-O .. ligand-N distance (A)                          [sanity: ~2.7-3.2 in crystals]
  base_sasa_rel        SASA of the base O atoms in context / isolated residue (0 = buried)
  base_apolar_frac     C/S fraction of protein heavy atoms (other residues) within `env_r` of the base O atoms
  base_sc_polar_4p5    side-chain polar atoms (N/O) of OTHER residues within `near_r` of a base O  (the thing to minimise)
  base_bb_polar_4p5    backbone N/O of other residues within `near_r` of a base O
  base_cation_6        Lys NZ / Arg NE,NH1,NH2 within `cat_r` of a base O (salt-bridge / field)
  nitro_hb_far         nitro O atoms (0-2) with a protein H-bond DONOR within `hb_cut` whose distance to every base O is >= `far_r`
  nitro_hb_any         nitro O atoms with any protein donor within `hb_cut` (regardless of where the donor is)
  nitro_hb_near_base   nitro O atoms whose only donors are closer than `far_r` to the base
  nitro_sasa_rel       mean relative SASA of the two nitro O atoms
  lig_clash_atoms      protein heavy-atom/ligand heavy-atom pairs closer than `lig_clash` (A)
  pp_clash_frac        protein-protein heavy-atom pairs (|dres|>1) closer than `pp_clash`, per protein atom
"""
from __future__ import annotations

import numpy as np
from biotite.structure import AtomArray, sasa
from scipy.spatial import cKDTree

from rfproteina.metrics.metric import Metric

_ACID_O = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
_NITRO = ("O21", "O11")
_DONOR = {("*", "N"), ("ASN", "ND2"), ("GLN", "NE2"), ("HIS", "ND1"), ("HIS", "NE2"), ("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"),
          ("ARG", "NH2"), ("TRP", "NE1"), ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("CYS", "SG")}
_CATION = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")}
_BB = ("N", "O", "OXT")


def _is_donor(res, atom):
    return (res, atom) in _DONOR or atom == "N" and res != "PRO"


class KempBaseMetric(Metric):
    def __init__(self, env_r: float = 6.0, near_r: float = 4.5, cat_r: float = 6.0, hb_cut: float = 3.3, far_r: float = 6.0,
                 lig_clash: float = 2.4, pp_clash: float = 2.2, **kwargs):
        super().__init__(**kwargs)
        self.p = dict(env_r=env_r, near_r=near_r, cat_r=cat_r, hb_cut=hb_cut, far_r=far_r, lig_clash=lig_clash, pp_clash=pp_clash)

    @property
    def kwargs_to_compute_args(self) -> dict:
        return {"atom_array": ("predicted_atom_array",)}

    def compute(self, atom_array: AtomArray) -> dict[str, float]:
        P = self.p
        a = atom_array[(atom_array.element != "H") & ~np.isnan(atom_array.coord).any(1)]
        lig, prot = a[a.hetero], a[~a.hetero]
        keys = ("base_lig_dist", "base_sasa_rel", "base_apolar_frac", "base_sc_polar_4p5", "base_bb_polar_4p5", "base_cation_6", "nitro_hb_far",
                "nitro_hb_any", "nitro_hb_near_base", "nitro_sasa_rel", "lig_clash_atoms", "pp_clash_frac")
        out = {k: float("nan") for k in keys}
        if len(lig) == 0 or len(prot) < 10:
            return out
        pt = cKDTree(prot.coord)
        out["lig_clash_atoms"] = float(sum(len(x) for x in pt.query_ball_point(lig.coord, P["lig_clash"])))
        pr = pt.query_pairs(P["pp_clash"], output_type="ndarray")
        out["pp_clash_frac"] = float((np.abs(prot.res_id[pr[:, 0]] - prot.res_id[pr[:, 1]]) > 1).sum() / len(prot)) if len(pr) else 0.0
        ring_n = lig[np.isin(lig.atom_name, ("N1", "N2", "N3"))]
        best = (np.inf, None)
        for rid in np.unique(prot.res_id[np.isin(prot.res_name, list(_ACID_O))]):
            r = prot[(prot.res_id == rid) & np.isin(prot.res_name, list(_ACID_O))]
            o = r[np.isin(r.atom_name, _ACID_O[r.res_name[0]])]
            if len(o) == 0 or len(ring_n) == 0:
                continue
            d = np.linalg.norm(o.coord[:, None] - ring_n.coord[None], axis=-1).min()
            if d < best[0]:
                best = (d, rid)
        if best[1] is None:
            return out
        rb = best[1]; sel_b = prot.res_id == rb
        base = prot[sel_b & np.isin(prot.res_name, list(_ACID_O))]; bo = base[np.isin(base.atom_name, _ACID_O[base.res_name[0]])]
        out["base_lig_dist"] = float(best[0])
        sa = sasa(a, vdw_radii="ProtOr", point_number=200)
        iso = sasa(base, vdw_radii="ProtOr", point_number=200)
        ctx = sa[(~a.hetero) & (a.res_id == rb) & np.isin(a.res_name, list(_ACID_O)) & np.isin(a.atom_name, _ACID_O[base.res_name[0]])]
        out["base_sasa_rel"] = float(np.nansum(ctx) / max(np.nansum(iso[np.isin(base.atom_name, _ACID_O[base.res_name[0]])]), 1e-6))
        other = prot[~sel_b]; ot = cKDTree(other.coord)
        env = other[np.unique(np.concatenate(ot.query_ball_point(bo.coord, P["env_r"]))).astype(int)]
        pol = np.isin(env.element, ("N", "O")).sum(); apo = (np.isin(env.element, ("C", "S")) & ~np.isin(env.atom_name, ("C", "CA"))).sum()
        out["base_apolar_frac"] = float(apo / max(pol + apo, 1))
        near = other[np.unique(np.concatenate(ot.query_ball_point(bo.coord, P["near_r"]))).astype(int)]
        npol = near[np.isin(near.element, ("N", "O"))]; isbb = np.isin(npol.atom_name, _BB)
        out["base_sc_polar_4p5"] = float((~isbb).sum()); out["base_bb_polar_4p5"] = float(isbb.sum())
        cat = other[np.unique(np.concatenate(ot.query_ball_point(bo.coord, P["cat_r"]))).astype(int)]
        out["base_cation_6"] = float(sum((r, n) in _CATION for r, n in zip(cat.res_name, cat.atom_name)))
        nit = lig[np.isin(lig.atom_name, _NITRO)]; far = anyc = nearb = 0
        for c in nit.coord:
            hits = [i for i in pt.query_ball_point(c, P["hb_cut"]) if _is_donor(prot.res_name[i], prot.atom_name[i]) and not sel_b[i]]
            if hits:
                anyc += 1
                dbase = [np.linalg.norm(bo.coord - prot.coord[i], axis=1).min() for i in hits]
                if max(dbase) >= P["far_r"]: far += 1
                else: nearb += 1
        out.update(nitro_hb_far=float(far), nitro_hb_any=float(anyc), nitro_hb_near_base=float(nearb))
        ns = sa[a.hetero][np.isin(lig.atom_name, _NITRO)]; iso_l = sasa(lig, vdw_radii="ProtOr", point_number=200)[np.isin(lig.atom_name, _NITRO)]
        out["nitro_sasa_rel"] = float(np.mean(ns / np.maximum(iso_l, 1e-6)))
        return out
