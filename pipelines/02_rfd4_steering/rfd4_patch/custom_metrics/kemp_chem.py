"""Kemp-eliminase chemistry descriptors around the fixed theozyme (base, H-bond partner, TS-analogue ligand).

The motif residues are unindexed, so they are located by superposing the reference theozyme on the design's ligand and picking the
same-type residue that best matches the reference side-chain atoms. Everything is geometric (heavy atoms, no hydrogens): it screens,
it does not replace TMol/QM/MD.
Outputs (floats; prefixed by the config key):
  base_rmsd / hb_rmsd      match of the located base / partner residue to the reference (A); large = motif not found
  base_sasa_rel            SASA of the base's catalytic O atoms / their SASA in the isolated residue (0 = fully buried)
  base_apolar_frac         fraction of C/S among protein heavy atoms (other residues) within `env_r` of the base O atoms (hydrophobic environment)
  base_n_polar_env         protein N/O atoms (other residues) within `env_r` of the base O atoms (competing polarity / water-like sites)
  base_hbond_partners      protein N/O within `hb_cut` of the base O atoms, excluding the base itself (each lowers the base's basicity)
  hb_partner_dist          distance of the partner atom to the nearest ligand N/O (A)
  lig_hbonds               ligand N/O atoms with >=1 protein N/O within `hb_cut` (partner and base included)
  nitro_hbonds             the same for the nitro O atoms only (TS charge delocalisation / leaving-group side)
  nitro_unsat_buried       nitro O atoms with no partner within `hb_cut` AND SASA < `buried_sasa` A^2 (desolvated without compensation)
  motif_sasa_rel           mean relative SASA of base O, partner atom and all ligand atoms
"""
from __future__ import annotations

import numpy as np
from biotite.structure import AtomArray, sasa, superimpose
from biotite.structure.io import pdb
from scipy.spatial import cKDTree

from rfproteina.metrics.metric import Metric

_NITRO = ("O21", "O11")


class KempChemMetric(Metric):
    def __init__(self, ref_pdb: str, base_res: str = "GLU", base_match: tuple = ("CG", "CD", "OE1", "OE2"),
                 base_o: tuple = ("OE1", "OE2"), hb_res: str = "SER", hb_match: tuple = ("CA", "CB", "OG"), hb_atom: str = "OG",
                 env_r: float = 6.0, hb_cut: float = 3.3, buried_sasa: float = 5.0, **kwargs):
        super().__init__(**kwargs)
        ref = pdb.PDBFile.read(ref_pdb).get_structure(model=1)
        self.ref = ref[ref.element != "H"]
        self.p = dict(base_res=base_res, base_match=tuple(base_match), base_o=tuple(base_o), hb_res=hb_res, hb_match=tuple(hb_match), hb_atom=hb_atom,
                      env_r=env_r, hb_cut=hb_cut, buried_sasa=buried_sasa)

    @property
    def kwargs_to_compute_args(self) -> dict:
        return {"atom_array": ("predicted_atom_array",)}

    def _locate(self, a, tr, res, match):
        ref_res = self.ref[(self.ref.res_name == res) & ~self.ref.hetero]
        ref_xyz = tr.apply(np.array([ref_res.coord[ref_res.atom_name == n][0] for n in match]))
        best, rid = np.inf, None
        for r in np.unique(a.res_id[(a.res_name == res) & ~a.hetero]):
            x = a[(a.res_id == r) & (a.res_name == res) & ~a.hetero]
            try:
                xyz = np.array([x.coord[x.atom_name == n][0] for n in match])
            except IndexError:
                continue
            d = float(np.sqrt(((xyz - ref_xyz) ** 2).sum(1).mean()))
            if d < best:
                best, rid = d, r
        return rid, best

    def compute(self, atom_array: AtomArray) -> dict[str, float]:
        P = self.p
        a = atom_array[(atom_array.element != "H") & ~np.isnan(atom_array.coord).any(1)]
        lig, prot = a[a.hetero], a[~a.hetero]
        keys = ("base_rmsd", "hb_rmsd", "base_sasa_rel", "base_apolar_frac", "base_n_polar_env", "base_hbond_partners", "hb_partner_dist",
                "lig_hbonds", "nitro_hbonds", "nitro_unsat_buried", "motif_sasa_rel")
        nan = {k: float("nan") for k in keys}
        rl = self.ref[self.ref.hetero]
        names = [n for n in rl.atom_name if n in set(lig.atom_name)]
        if len(names) < 6 or len(prot) < 10:
            return nan
        A = np.array([lig.coord[lig.atom_name == n][0] for n in names]); B = np.array([rl.coord[rl.atom_name == n][0] for n in names])
        _, tr = superimpose(A, B)                      # maps reference -> design frame
        rb, db = self._locate(prot, tr, P["base_res"], P["base_match"]); rh, dh = self._locate(prot, tr, P["hb_res"], P["hb_match"])
        if rb is None or rh is None:
            return nan
        base = prot[(prot.res_id == rb) & (prot.res_name == P["base_res"])]
        bo = base[np.isin(base.atom_name, P["base_o"])]
        # SASA in context vs isolated residue
        sa = sasa(a, vdw_radii="ProtOr", point_number=200)
        iso = sasa(base, vdw_radii="ProtOr", point_number=200)
        in_ctx = sa[(a.res_id == rb) & (a.res_name == P["base_res"]) & ~a.hetero & np.isin(a.atom_name, P["base_o"])]
        iso_o = iso[np.isin(base.atom_name, P["base_o"])]
        base_rel = float(np.nansum(in_ctx) / max(np.nansum(iso_o), 1e-6))
        # environment of the base O atoms
        other = prot[prot.res_id != rb]; ot = cKDTree(other.coord)
        idx = np.unique(np.concatenate(ot.query_ball_point(bo.coord, P["env_r"]))).astype(int)
        env = other[idx]; polar = np.isin(env.element, ("N", "O")); apolar = np.isin(env.element, ("C", "S")) & ~np.isin(env.atom_name, ("C", "CA"))
        n_env = max(int(polar.sum() + apolar.sum()), 1)
        close = np.unique(np.concatenate(ot.query_ball_point(bo.coord, P["hb_cut"]))).astype(int)
        base_hb = int(np.isin(other.element[close], ("N", "O")).sum())
        # ligand H-bonds (protein N/O within cutoff)
        pol = prot[np.isin(prot.element, ("N", "O"))]; pt = cKDTree(pol.coord)
        lp = lig[np.isin(lig.element, ("N", "O"))]
        hb_n = np.array([len(pt.query_ball_point(c, P["hb_cut"])) for c in lp.coord])
        nit = np.isin(lp.atom_name, _NITRO)
        sa_l = sa[a.hetero][np.isin(lig.element, ("N", "O"))]
        unsat = int(((hb_n == 0) & nit & (sa_l < P["buried_sasa"])).sum())
        hb_res = prot[(prot.res_id == rh) & (prot.res_name == P["hb_res"]) & (prot.atom_name == P["hb_atom"])]
        hbd = float(np.linalg.norm(lp.coord - hb_res.coord[0], axis=1).min())
        sa_lig = sa[a.hetero]; iso_l = sasa(lig, vdw_radii="ProtOr", point_number=200)
        hb_sa = sa[(a.res_id == rh) & (a.res_name == P["hb_res"]) & ~a.hetero & (a.atom_name == P["hb_atom"])]
        iso_h = sasa(prot[(prot.res_id == rh) & (prot.res_name == P["hb_res"])], vdw_radii="ProtOr", point_number=200)
        hb_iso = iso_h[prot[(prot.res_id == rh) & (prot.res_name == P["hb_res"])].atom_name == P["hb_atom"]]
        motif_rel = np.nanmean(np.concatenate([[base_rel], hb_sa / np.maximum(hb_iso, 1e-6), sa_lig / np.maximum(iso_l, 1e-6)]))
        return {"base_rmsd": db, "hb_rmsd": dh, "base_sasa_rel": base_rel, "base_apolar_frac": float(apolar.sum() / n_env),
                "base_n_polar_env": float(polar.sum()), "base_hbond_partners": float(base_hb), "hb_partner_dist": hbd,
                "lig_hbonds": float((hb_n > 0).sum()), "nitro_hbonds": float(((hb_n > 0) & nit).sum()), "nitro_unsat_buried": float(unsat),
                "motif_sasa_rel": float(motif_rel)}
