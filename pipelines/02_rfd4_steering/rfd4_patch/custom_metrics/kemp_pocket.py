"""Pocket volume match and active-site loop content around a fixed ligand (Kemp project).

Outputs (flat floats, prefixed by the config key):
  pocket_volume       A^3 of grid points (0.7 A grid) that are free of protein heavy atoms (>= clearance from any protein atom),
                      enclosed (ray test) and within `shell_r` of the ligand. The ligand's own space counts as pocket.
  ligand_volume       A^3 of grid points within `lig_vdw` of a ligand heavy atom (same grid)
  volume_ratio        pocket_volume / ligand_volume
  volume_ratio_error  |volume_ratio - target_ratio|   (minimise)
  ligand_sasa_frac    fraction of ligand surface (Shrake-Rupley, probe 1.4) not covered by protein heavy atoms
  n_shell             protein residues with a heavy atom within `shell_r` of the ligand
  loop_frac_shell     fraction of those residues annotated coil by P-SEA (CA trace) -- coarse, intended for ranking
"""
from __future__ import annotations

import numpy as np
from biotite.structure import AtomArray, annotate_sse
from scipy.spatial import cKDTree

from rfproteina.metrics.metric import Metric


def _dirs(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


class KempPocketMetric(Metric):
    def __init__(self, target_ratio: float = 1.0, shell_r: float = 8.0, clearance: float = 3.0, grid: float = 0.7,
                 lig_vdw: float = 1.7, enclosure_min: float = 0.75, tol: float = 0.0, **kwargs):
        super().__init__(**kwargs)
        self.target_ratio, self.shell_r, self.clearance, self.tol = target_ratio, shell_r, clearance, tol
        self.grid, self.lig_vdw, self.enclosure_min = grid, lig_vdw, enclosure_min
        self._d = _dirs(40)

    @property
    def kwargs_to_compute_args(self) -> dict:
        return {"atom_array": ("predicted_atom_array",)}

    def compute(self, atom_array: AtomArray) -> dict[str, float]:
        a = atom_array[(atom_array.element != "H") & ~np.isnan(atom_array.coord).any(1)]
        lig, prot = a[a.hetero], a[~a.hetero]
        nan = float("nan")
        if len(lig) == 0 or len(prot) < 10:
            return {k: nan for k in ("pocket_volume", "ligand_volume", "volume_ratio", "volume_ratio_error", "volume_ratio_band_error",
                                     "ligand_sasa_frac", "n_shell", "loop_frac_shell")}
        pt, lt = cKDTree(prot.coord), cKDTree(lig.coord)
        c = lig.coord.mean(0)
        r = self.shell_r + float(np.linalg.norm(lig.coord - c, axis=1).max())
        ax = np.arange(-r, r + 1e-6, self.grid)
        g = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3) + c
        dl, _ = lt.query(g)
        g, dl = g[dl <= self.shell_r], dl[dl <= self.shell_r]
        dp, _ = pt.query(g)
        free = dp >= self.clearance
        g, dl = g[free], dl[free]
        # enclosure: fraction of 40 rays hitting a protein atom (within 1.8 A of the ray) before 14 A
        steps = np.arange(1.5, 14.0, 1.0)
        enc = np.zeros(len(g))
        for d in self._d:
            p = g[:, None, :] + steps[None, :, None] * d[None, None, :]
            hit = (pt.query(p.reshape(-1, 3))[0].reshape(len(g), -1) < 1.8).any(1)
            enc += hit
        enc /= len(self._d)
        vox = self.grid ** 3
        pocket = float((enc >= self.enclosure_min).sum() * vox)
        # ligand volume on the same grid
        gl = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3) + c
        lig_vol = float((lt.query(gl)[0] <= self.lig_vdw).sum() * vox)
        ratio = pocket / lig_vol if lig_vol > 0 else nan
        # ligand SASA fraction
        pts = lig.coord[:, None, :] + (self.lig_vdw + 1.4) * _dirs(64)[None]
        free_pt = pt.query(pts.reshape(-1, 3))[0].reshape(len(lig), -1) > (1.7 + 1.4)
        # exclude points buried in the ligand itself
        own = lt.query(pts.reshape(-1, 3))[0].reshape(len(lig), -1) > (self.lig_vdw + 1.4 - 1e-3)
        sasa = float((free_pt & own).sum() / max(own.sum(), 1))
        # loops in the shell
        shell_res = set()
        for i in np.unique(np.concatenate(pt.query_ball_point(lig.coord, self.shell_r))).astype(int):
            shell_res.add((prot.chain_id[i], int(prot.res_id[i])))
        loop = nan
        try:
            sse = annotate_sse(prot)
            ca = prot[prot.atom_name == "CA"]
            keys = list(zip(ca.chain_id, ca.res_id.astype(int)))
            m = np.array([k in shell_res for k in keys])
            loop = float((np.array(sse)[m] == "c").mean()) if m.any() else nan
        except Exception:
            pass
        return {"pocket_volume": pocket, "ligand_volume": lig_vol, "volume_ratio": ratio,
                "volume_ratio_error": abs(ratio - self.target_ratio) if ratio == ratio else nan,
                "volume_ratio_band_error": max(0.0, abs(ratio - self.target_ratio) - self.tol) if ratio == ratio else nan,
                "ligand_sasa_frac": sasa, "n_shell": float(len(shell_res)), "loop_frac_shell": loop}
