"""Solvent/substrate channel from the ligand pocket to bulk (ligand removed: it is replaced by substrate).

Grid (`grid` A) in a cube of half-width `box` A around the ligand; a point is free for probe radius r if it is >= 1.7 + r from every protein heavy atom.
Start = free points within `start_r` of any ligand atom; a channel of radius r exists if that region is 26-connected to the cube boundary.
Outputs: bottleneck_radius (largest tested r with a path; 0 if none even at the smallest), channel_band_error = max(0, lo - r) + max(0, r - hi)  (lo=water 1.4)."""
from __future__ import annotations

import numpy as np
from biotite.structure import AtomArray
from scipy import ndimage
from scipy.spatial import cKDTree

from rfproteina.metrics.metric import Metric


class KempChannelMetric(Metric):
    def __init__(self, lo: float = 1.4, hi: float = 1.9, grid: float = 0.8, box: float = 24.0, start_r: float = 3.0,
                 radii: tuple = (0.5, 0.7, 0.9, 1.1, 1.25, 1.4, 1.55, 1.7, 1.9, 2.1), **kwargs):
        super().__init__(**kwargs)
        self.lo, self.hi, self.grid, self.box, self.start_r, self.radii = lo, hi, grid, box, start_r, tuple(float(x) for x in radii)

    @property
    def kwargs_to_compute_args(self) -> dict:
        return {"atom_array": ("predicted_atom_array",)}

    def compute(self, atom_array: AtomArray) -> dict[str, float]:
        a = atom_array[(atom_array.element != "H") & ~np.isnan(atom_array.coord).any(1)]
        lig, p = a[a.hetero], a[~a.hetero]
        if len(lig) == 0 or len(p) < 10:
            return {"bottleneck_radius": float("nan"), "channel_band_error": float("nan")}
        c = lig.coord.mean(0); ax = np.arange(-self.box, self.box + 1e-6, self.grid) + 0.0
        G = np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1) + c; sh = G.shape[:3]; pts = G.reshape(-1, 3)
        dp = cKDTree(p.coord).query(pts, distance_upper_bound=4.0)[0].reshape(sh); dl = cKDTree(lig.coord).query(pts)[0].reshape(sh)
        best = 0.0
        for r in self.radii:
            free = dp >= 1.7 + r; st = free & (dl <= self.start_r)
            if not st.any():
                break
            lab, _ = ndimage.label(free, structure=np.ones((3, 3, 3)))
            ids = set(np.unique(lab[st])) - {0}
            edge = set(np.unique(np.concatenate([lab[0].ravel(), lab[-1].ravel(), lab[:, 0].ravel(), lab[:, -1].ravel(), lab[:, :, 0].ravel(), lab[:, :, -1].ravel()]))) - {0}
            if ids & edge:
                best = r
            else:
                break
        return {"bottleneck_radius": best, "channel_band_error": max(0.0, self.lo - best) + max(0.0, best - self.hi)}
