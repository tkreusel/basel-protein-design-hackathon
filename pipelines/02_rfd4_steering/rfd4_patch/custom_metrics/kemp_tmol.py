"""TMol protein<->ligand interaction energy that tolerates in-flight (search) structures.

TMol's atom-array canonicalisation trips over RFD4's conditioning annotations ("Coordinate group 0 must match absolute
coordinate conditioning") when the structure comes from a search checkpoint. Those annotations carry no chemistry, so they
are removed from a copy before scoring."""
from __future__ import annotations

import numpy as np
import torch
from biotite.structure import AtomArray

from rfproteina.metrics.physics.tmol_block_pair_metric import TMolBlockPairMetric


class KempTMolMetric(TMolBlockPairMetric):
    def compute(self, predicted_atom_array: AtomArray) -> dict[str, float]:
        a = predicted_atom_array.copy()
        for cat in list(a.get_annotation_categories()):
            if cat.startswith(("C_", "condition_")):
                a.del_annotation(cat)
        a.coord = np.asarray(a.coord, dtype=np.float32)
        # search runs inside bf16 autocast / a half-precision default dtype; TMol's packer needs float32/64
        prev = torch.get_default_dtype()
        torch.set_default_dtype(torch.float32)
        try:
            with torch.autocast(device_type="cuda", enabled=False):
                return super().compute(a)
        finally:
            torch.set_default_dtype(prev)
