"""Measure the 'ideal approach geometry' from the HG3-series crystals.

The hosts' slide scores base / H-bond-donor geometry "against the ideal approach
geometry" but gives no numbers, so the ideal is taken from every 6NT-bound
crystal copy that has the Asp127 base AND a Gln50 donor (the HG3.7+ family:
7K4Z, 4BS0, 7K4X, 8RD5, 7K4U). Windows are mean +/- max(sd, floor): eight copies
of one scaffold have tiny sd, and a de novo scaffold should not be punished for
being a few tenths of an Angstrom away from an evolved one.

Dihedrals and out-of-plane angles are scored by magnitude: their SIGN flips
between crystal copies (b_dih ranges -13..+11 deg) so the sign carries no
information about productivity, only the in-plane-ness does.

Writes configs/host_reference.yaml.   Usage: python -I calibrate_host_reference.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import yaml
from biotite.structure.io.pdbx import CIFFile, get_structure

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from host_geometry import BASE_METRICS, DONOR_METRICS, DIHEDRAL, burial_count, geometry  # noqa: E402
from kemp import geometry as g  # noqa: E402

FLOOR = {"dist": 0.25, "angle": 10.0}   # Angstrom / degrees
amap = g.load_atom_map(ROOT / "configs" / "atom_map.yaml")

rows, bur, src = [], [], []
for pdb in ("7K4Z", "4BS0", "7K4X", "8RD5", "7K4U"):
    s = get_structure(CIFFile.read(ROOT / "benchmark" / "pdb" / f"{pdb}.cif"), model=1)
    for lc, lr in g.iter_ligand_copies(s, amap):
        try:
            m = geometry(s, chain=lc, base_resid=127, donor_resid=50, lig_chain=lc, lig_resid=lr)
        except ValueError:
            continue
        for k in DIHEDRAL:
            m[k] = abs(m[k])
        rows.append(m); bur.append(burial_count(s, chain=lc, base_resid=127)); src.append(f"{pdb}/{lc}")

ref = {"source_copies": src, "n_copies": len(rows), "metrics": {}}
for k in BASE_METRICS + DONOR_METRICS:
    v = np.array([r[k] for r in rows])
    floor = FLOOR["dist"] if k.endswith("dist") else FLOOR["angle"]
    ref["metrics"][k] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
                         "floor": floor, "sigma": float(max(v.std(ddof=1), floor)),
                         "min": float(v.min()), "max": float(v.max())}
for key in ("burial_CG", "burial_O"):
    v = np.array([b[key] for b in bur])
    ref["metrics"][key] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
                           "min": float(v.min()), "max": float(v.max())}
ref["notes"] = ("ideal = crystal mean; geometry score per metric = exp(-0.5 z^2) with "
                "z = (x - mean)/sigma, sigma = max(sd, floor); dihedral/oop by magnitude; "
                "burial score = min(1, count / mean)")
out = ROOT / "configs" / "host_reference.yaml"
out.write_text(yaml.safe_dump(ref, sort_keys=False))
print(f"{len(rows)} crystal copies -> {out}")
print(f"{'metric':10}{'mean':>8}{'sd':>7}{'sigma':>8}   range")
for k, v in ref["metrics"].items():
    print(f"{k:10}{v['mean']:8.2f}{v['sd']:7.2f}{v.get('sigma', float('nan')):8.2f}   {v['min']:.1f} .. {v['max']:.1f}")
