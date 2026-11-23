"""Derive theozyme reference constraints from the crystal TSA complexes (HG3 lineage).

Writes config/theozyme_reference.yaml. These are *measured* values from the deposited 6NT
complexes (TSA, not the true TS) and are the default constraint set whenever no Rosetta-style
CST file is supplied (--cst). Tolerance = max(2 sd, floor); floors stop a 2-copy 'sd' of 0.01 A
from becoming an impossible constraint.
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd, yaml
from . import geometry as G

warnings.filterwarnings("ignore")
# metric key -> (kind, floor, force_k)
SPEC = {
    "geo.d_O_X":          ("linear", 0.25, 100.0),   # heavy-atom O...X (abstraction centre)
    "geo.d_O_H":          ("linear", 0.25, 100.0),
    "geo.ang_O_H_X":      ("linear", 15.0, 50.0),
    "geo.ang_C_O_H":      ("linear", 12.0, 50.0),
    "geo.tors_syn":       ("circular", 25.0, 30.0),
    "geo.oop_H":          ("linear", 15.0, 30.0),
    "geo.tors_C_O_H_X":   ("circular", 30.0, 20.0),
    "geo.tors_O_H_X_C3A": ("circular", 30.0, 20.0),
    "geo.wall_d_centroid": ("linear", 0.5, 20.0),
    "geo.wall_plane_angle": ("linear", 10.0, 10.0),
}
DONOR_SPEC = {"geo.donor_best_d": ("linear", 0.25, 100.0)}


def _stats(vals, kind):
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return None
    if kind == "circular":
        a = np.radians(v)
        mean = float(np.degrees(np.arctan2(np.sin(a).mean(), np.cos(a).mean())))
        d = np.array([G.wrap180(x - mean) for x in v])
        # circular: min/max are *deviations from the circular mean*
        return dict(mean=mean, sd=float(d.std()), dev_min=float(d.min()), dev_max=float(d.max()), n=len(v))
    return dict(mean=float(v.mean()), sd=float(v.std()), min=float(v.min()), max=float(v.max()), n=len(v))


def build(csv: Path, out: Path):
    df = pd.read_csv(csv, index_col=0)
    groups = {"hg3_lineage": df[df.index.str.startswith(("4BS0", "5RGA", "7K4Q", "8RD5"))],
              "kable2_5_nmr": df[df.index.str.startswith("29SB")],
              "hg3_with_gln50": df[df.index.str.startswith(("4BS0", "8RD5"))]}
    ref = {"provenance": f"derived from {csv.name}; TSA (6NT) crystal/NMR complexes, ideal-H geometry",
           "groups": {}}
    for g, sub in groups.items():
        spec = dict(SPEC)
        if g == "hg3_with_gln50":
            spec.update(DONOR_SPEC)
        ref["groups"][g] = {"n_copies": int(len(sub)), "constraints": {}}
        for k, (kind, floor, kf) in spec.items():
            s = _stats(sub[k], kind)
            if s is None:
                continue
            s.update(kind=kind, tol=float(max(2 * s["sd"], floor)), force_k=kf)
            ref["groups"][g]["constraints"][k] = {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in s.items()}
    out.write_text(yaml.safe_dump(ref, sort_keys=False))
    return ref


if __name__ == "__main__":
    r = build(Path(sys.argv[1]), Path(sys.argv[2]))
    for g, d in r["groups"].items():
        print(g, d["n_copies"])
        for k, s in d["constraints"].items():
            print(f"   {k:24s} mean {s['mean']:8.2f}  sd {s['sd']:6.2f}  tol {s['tol']:6.2f}  [{s.get('min', s.get('dev_min')):.2f}, {s.get('max', s.get('dev_max')):.2f}] {s['kind']}")
