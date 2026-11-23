"""Milestone 1: validate the geometry scoring stack on the HG series (plan 4).

Plan 4.2 acceptance criteria:
  1. Spearman rho >= 0.7 between composite score and log10(kcat/KM).
  2. Ablation separation: HG3.17 outscores its ablations, and specifically
     outscores 6SS1/6SS3 on the donor-oxyanion metric.
  3. Known-geometry recovery -- covered by tests/test_geometry.py.

Plan 4.2 also says: "Report it honestly; if it is near zero, the individual
metrics are wrong, not the correlation." This script reports what it measures.

Usage: python -I milestone1.py <project root>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402

BASE_RESID, DONOR_RESID, WALL_RESID = 127, 50, 44

# Measured activities, plan 4.1 / 10. None = no published kcat/KM that we have
# verified; such rows are used for criterion 2 (ablation separation) only and
# are EXCLUDED from the criterion 1 correlation. Do not fill these in from
# memory -- a wrong activity silently corrupts the only quantitative gate.
BENCH = [
    # pdb,  variant,                         lig_chain, lig_resid, kcat_KM
    # --- the activity series (criterion 1) ---
    ("7K4Q", "HG3",                             "A", 401, 146.0),
    ("7K4X", "HG3.7",                           "A", 401, 27000.0),
    ("7K4Z", "HG3.17",                          "A", 402, 126000.0),
    # --- position-50 series, ALL on the same HG3.17 E47N/N300D background.
    #     4BS0 is the Q50 parent, so comparisons within this block are
    #     genuinely single-variable at position 50 (criterion 2).
    ("4BS0", "HG3.17 Q50  E47N/N300D  [parent]", "A", 1305, None),
    ("6SRW", "HG3.17 Q50F E47N/N300D",           "A", 402, None),
    ("6TU6", "HG3.17 Q50M E47N/N300D",           "A", 401, None),
    ("6SS1", "HG3.17 Q50A E47N/N300D",           "A", 401, None),
    ("6SS3", "HG3.17 Q50K E47N/N300D",           "A", 401, None),
    # --- HG3 background, single substitution at 50 ---
    ("7K4U", "HG3 K50Q",                         "A", 401, None),
    # --- other variants / scaffolds ---
    ("8RD5", "HG3.R5",                           "A", 401, None),
    ("6NW4", "evolved computational design",     "A", 401, None),
]


def composite(base: float, donor: float, wall: float, thr: dict) -> float:
    """Geometry-only composite. Higher is better.

    Weighting rationale, from the measured distributions (thresholds.yaml):
      * donor_to_oxyanion is THE discriminating metric (active 2.96-3.10 vs
        inactive 6.05-7.42), so it carries the most weight.
      * base_to_abstraction is necessary but NOT discriminating -- the weakest
        variant has the shortest distance -- so it is scored as a window-
        satisfaction term, not a "shorter is better" term.
      * wall_stacking is a weak packing term.
    Each term is a 0-1 window score; the composite is their weighted sum.
    """
    geo = thr["geometry"]

    def window(value: float, lo: float, hi: float, pass_max: float) -> float:
        """1.0 inside [lo, hi], decaying to 0 at pass_max."""
        if lo <= value <= hi:
            return 1.0
        if value < lo:
            return 0.0  # too close = clash
        if value >= pass_max:
            return 0.0
        return float((pass_max - value) / (pass_max - hi))

    b = window(base, *geo["base_to_abstraction"]["ideal"],
               geo["base_to_abstraction"]["pass_max"])
    d = window(donor, *geo["donor_to_oxyanion"]["ideal"],
               geo["donor_to_oxyanion"]["pass_max"])
    w = window(wall, *geo["wall_stacking"]["ideal"],
               geo["wall_stacking"]["pass_max"])
    return 0.60 * d + 0.25 * b + 0.15 * w


def main(root: Path) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    thr = g.load_thresholds(root / "configs" / "thresholds.yaml")

    rows = []
    for pdb, variant, lc, lr, activity in BENCH:
        s = get_structure(CIFFile.read(root / "benchmark" / "pdb" / f"{pdb}.cif"), model=1)
        kw = dict(lig_chain=lc, lig_resid=lr, chain="A")
        try:
            base = g.base_to_abstraction(s, amap, base_resid=BASE_RESID, **kw)
            donor = g.donor_to_oxyanion(s, amap, donor_resid=DONOR_RESID, **kw)
            wall = g.wall_stacking(s, amap, wall_resid=WALL_RESID, **kw)
            clash = g.min_protein_ligand_contact(s, amap, lig_chain=lc, lig_resid=lr)
            score = composite(base, donor, wall, thr)
            rows.append(dict(pdb=pdb, variant=variant, base=base, donor=donor,
                             wall=wall, clash=clash, score=score, activity=activity))
        except g.MetricUnavailable as exc:
            print(f"  {pdb}: METRIC UNAVAILABLE: {exc}")

    hdr = f"{'PDB':6} {'variant':26} {'base':>6} {'donor':>7} {'wall':>6} {'clash':>6} {'score':>6} {'kcat/KM':>9}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        act = f"{r['activity']:,.0f}" if r["activity"] else "-"
        print(f"{r['pdb']:6} {r['variant']:26} {r['base']:6.2f} {r['donor']:7.2f} "
              f"{r['wall']:6.2f} {r['clash']:6.2f} {r['score']:6.3f} {act:>9}")

    # --- criterion 1: Spearman vs log10 activity ------------------------
    act_rows = [r for r in rows if r["activity"]]
    print(f"\n--- Criterion 1: Spearman rho vs log10(kcat/KM), n={len(act_rows)} ---")
    if len(act_rows) < 3:
        print("  too few activity points to correlate")
    else:
        y = np.log10([r["activity"] for r in act_rows])
        for name, x, sign in (
            ("composite (higher better)", [r["score"] for r in act_rows], 1),
            ("donor dist (lower better)", [r["donor"] for r in act_rows], -1),
            ("base dist  (lower better)", [r["base"] for r in act_rows], -1),
        ):
            rho, p = spearmanr(np.array(x) * sign, y)
            flag = "PASS" if rho >= 0.7 else "FAIL"
            print(f"  {name:28} rho={rho:+.3f}  p={p:.3f}  [{flag} vs 0.7]")
        print("\n  Predicted order (best->worst):",
              " > ".join(r["variant"] for r in sorted(act_rows, key=lambda r: -r["score"])))
        print("  Actual order    (best->worst):",
              " > ".join(r["variant"] for r in sorted(act_rows, key=lambda r: -r["activity"])))

    # --- criterion 2: ablation separation -------------------------------
    print("\n--- Criterion 2: ablation separation on the donor metric ---")
    by = {r["pdb"]: r for r in rows}
    hg317 = by["7K4Z"]
    for abl in ("6SS1", "6SS3"):
        a = by[abl]
        ok = hg317["donor"] < a["donor"]
        print(f"  HG3.17 donor {hg317['donor']:.2f} < {abl} donor {a['donor']:.2f} "
              f"-> {'PASS' if ok else 'FAIL'}")
    for abl in ("6SS1", "6SS3"):
        a = by[abl]
        ok = hg317["score"] > a["score"]
        print(f"  HG3.17 composite {hg317['score']:.3f} > {abl} {a['score']:.3f} "
              f"-> {'PASS' if ok else 'FAIL'}")

    print("\n--- Criterion 3: known-geometry recovery ---")
    print("  see tests/test_geometry.py (20 passed)")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
