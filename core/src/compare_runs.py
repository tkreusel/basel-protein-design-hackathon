"""Compare finalized run folders side by side.

Built to answer two questions:

  1. Does conditioning the ligand as BURIED raise the pocket hit rate above
     the 37.5% the unconditioned spec achieved?
  2. Does supplying only tip atoms (letting the model infer rotamers and
     backbone placement) beat supplying whole residues?

Reads each run's manifest, so it reports exactly what the gate decided and
cannot disagree with the run's own summary.csv.

Usage: python -I compare_runs.py <run dir> [<run dir> ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp.record import read_manifest  # noqa: E402

GATE = "stage1_motif"
# Crystal reference, measured by src/check_burial.py on 7K4Z/4BS0/7K4X.
CRYSTAL = {"rel_sasa": 0.008, "enclosure": 0.996, "n_contacts": 66.0,
           "base": 2.53, "donor": 3.10, "wall": 5.53}


def stats(run: Path) -> dict | None:
    mf = run / "manifest.jsonl"
    if not mf.exists():
        return None
    recs = list(read_manifest(mf))
    if not recs:
        return None
    ok = [r for r in recs if r.status == "ok" and r.gates.get(GATE) == "pass"]

    def col(rs, key):
        v = [r.scores.get(key) for r in rs]
        return np.array([x for x in v if isinstance(x, (int, float))], float)

    # Why each failure happened, so the two experiments are comparable on
    # mechanism and not just on a single pass rate.
    reasons: dict[str, int] = {}
    for r in recs:
        if r.status == "ok" and r.gates.get(GATE) == "pass":
            continue
        rl = r.gates.get(f"{GATE}_reasons") or ["unknown"]
        t = rl[0]
        k = ("ligand-on-surface" if any(s in t for s in
             ("rel_sasa", "enclosure", "contacts"))
             else "base-too-far" if "base->" in t
             else "donor-too-far" if "donor->" in t
             else "wall-too-far" if "wall" in t
             else "too-long" if "length" in t
             else "other")
        reasons[k] = reasons.get(k, 0) + 1

    all_sasa = col(recs, "geom.rel_sasa")
    return {
        "run": run.name,
        "n": len(recs),
        "n_ok": len(ok),
        "rate": len(ok) / len(recs) if recs else 0.0,
        "reasons": reasons,
        # geometry over ALL attempts: did the motif survive regardless of pocket?
        "base_all": col(recs, "geom.base_to_abstraction"),
        "donor_all": col(recs, "geom.donor_to_oxyanion"),
        "wall_all": col(recs, "geom.wall_stacking"),
        "sasa_all": all_sasa,
        "encl_all": col(recs, "geom.enclosure"),
        # and over the successes only
        "sasa_ok": col(ok, "geom.rel_sasa"),
        "contacts_ok": col(ok, "geom.n_contacts"),
        "len_ok": np.array([r.length for r in ok if r.length], float),
    }


def rng(a):
    return f"{a.min():.2f}-{a.max():.2f}" if a.size else "-"


def main(runs: list[Path]) -> None:
    rows = [s for s in (stats(r) for r in runs) if s]
    if not rows:
        raise SystemExit("no finalized runs found (need manifest.jsonl)")

    print("=" * 78)
    print("POCKET HIT RATE  (success = active site survived AND ligand buried)")
    print("=" * 78)
    print(f"{'run':46} {'n':>3} {'pass':>5} {'rate':>7}")
    print("-" * 78)
    for s in rows:
        print(f"{s['run'][:46]:46} {s['n']:>3} {s['n_ok']:>5} {s['rate']*100:>6.1f}%")

    print("\n" + "=" * 78)
    print("WHY THE FAILURES FAILED")
    print("=" * 78)
    keys = sorted({k for s in rows for k in s["reasons"]})
    print(f"{'run':46} " + " ".join(f"{k[:17]:>18}" for k in keys))
    print("-" * 78)
    for s in rows:
        print(f"{s['run'][:46]:46} " +
              " ".join(f"{s['reasons'].get(k, 0):>18}" for k in keys))

    print("\n" + "=" * 78)
    print("DID THE MOTIF SURVIVE?  (all attempts, not just successes)")
    print("=" * 78)
    print(f"{'run':46} {'base':>12} {'donor':>12} {'wall':>12}")
    print("-" * 78)
    for s in rows:
        print(f"{s['run'][:46]:46} {rng(s['base_all']):>12} "
              f"{rng(s['donor_all']):>12} {rng(s['wall_all']):>12}")
    print(f"{'HG3.17 crystal (7K4Z)':46} {CRYSTAL['base']:>12.2f} "
          f"{CRYSTAL['donor']:>12.2f} {CRYSTAL['wall']:>12.2f}")

    print("\n" + "=" * 78)
    print("BURIAL  (mean over all attempts; and over successes only)")
    print("=" * 78)
    print(f"{'run':46} {'sasa_all':>9} {'encl_all':>9} {'sasa_ok':>8} {'cont_ok':>8}")
    print("-" * 78)
    for s in rows:
        m = lambda a: f"{a.mean():.3f}" if a.size else "-"  # noqa: E731
        mi = lambda a: f"{a.mean():.1f}" if a.size else "-"  # noqa: E731
        print(f"{s['run'][:46]:46} {m(s['sasa_all']):>9} {m(s['encl_all']):>9} "
              f"{m(s['sasa_ok']):>8} {mi(s['contacts_ok']):>8}")
    print(f"{'crystals (7K4Z/4BS0/7K4X)':46} {CRYSTAL['rel_sasa']:>9.3f} "
          f"{CRYSTAL['enclosure']:>9.3f} {CRYSTAL['rel_sasa']:>8.3f} "
          f"{CRYSTAL['n_contacts']:>8.1f}")

    print("\nlower sasa = better buried; higher enclosure = better surrounded")


if __name__ == "__main__":
    main([Path(p) for p in sys.argv[1:]])
