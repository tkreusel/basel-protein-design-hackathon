"""Confidence metrics for our designs, joined to the structural verdicts.

Answers two things:
  1. How do the standard confidence metrics actually look (pLDDT, ipTM,
     ligand ipTM, pTM) and the directional ipSAE?
  2. Do they PREDICT whether the active site survived? This is the question
     that matters -- if confidence is uncorrelated with active-site survival,
     then ranking on confidence (what most pipelines do) selects for the
     wrong thing.

ipSAE direction: ipsae2 emits chn1->chn2 as `ipsae12`. Our Boltz inputs put
protein on chain A and ligand on chain B, so `ipsae12` is the
**protein->ligand** direction that KEMP_PIPELINE_PLAN.md section 6.6 asks for.
Plain `iptm` is reported alongside to show how much more optimistic it is --
the plan warns a single ligand atom can dominate it.

Usage: python -I confidence_report.py <decomposed.csv> <ipsae.tsv> [<ipsae.tsv> ...]
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path


def fnum(v):
    try:
        f = float(v)
        return f if f == f else None
    except (TypeError, ValueError):
        return None


def load_ipsae(paths: list[Path]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for p in paths:
        with open(p) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                # protein->ligand rows only
                if row.get("type1") != "prot" or row.get("type2") != "lig":
                    continue
                model = row["model"]
                # ".../predictions/<name>/<name>_model_0.cif" -> <name>
                name = Path(model).parent.name
                out[name] = {
                    "ipsae_p2l": fnum(row.get("ipsae12")),
                    "ipsae_l2p": fnum(row.get("ipsae21")),
                    "iptm_pair": fnum(row.get("iptm")),
                    "ipae_p2l": fnum(row.get("iPAE12")),
                    "n_contact": fnum(row.get("n_scored_contact12")),
                }
    return out


def spearman(xs, ys):
    n = len(xs)
    if n < 4:
        return None
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else None


def main(decomposed: Path, ipsae_paths: list[Path]) -> None:
    rows = list(csv.DictReader(open(decomposed)))
    ips = load_ipsae(ipsae_paths)
    for r in rows:
        r.update({k: v for k, v in ips.get(r["design_id"], {}).items()})

    print(f"designs: {len(rows)}   with ipSAE: "
          f"{sum(1 for r in rows if r.get('ipsae_p2l') is not None)}")

    print("\n" + "=" * 72)
    print("CONFIDENCE METRICS")
    print("=" * 72)
    print(f"{'metric':22}{'median':>9}{'min':>9}{'max':>9}  {'note':<22}")
    print("-" * 72)
    for key, label, note in (
        ("plddt", "pLDDT (complex)", "fold confidence"),
        ("iptm", "ipTM", "optimistic for ligands"),
        ("ligand_iptm", "ligand ipTM", "Boltz's own"),
        ("ipsae_p2l", "ipSAE protein->lig", "plan 6.6 metric"),
        ("ipsae_l2p", "ipSAE lig->protein", "the other direction"),
        ("ipae_p2l", "iPAE protein->lig", "Angstrom, lower better"),
    ):
        v = [fnum(r.get(key)) for r in rows]
        v = [x for x in v if x is not None]
        if not v:
            print(f"{label:22}{'-':>9}{'-':>9}{'-':>9}  {note:<22}")
            continue
        print(f"{label:22}{st.median(v):>9.3f}{min(v):>9.3f}{max(v):>9.3f}  {note:<22}")

    print("\n" + "=" * 72)
    print("DOES CONFIDENCE PREDICT ACTIVE-SITE SURVIVAL?")
    print("=" * 72)
    print("Spearman rho vs each structural outcome. Near 0 = the metric is")
    print("blind to whether the catalytic site survived.\n")
    targets = [("tm_score", "fold fidelity", +1),
               ("triad_rmsd", "constellation (lower better)", -1),
               ("donor_tip_shift", "donor rotamer (lower better)", -1),
               ("rf_rel_sasa", "ligand burial (lower better)", -1)]
    conf = ["plddt", "iptm", "ligand_iptm", "ipsae_p2l"]
    print(f"{'':30}" + "".join(f"{c[:11]:>13}" for c in conf))
    print("-" * 72)
    for tkey, tlabel, sign in targets:
        cells = []
        for ckey in conf:
            pairs = [(fnum(r.get(ckey)), fnum(r.get(tkey))) for r in rows]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
            rho = spearman([a for a, _ in pairs], [sign * b for _, b in pairs])
            cells.append(f"{rho:>+13.3f}" if rho is not None else f"{'-':>13}")
        print(f"{tlabel:30}" + "".join(cells))

    print("\n" + "=" * 72)
    print("TOP 10 BY ipSAE (protein->ligand), with what actually survived")
    print("=" * 72)
    hdr = (f"{'design':22}{'ipSAE':>7}{'ipTM':>7}{'pLDDT':>7}"
           f"{'TM':>7}{'triad':>7}{'donor':>7}{'sasa':>7}")
    print(hdr); print("-" * len(hdr))
    ranked = sorted((r for r in rows if fnum(r.get("ipsae_p2l")) is not None),
                    key=lambda r: -fnum(r["ipsae_p2l"]))
    for r in ranked[:10]:
        g = lambda k, w=7, p=2: (f"{fnum(r.get(k)):{w}.{p}f}"  # noqa: E731
                                 if fnum(r.get(k)) is not None else f"{'-':>{w}}")
        print(f"{r['design_id'][:22]:22}{g('ipsae_p2l',7,3)}{g('iptm',7,3)}"
              f"{g('plddt',7,3)}{g('tm_score',7,3)}{g('triad_rmsd')}"
              f"{g('donor_tip_shift')}{g('rf_rel_sasa',7,3)}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), [Path(p) for p in sys.argv[2:]])
