"""Compare design arm A (RFD4 co-design) against arm B (LigandMPNN).

Both arms are refolded with Boltz-2 from sequence alone, with 6NT bound. The
question this answers is not "which sequence scores better" but:

    does the SEQUENCE alone encode the active site?

RFD4 built the geometry, so measuring geometry on RFD4's own output is close to
circular. Measuring it on a Boltz-2 refold is not: Boltz sees only the sequence
and the ligand, so if the catalytic residues still land on the ligand in the
refold, the sequence genuinely specifies the site. That is the real comparison
between the two sequence-design strategies.

Reported per design/arm:
  geom.*        active-site geometry on RFD4's output structure
  refold.*      the same metrics on the Boltz-2 predicted structure
  boltz.*       Boltz confidence (plddt, iptm, ligand_iptm, ptm)
  refold_rmsd   not computed here; see plan 6.6 seed-consistency arm

Usage: python -I compare_arms.py <project root> <stage dir>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml
from biotite.structure.io.pdbx import CIFFile, get_structure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402
from kemp.record import read_manifest  # noqa: E402

# In the Boltz input YAML the protein is chain A and the ligand chain B.
REFOLD_PROT_CHAIN = "A"


def find_boltz_outputs(stage_dir: Path) -> dict[str, dict]:
    """Map Boltz prediction name -> {cif, confidence} paths."""
    out: dict[str, dict] = {}
    root = stage_dir / "boltz"
    if not root.exists():
        return out
    for pred_dir in root.rglob("predictions/*"):
        if not pred_dir.is_dir():
            continue
        name = pred_dir.name
        cifs = sorted(pred_dir.glob(f"{name}_model_*.cif"))
        confs = sorted(pred_dir.glob(f"confidence_{name}_model_*.json"))
        if cifs:
            out[name] = {
                "cif": cifs[0],
                "confidence": confs[0] if confs else None,
            }
    return out


def geometry_on_refold(cif: Path, amap: g.AtomMap, cat: dict) -> dict:
    """Our calibrated metrics on a Boltz-2 prediction.

    Boltz renumbers the protein 1..N, which is exactly the 1-indexed sequence
    positions we stored in the DesignRecord, so the catalytic positions carry
    over directly.
    """
    struct = get_structure(CIFFile.read(cif), model=1)
    copies = g.iter_ligand_copies(struct, amap)
    if not copies:
        raise g.MetricUnavailable(f"{cif.name}: ligand {amap.ligand_ccd_id} absent")
    lchain, lresid = copies[0]
    kw = dict(lig_chain=lchain, lig_resid=lresid, chain=REFOLD_PROT_CHAIN)
    wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat["wall"]
    return {
        "base_to_abstraction": g.base_to_abstraction(
            struct, amap, base_resid=cat["base"], **kw),
        "donor_to_oxyanion": g.donor_to_oxyanion(
            struct, amap, donor_resid=cat["donor"], **kw),
        "wall_stacking": g.wall_stacking(
            struct, amap, wall_resid=wall, **kw),
        "min_contact": g.min_protein_ligand_contact(
            struct, amap, lig_chain=lchain, lig_resid=lresid),
    }


def main(root: Path, stage_dir: Path) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    thr = g.load_thresholds(root / "configs" / "thresholds.yaml")
    geo_thr = thr["geometry"]

    records = {r.design_id: r for r in read_manifest(stage_dir / "manifest.jsonl")}
    boltz = find_boltz_outputs(stage_dir)
    index_path = stage_dir / "boltz_in" / "index.yaml"
    index = yaml.safe_load(index_path.read_text()) if index_path.exists() else []

    print(f"designs in manifest: {len(records)}")
    print(f"boltz predictions:   {len(boltz)}")
    if not index:
        print("no boltz_in/index.yaml -- did stage 5 run?")

    rows = []
    for item in index:
        name, arm, tag = item["name"], item["arm"], item["design"]
        did = tag.split("_mpnn")[0]
        rec = records.get(did)
        out = boltz.get(name)
        row = {
            "name": name, "arm": arm, "design": did,
            "length": item["length"],
            "geom_base": rec.scores.get("geom.base_to_abstraction") if rec else None,
            "geom_donor": rec.scores.get("geom.donor_to_oxyanion") if rec else None,
        }
        if out and out["confidence"] and out["confidence"].exists():
            c = json.loads(out["confidence"].read_text())
            row["plddt"] = c.get("complex_plddt")
            row["iptm"] = c.get("iptm")
            row["ligand_iptm"] = c.get("ligand_iptm")
            row["ptm"] = c.get("ptm")
        if out and rec and rec.catalytic_residues:
            try:
                rg = geometry_on_refold(out["cif"], amap, rec.catalytic_residues)
                row.update({f"refold_{k}": v for k, v in rg.items()})
            except Exception as exc:  # noqa: BLE001 - report, do not drop
                row["refold_error"] = f"{type(exc).__name__}: {exc}"
        rows.append(row)

    if not rows:
        raise SystemExit("nothing to compare")

    # --- table -----------------------------------------------------------
    hdr = (f"{'arm':5} {'design':13} {'len':>4} {'pLDDT':>6} {'lig_ipTM':>8} "
           f"{'rfd4_b':>7} {'rfd4_d':>7} {'refold_b':>8} {'refold_d':>8} {'site?':>6}")
    print("\n" + hdr)
    print("-" * len(hdr))

    def f(v, w=7, p=2):
        return f"{v:{w}.{p}f}" if isinstance(v, (int, float)) else f"{'-':>{w}}"

    pass_b = geo_thr["base_to_abstraction"]["pass_max"]
    pass_d = geo_thr["donor_to_oxyanion"]["pass_max"]

    for r in sorted(rows, key=lambda x: (x["arm"], x["design"])):
        rb, rd = r.get("refold_base_to_abstraction"), r.get("refold_donor_to_oxyanion")
        kept = (
            "YES" if isinstance(rb, float) and isinstance(rd, float)
            and rb <= pass_b and rd <= pass_d else
            ("no" if isinstance(rb, float) else "?")
        )
        print(f"{r['arm']:5} {r['design']:13} {r['length']:>4} "
              f"{f(r.get('plddt'), 6)} {f(r.get('ligand_iptm'), 8)} "
              f"{f(r.get('geom_base'))} {f(r.get('geom_donor'))} "
              f"{f(rb, 8)} {f(rd, 8)} {kept:>6}")
        if r.get("refold_error"):
            print(f"      refold_error: {r['refold_error']}")

    # --- per-arm summary -------------------------------------------------
    print("\n--- per-arm summary ---")
    for arm in sorted({r["arm"] for r in rows}):
        sub = [r for r in rows if r["arm"] == arm]
        def mean(key):
            vals = [r[key] for r in sub if isinstance(r.get(key), float)]
            return float(np.mean(vals)) if vals else float("nan")
        kept = [
            r for r in sub
            if isinstance(r.get("refold_base_to_abstraction"), float)
            and isinstance(r.get("refold_donor_to_oxyanion"), float)
            and r["refold_base_to_abstraction"] <= pass_b
            and r["refold_donor_to_oxyanion"] <= pass_d
        ]
        print(f"  {arm}: n={len(sub)}  mean pLDDT={mean('plddt'):.3f}  "
              f"mean lig_ipTM={mean('ligand_iptm'):.3f}  "
              f"mean refold base={mean('refold_base_to_abstraction'):.2f}  "
              f"donor={mean('refold_donor_to_oxyanion'):.2f}  "
              f"site kept in refold: {len(kept)}/{len(sub)}")

    out_csv = stage_dir / "comparison.csv"
    keys = sorted({k for r in rows for k in r})
    with open(out_csv, "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(str(r.get(k, "")) for k in keys) + "\n")
    print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
