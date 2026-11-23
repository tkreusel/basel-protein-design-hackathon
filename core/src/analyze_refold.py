"""Decompose refold outcome into the four things it actually measures.

The original refold gate asked one question -- "do the stage-1 criteria still
pass on the Boltz structure?" -- and answered 0/31. That single pass/fail
conflates four independent failures, and the conflation hid a good result:

  1. FOLD        does the sequence fold to the designed backbone?
                 (TM-score, design vs refold, protein only)
  2. POCKET      did Boltz put the ligand in the designed site?
                 (ligand centroid displacement after superposition, and
                  whether the ligand is buried at all)
  3. CONSTELLATION  are the three catalytic groups still arranged correctly
                 relative to each other? (9 tip atoms, superposed RMSD)
  4. ROTAMER     is the oxyanion donor side chain pointing the right way?
                 (donor tip displacement after superposition)

Measured on run 27: fold was excellent (25/31 at TM>=0.85) and the ligand was
buried in the designed pocket (rel_sasa ~0.000), while the donor rotamer was
wrong. Reporting that as "0/31 failed" was actively misleading.

Superposition uses the protein backbone only, so ligand and side-chain
displacements are independent of how well the fold matched.

Usage: python -I analyze_refold.py <project root> <stage1 run> <raw refold dir> [out.csv]
"""

from __future__ import annotations

import glob
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402
from kemp.record import read_manifest  # noqa: E402

TIPS = {"ASP": ["CG", "OD1", "OD2"],
        "GLN": ["CD", "OE1", "NE2"],
        "TRP": ["CD1", "NE1", "CE2"]}
DONOR_TIP = "NE2"
USALIGN = Path(__file__).resolve().parents[1] / "bin" / "USalign"


def _bb(struct, chain="A"):
    m = (~struct.hetero) & (struct.chain_id == chain) & (struct.atom_name == "CA")
    return struct[m]


def _kabsch_transform(mobile: np.ndarray, target: np.ndarray):
    """Return (rotation, mobile_centroid, target_centroid) aligning mobile->target."""
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    p, q = mobile - mc, target - tc
    v, _, wt = np.linalg.svd(p.T @ q)
    if np.linalg.det(v @ wt) < 0:
        v[:, -1] *= -1
    return (v @ wt), mc, tc


def _apply(xyz, rot, mc, tc):
    return (xyz - mc) @ rot + tc


def tm_score(design_pdb: Path, refold_pdb: Path) -> tuple[float, float]:
    if not USALIGN.exists():
        return float("nan"), float("nan")
    p = subprocess.run([str(USALIGN), str(refold_pdb), str(design_pdb)],
                       capture_output=True, text=True)
    tm = rms = float("nan")
    for line in p.stdout.splitlines():
        if line.startswith("TM-score=") and np.isnan(tm):
            try:
                tm = float(line.split("=")[1].split()[0])
            except Exception:
                pass
        if "RMSD=" in line and np.isnan(rms):
            try:
                rms = float(line.split("RMSD=")[1].split(",")[0])
            except Exception:
                pass
    return tm, rms


def tip_points(struct, cat, chain="A"):
    pts, names = [], []
    for role, rn in (("base", "ASP"), ("donor", "GLN"), ("wall", "TRP")):
        for an in TIPS[rn]:
            m = ((struct.chain_id == chain) & (struct.res_id == cat[role])
                 & (struct.atom_name == an) & (~struct.hetero))
            if not m.any():
                return None, None
            pts.append(struct.coord[m][0])
            names.append(f"{role}.{an}")
    return np.array(pts), names


def main(root: Path, stage1_run: Path, raw_refold: Path,
         out_csv: Path | None = None) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    recs = {r.design_id: r for r in read_manifest(stage1_run / "manifest.jsonl")}
    meta = json.loads((raw_refold / "designs.json").read_text())
    tmp = Path("/tmp") / f"kemp_tm_{raw_refold.name}"
    tmp.mkdir(parents=True, exist_ok=True)

    rows = []
    for name, m in meta.items():
        did = m["design_id"]
        # Arm B (LigandMPNN) ids are "<parent>_s<n>" -- the reference STRUCTURE
        # is the parent RFD4 backbone, since that is what the sequence was
        # designed onto. Arm A ids are their own parent.
        parent = m.get("parent_design", did)
        rec = recs.get(parent)
        hits = glob.glob(f"{raw_refold}/boltz/**/predictions/{name}/{name}_model_0.cif",
                         recursive=True)
        if not rec or not hits:
            continue
        ds = get_structure(CIFFile.read(rec.paths["rfd4_cif"]), model=1)
        rs = get_structure(CIFFile.read(hits[0]), model=1)

        dcat = {"base": rec.catalytic_residues["base"],
                "donor": rec.catalytic_residues["donor"],
                "wall": rec.catalytic_residues["wall"][0]}
        rcat = m["catalytic"]

        # --- 1. FOLD -----------------------------------------------------
        dp, rp = tmp / f"{did}_d.pdb", tmp / f"{did}_r.pdb"
        for p, st in ((dp, ds), (rp, rs)):
            f = PDBFile()
            f.set_structure(st[~st.hetero])
            f.write(p)
        tm, bb_rms = tm_score(dp, rp)

        # backbone superposition refold -> design
        dca, rca = _bb(ds), _bb(rs)
        n = min(dca.array_length(), rca.array_length())
        rot, mc, tc = _kabsch_transform(rca.coord[:n], dca.coord[:n])

        row = dict(design_id=did, stage1_rank=m["stage1_rank"],
                   length=m["length"], tm_score=tm, bb_rmsd=bb_rms)

        # Boltz's own confidence, so it can be correlated against the
        # structural outcomes rather than trusted blind.
        conf = Path(hits[0]).parent / f"confidence_{name}_model_0.json"
        if conf.exists():
            try:
                c = json.loads(conf.read_text())
                row["plddt"] = c.get("complex_plddt")
                row["iptm"] = c.get("iptm")
                row["ligand_iptm"] = c.get("ligand_iptm")
                row["ptm"] = c.get("ptm")
                row["complex_iplddt"] = c.get("complex_iplddt")
            except Exception:
                pass

        # --- 3. CONSTELLATION --------------------------------------------
        dpts, _ = tip_points(ds, dcat)
        rpts, _ = tip_points(rs, rcat)
        row["triad_rmsd"] = (g.kabsch_rmsd(rpts, dpts)
                             if dpts is not None and rpts is not None else np.nan)

        # --- 4. ROTAMER: donor tip displacement after bb superposition ----
        def tip(struct, cat, rn, an):
            mm = ((struct.chain_id == "A") & (struct.res_id == cat["donor"])
                  & (struct.atom_name == an) & (~struct.hetero))
            return struct.coord[mm][0] if mm.any() else None
        dt, rt = tip(ds, dcat, "GLN", DONOR_TIP), tip(rs, rcat, "GLN", DONOR_TIP)
        row["donor_tip_shift"] = (float(np.linalg.norm(_apply(rt, rot, mc, tc) - dt))
                                  if dt is not None and rt is not None else np.nan)

        # --- 2. POCKET: ligand centroid displacement ----------------------
        dl = ds[ds.res_name == amap.ligand_ccd_id]
        rl = rs[rs.res_name == amap.ligand_ccd_id]
        if dl.array_length() and rl.array_length():
            row["ligand_shift"] = float(np.linalg.norm(
                _apply(rl.coord.mean(axis=0), rot, mc, tc) - dl.coord.mean(axis=0)))
            cps = g.iter_ligand_copies(rs, amap)
            lc, lr = cps[0]
            try:
                b = g.ligand_burial(rs, amap, lig_chain=lc, lig_resid=lr)
                row["rf_rel_sasa"] = b["rel_sasa"]
                row["rf_enclosure"] = b["enclosure"]
            except g.MetricUnavailable:
                pass
            kw = dict(lig_chain=lc, lig_resid=lr, chain="A")
            for key, fn, kwarg in (("rf_base", g.base_to_abstraction, "base_resid"),
                                   ("rf_donor", g.donor_to_oxyanion, "donor_resid"),
                                   ("rf_wall", g.wall_stacking, "wall_resid")):
                try:
                    row[key] = fn(rs, amap, **{kwarg: rcat[key.split("_")[1]]}, **kw)
                except Exception:
                    row[key] = np.nan
        rows.append(row)

    if not rows:
        raise SystemExit("no refolds to analyse")

    def frac(key, op, thr):
        v = np.array([r.get(key, np.nan) for r in rows], float)
        v = v[~np.isnan(v)]
        return (int((v < thr).sum()) if op == "<" else int((v > thr).sum())), len(v)

    print(f"\n{'=' * 78}")
    print(f"REFOLD DECOMPOSED  --  {raw_refold.name}   n={len(rows)}")
    print("=" * 78)
    hdr = (f"{'design':14}{'rank':>5}{'TM':>7}{'triad':>8}{'donor_sh':>9}"
           f"{'lig_sh':>8}{'rf_base':>8}{'rf_donor':>9}{'sasa':>7}")
    print(hdr); print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: -(x["tm_score"] if x["tm_score"] == x["tm_score"] else 0)):
        f = lambda k, w, p=2: (f"{r[k]:{w}.{p}f}" if isinstance(r.get(k), float)
                               and r[k] == r[k] else f"{'-':>{w}}")  # noqa: E731
        print(f"{r['design_id']:14}{r['stage1_rank']:>5}{f('tm_score',7,3)}"
              f"{f('triad_rmsd',8)}{f('donor_tip_shift',9)}{f('ligand_shift',8)}"
              f"{f('rf_base',8)}{f('rf_donor',9)}{f('rf_rel_sasa',7,3)}")

    print("\n" + "=" * 78)
    print("WHAT ACTUALLY SURVIVED")
    print("=" * 78)
    for label, key, op, thr, why in (
        ("1. FOLD matches design (TM>=0.85)", "tm_score", ">", 0.85,
         "the sequence encodes the designed backbone"),
        ("2. ligand in designed pocket (<3A)", "ligand_shift", "<", 3.0,
         "Boltz found the right site"),
        ("2b. ligand buried (rel_sasa<0.15)", "rf_rel_sasa", "<", 0.15,
         "it is in a pocket at all"),
        ("3. constellation kept (triad<1.5A)", "triad_rmsd", "<", 1.5,
         "catalytic groups still arranged right"),
        ("4. donor rotamer kept (<1.5A)", "donor_tip_shift", "<", 1.5,
         "the Gln tip points the right way"),
    ):
        n, tot = frac(key, op, thr)
        pct = 100 * n / tot if tot else 0
        print(f"  {label:36} {n:>3}/{tot:<3} ({pct:>5.1f}%)  {why}")

    if out_csv:
        cols = sorted({k for r in rows for k in r})
        cols = ["design_id", "stage1_rank"] + [c for c in cols
                                               if c not in ("design_id", "stage1_rank")]
        with open(out_csv, "w") as fh:
            fh.write(",".join(cols) + "\n")
            for r in rows:
                fh.write(",".join(str(r.get(c, "")) for c in cols) + "\n")
        print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]),
         Path(sys.argv[4]) if len(sys.argv) > 4 else None)
