"""Count polar side chains and cations around the catalytic base tip.

Rationale: a polar residue next to the carboxylate stabilises (and so
deactivates) the base by solvating or H-bonding it, and a cation (Lys/Arg)
within a few Angstrom forms a salt bridge that lowers its basicity. Both can
reduce activity, so we want to know how many each design has.

Base tip = CG, OD1, OD2 of the catalytic Asp (the organisers' Kemp.json tip
selection). A residue is counted if ANY qualifying side-chain atom is within
the cutoff of ANY tip atom. Backbone atoms never count.

  polar    side-chain polar atoms within 4.5 A:
           Ser OG, Thr OG1, Tyr OH, Asn OD1/ND2, Gln OE1/NE2,
           His ND1/NE2, Trp NE1
  cation   within 6.0 A (<=):
           Lys NZ, Arg NE/NH1/NH2.  His is reported SEPARATELY (its charge
           depends on protonation, which we cannot know).
  anion    Asp/Glu carboxylate O within 4.5 A, reported as extra context.

Excluded from every count: the base residue itself, the 6NT ligand, and the
catalytic donor Gln and wall Trp (designed to sit there; reported in their own
columns).

Usage: python -I base_environment.py <project root> <out csv>
"""

from __future__ import annotations

import csv
import glob
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp.record import read_manifest  # noqa: E402

TIP = ["CG", "OD1", "OD2"]
POLAR = {"SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"], "ASN": ["OD1", "ND2"],
         "GLN": ["OE1", "NE2"], "HIS": ["ND1", "NE2"], "TRP": ["NE1"]}
CATION = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"]}
ANION = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
R_POLAR, R_CATION = 4.5, 6.0

RUNS = {
    "A_full_RFD4": ("runs/2026-11-21_rfd4-backbones-BURIED_job21_32designs_seed201",
                    "runs/raw_full-motif-buried_job27",
                    "runs/2026-11-21_refold-full-motif-buried_job27/decomposed.csv"),
    "B_full_MPNN": ("runs/2026-11-21_rfd4-backbones-BURIED_job21_32designs_seed201",
                    "runs/raw_armB-mpnn_job32",
                    "runs/raw_armB-mpnn_job32/decomposed.csv"),
    "D_tip_MPNN": ("runs/2026-11-21_rfd4-backbones-TIPATOMS-BURIED_job23_32designs_seed202",
                   "runs/raw_tip-mpnn_job44",
                   "runs/raw_tip-mpnn_job44/decomposed.csv"),
    "C_tip_RFD4": ("runs/2026-11-21_rfd4-backbones-TIPATOMS-BURIED_job23_32designs_seed202",
                   "runs/raw_tip-atoms-buried_job28",
                   "runs/raw_tip-atoms-buried_job28/decomposed.csv"),
}


def fnum(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


def environment(struct, base_resid: int, donor_resid: int, wall_resid: int | None = None,
                chain="A") -> dict | None:
    prot = struct[(~struct.hetero) & (struct.chain_id == chain)]
    tipm = (prot.res_id == base_resid) & np.isin(prot.atom_name, TIP)
    if tipm.sum() < 3 or str(prot.res_name[tipm][0]) != "ASP":
        return None
    tip = prot.coord[tipm]

    def hits(table, cutoff):
        found = {}
        for rn, atoms in table.items():
            m = (prot.res_name == rn) & np.isin(prot.atom_name, atoms) \
                & (prot.res_id != base_resid)
            if not m.any():
                continue
            d = np.linalg.norm(prot.coord[m][:, None, :] - tip[None, :, :], axis=-1)
            dmin = d.min(axis=1)
            for rid, dd, an in zip(prot.res_id[m][dmin <= cutoff],
                                   dmin[dmin <= cutoff], prot.atom_name[m][dmin <= cutoff]):
                key = (int(rid), rn)
                if key not in found or dd < found[key][0]:
                    found[key] = (float(dd), str(an))
        return found

    polar = hits(POLAR, R_POLAR)
    # The catalytic donor Gln and wall Trp are DESIGNED to sit in the pocket
    # (the HG3.17 crystals have TRP44 at 3.2-3.3 A from the base tip), so they
    # are reported separately and not counted as unwanted polar neighbours.
    donor = {k: v for k, v in polar.items() if k[0] == donor_resid}
    wall = {k: v for k, v in polar.items() if wall_resid and k[0] == wall_resid}
    polar = {k: v for k, v in polar.items()
             if k[0] != donor_resid and not (wall_resid and k[0] == wall_resid)}
    allc = hits({**CATION, "HIS": POLAR["HIS"]}, R_CATION)
    cat = {k: v for k, v in allc.items() if k[1] in CATION}
    his6 = {k: v for k, v in allc.items() if k[1] == "HIS"}
    anion = hits(ANION, R_POLAR)

    fmt = lambda d: ";".join(f"{rn}{rid}@{v[0]:.1f}" for (rid, rn), v in sorted(d.items()))  # noqa: E731
    return {"n_polar": len(polar), "polar": fmt(polar),
            "n_cation": len(cat), "cation": fmt(cat),
            "n_his_6A": len(his6), "n_anion": len(anion),
            "donor_in_4p5": int(bool(donor)), "wall_in_4p5": int(bool(wall))}


def main(root: Path, out_csv: Path) -> None:
    out = []
    for arm, (s1, raw, dec) in RUNS.items():
        recs = {r.design_id: r for r in read_manifest(root / s1 / "manifest.jsonl")}
        meta = json.loads((root / raw / "designs.json").read_text())
        rows = {r["design_id"]: r for r in csv.DictReader(open(root / dec))}
        for name, m in meta.items():
            r = rows.get(m["design_id"])
            if r is None:
                continue
            tm, sasa = fnum(r.get("tm_score")), fnum(r.get("rf_rel_sasa"))
            triad, dsh = fnum(r.get("triad_rmsd")), fnum(r.get("donor_tip_shift"))
            promising = (tm is not None and tm >= 0.85 and sasa is not None
                         and sasa < 0.15 and ((triad is not None and triad < 1.5)
                                              or (dsh is not None and dsh < 1.5)))
            hits_ = glob.glob(f"{root}/{raw}/**/predictions/{name}/{name}_model_0.cif",
                              recursive=True)
            if not hits_:
                continue
            cat = m["catalytic"]
            row = dict(arm=arm, design_id=name, promising=int(promising), tm=tm,
                       triad=triad, donor_shift=dsh, sasa=sasa,
                       plddt=fnum(r.get("plddt")))
            env = environment(get_structure(CIFFile.read(hits_[0]), model=1),
                              cat["base"], cat["donor"], cat["wall"])
            if env is None:
                continue
            row.update({f"rf_{k}": v for k, v in env.items()})
            # Cross-check on the designed structure where the sequence is RFD4's
            parent = m.get("parent_design", m["design_id"])
            if m["design_id"] == parent and parent in recs:
                dcat = recs[parent].catalytic_residues
                denv = environment(get_structure(CIFFile.read(recs[parent].paths["rfd4_cif"]),
                                                 model=1), dcat["base"], dcat["donor"], dcat["wall"][0])
                if denv:
                    row.update({f"design_{k}": v for k, v in denv.items()})
            out.append(row)

    cols = sorted({k for r in out for k in r})
    first = ["arm", "design_id", "promising"]
    cols = first + [c for c in cols if c not in first]
    with open(out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in out:
            fh.write(",".join("" if r.get(c) is None else str(r.get(c, "")) for c in cols) + "\n")

    prom = [r for r in out if r["promising"]]
    print(f"measured {len(out)} refolds; {len(prom)} 'promising' "
          "(TM>=0.85, buried, and constellation OR donor <1.5 A)\n")
    for arm in RUNS:
        sub = [r for r in prom if r["arm"] == arm]
        print(f"{arm}: {len(sub)} promising")
    print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
