"""Join every per-design metric across all four arms and summarise.

Arms:
  A_full_RFD4   full-motif backbones, RFD4 co-designed sequence
  B_full_MPNN   same backbones, LigandMPNN sequences
  C_tip_RFD4    tip-atom backbones, RFD4 co-designed sequence
  D_tip_MPNN    same backbones, LigandMPNN sequences

Per design: refold decomposition (fold / pocket / constellation / donor), full
length CA scRMSD vs the designed backbone, protein->ligand ipSAE, and the
base-environment counts from base_environment.py.

Usage: python -I final_summary.py <project root> <out csv>
"""

from __future__ import annotations

import csv
import glob
import json
import statistics as st
import sys
import warnings
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from base_environment import RUNS, fnum  # noqa: E402
from kemp import geometry as g  # noqa: E402
from kemp.record import read_manifest  # noqa: E402

IPSAE = {"A_full_RFD4": ["armA_ipsae.tsv"],
         "B_full_MPNN": ["armB_s0_ipsae.tsv", "armB_s1_ipsae.tsv"],
         "C_tip_RFD4": ["armC_tipRFD4_ipsae.tsv"],
         "D_tip_MPNN": ["armD_tipMPNN_ipsae.tsv"]}


def ca(path):
    s = get_structure(CIFFile.read(path), model=1)
    return s[(~s.hetero) & (s.chain_id == "A") & (s.atom_name == "CA")].coord


def load_ipsae(root: Path, files):
    out = {}
    for f in files:
        for r in csv.DictReader(open(root / "results_ipsae" / f), delimiter="\t"):
            if r.get("type1") == "prot" and r.get("type2") == "lig":
                out[Path(r["model"]).parent.name] = fnum(r.get("ipsae12"))
    return out


def main(root: Path, out_csv: Path) -> None:
    env = {(r["arm"], r["design_id"]): r
           for r in csv.DictReader(open(root / "results_base_environment.csv"))}
    rows = []
    for arm, (s1, raw, dec) in RUNS.items():
        recs = {r.design_id: r for r in read_manifest(root / s1 / "manifest.jsonl")}
        meta = json.loads((root / raw / "designs.json").read_text())
        dr = {r["design_id"]: r for r in csv.DictReader(open(root / dec))}
        ips = load_ipsae(root, IPSAE[arm])
        for name, m in meta.items():
            d = dr.get(m["design_id"])
            if d is None:
                continue
            parent = m.get("parent_design", m["design_id"])
            hits = glob.glob(f"{root}/{raw}/**/predictions/{name}/{name}_model_0.cif",
                             recursive=True)
            sc = None
            if hits and parent in recs:
                a, b = ca(recs[parent].paths["rfd4_cif"]), ca(hits[0])
                if len(a) == len(b):
                    sc = g.kabsch_rmsd(b, a)
            e = env.get((arm, name), {})
            rows.append(dict(
                arm=arm, design_id=name, parent=parent,
                tm=fnum(d.get("tm_score")), scrmsd=sc,
                triad=fnum(d.get("triad_rmsd")), donor_shift=fnum(d.get("donor_tip_shift")),
                ligand_shift=fnum(d.get("ligand_shift")), sasa=fnum(d.get("rf_rel_sasa")),
                rf_base=fnum(d.get("rf_base")), rf_donor=fnum(d.get("rf_donor")),
                plddt=fnum(d.get("plddt")), iptm=fnum(d.get("iptm")),
                ipsae=ips.get(name),
                n_polar=fnum(e.get("rf_n_polar")), n_cation=fnum(e.get("rf_n_cation")),
                n_his6=fnum(e.get("rf_n_his_6A")), n_anion=fnum(e.get("rf_n_anion")),
                polar=e.get("rf_polar", ""), cation=e.get("rf_cation", ""),
            ))

    def ok(r, k, thr, lt=True):
        v = r.get(k)
        return v is not None and (v < thr if lt else v >= thr)

    for r in rows:
        r["fold_ok"] = ok(r, "tm", 0.85, False) and ok(r, "scrmsd", 2.0)
        r["buried"] = ok(r, "sasa", 0.15)
        r["const_ok"] = ok(r, "triad", 1.5)
        r["donor_ok"] = ok(r, "donor_shift", 1.5)
        r["strict"] = r["fold_ok"] and r["buried"] and r["const_ok"] and r["donor_ok"]
        r["promising"] = (r["fold_ok"] and r["buried"]
                          and (r["const_ok"] or r["donor_ok"]))
        r["clean_base"] = (r["n_cation"] == 0 and r["n_polar"] is not None
                           and r["n_polar"] <= 1)

    cols = list(rows[0])
    with open(out_csv, "w") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")

    arms = list(RUNS)
    print(f"{len(rows)} refolded sequences\n")
    print("FUNNEL per arm  (fold = TM>=0.85 AND scRMSD<2 A)")
    hdr = f"{'arm':14}{'n':>5}{'fold':>6}{'+bur':>6}{'+const':>7}{'+donor':>7}{'strict':>7}{'promis.':>8}"
    print(hdr); print("-" * len(hdr))
    for a in arms:
        s = [r for r in rows if r["arm"] == a]
        f1 = [r for r in s if r["fold_ok"]]
        b1 = [r for r in f1 if r["buried"]]
        c1 = [r for r in b1 if r["const_ok"]]
        d1 = [r for r in c1 if r["donor_ok"]]
        print(f"{a:14}{len(s):>5}{len(f1):>6}{len(b1):>6}{len(c1):>7}{len(d1):>7}"
              f"{sum(r['strict'] for r in s):>7}{sum(r['promising'] for r in s):>8}")

    print("\nMEDIANS per arm")
    print(f"{'arm':14}{'scRMSD':>8}{'pLDDT':>7}{'ipTM':>7}{'ipSAE':>7}{'donor_sh':>9}{'triad':>7}")
    for a in arms:
        s = [r for r in rows if r["arm"] == a]
        m = lambda k: st.median([r[k] for r in s if r[k] is not None])  # noqa: E731
        print(f"{a:14}{m('scrmsd'):>8.2f}{m('plddt'):>7.3f}{m('iptm'):>7.3f}"
              f"{m('ipsae'):>7.3f}{m('donor_shift'):>9.2f}{m('triad'):>7.2f}")

    print("\nBASE ENVIRONMENT of 'promising' designs (polar<=4.5A, cation<=6A of base tip)")
    print("crystal reference: HG3.17 / HG3.7 / HG3 all 0 polar (excl. wall Trp), 0 cation")
    for a in arms + ["ALL"]:
        s = [r for r in rows if r["promising"] and (a == "ALL" or r["arm"] == a)]
        if not s:
            continue
        print(f"  {a:12} n={len(s):>3}  0 cation: {sum(r['n_cation']==0 for r in s):>3}"
              f"  0 polar: {sum(r['n_polar']==0 for r in s):>3}"
              f"  clean(0 cation, <=1 polar): {sum(r['clean_base'] for r in s):>3}"
              f"  mean polar {np.mean([r['n_polar'] for r in s]):.2f}"
              f"  mean cation {np.mean([r['n_cation'] for r in s]):.2f}")
    print(f"\nwrote {out_csv}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
