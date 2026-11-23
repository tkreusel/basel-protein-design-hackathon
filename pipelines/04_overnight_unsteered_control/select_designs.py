"""Stage 2 of the overnight run: measure every unsteered RFD4 design and pick
the ones that go on to LigandMPNN + Boltz.

Criteria follow Maxi's `select_passers.py` / `eval_r5.py` where they can be
computed faithfully here (base-environment metrics are our reimplementation,
validated 198/198 against his published numbers):

  lig_sasa     ligand SASA fraction            <= 0.10
  base_buried  base rel SASA                   <= 0.10
  no_sc_polar  base polar side-chain atoms     <= 1   (with-Trp: the fixed Trp NE1)
  no_cation    Lys/Arg within 6 A of base O    == 0
  nitro        nitro H-bond far >=1 and none near the base
  clash        ligand clash atoms              == 0
  motif        base AND partner tip RMSD       <= 0.3 A  (ligand-anchored locator)

NOT applied here (need his custom metrics / TMol, not run in this pipeline):
pocket volume ratio 0.4-1.4, pp clash fraction, TMol dG < 0. This is stated in
every output so the pass counts are not mistaken for his full `pass_core`.

Selection modes (decided with Maxi's session, see overnight/README.md):
  all        every design
  topk       the K designs passing the most criteria, ties broken by fewer
             base polar atoms + cations, then lower base SASA
  core       only designs passing every applied criterion

Usage: python -I select_designs.py <raw rfd4 out dir> <out dir> [--mode topk --k 100]
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import shutil
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from motif_locate import locate  # noqa: E402
from steering_criteria import measure  # noqa: E402

MOTIF_CIF = ROOT / "inputs" / "maxi" / "kemp_motif.cif"
ROLES = {  # Maxi's tip sets (full Trp ring, as in his fixed-residue list)
    "base": dict(resid=127, resname="ASP", tips=["CB", "CG", "OD1", "OD2"]),
    "partner": dict(resid=50, resname="GLN", tips=["CG", "CD", "OE1", "NE2"]),
    "wall": dict(resid=44, resname="TRP",
                 tips=["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"]),
}
CRITERIA = ["lig_sasa", "base_buried", "no_sc_polar", "no_cation", "nitro", "clash", "motif"]


def _one(path: str) -> dict:
    motif = get_structure(CIFFile.read(MOTIF_CIF), model=1)
    s = get_structure(CIFFile.read(path), model=1)
    r = {"cif": path, "design": Path(path).parent.name}
    try:
        m = measure(s)
        loc = locate(s, motif, ROLES, tol=8.0)
        r.update(m)
        r.update(length=len({int(x) for x in s.res_id[~s.hetero]}),
                 rmsd_base=loc["base"]["rmsd"], rmsd_partner=loc["partner"]["rmsd"],
                 rmsd_wall=loc["wall"]["rmsd"], base_resid=loc["base"]["resid"],
                 partner_resid=loc["partner"]["resid"], wall_resid=loc["wall"]["resid"])
        r["ok_lig_sasa"] = m["lig_sasa_frac"] <= 0.10
        r["ok_base_buried"] = m["base_sasa_rel"] <= 0.10
        r["ok_no_sc_polar"] = m["base_sc_polar_4p5"] <= 1
        r["ok_no_cation"] = m["base_cation_6"] == 0
        r["ok_nitro"] = m["nitro_hb_far"] >= 1 and m["nitro_hb_near_base"] == 0
        r["ok_clash"] = m["lig_clash"] == 0
        r["ok_motif"] = r["rmsd_base"] <= 0.3 and r["rmsd_partner"] <= 0.3
        r["n_pass"] = sum(bool(r[f"ok_{c}"]) for c in CRITERIA)
        r["pass_applied"] = r["n_pass"] == len(CRITERIA)
        r["status"] = "ok"
    except Exception as exc:  # noqa: BLE001 - record, never drop
        r["status"] = f"failed: {type(exc).__name__}: {exc}"
    return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("raw"); ap.add_argument("out")
    ap.add_argument("--mode", default="topk", choices=["all", "topk", "core"])
    ap.add_argument("--k", type=int, default=100)
    ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    out = Path(a.out); (out / "designs_in").mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(f"{a.raw}/inference_samples/*/*_generated.cif"))
    print(f"{len(files)} generated designs under {a.raw}")
    with Pool(a.procs) as p:
        rows = p.map(_one, files)
    ok = [r for r in rows if r["status"] == "ok"]
    print(f"measured {len(ok)}, failed {len(rows) - len(ok)}")

    ok.sort(key=lambda r: (-r["n_pass"], r["base_sc_polar_4p5"] + r["base_cation_6"],
                           r["base_sasa_rel"]))
    if a.mode == "all":
        sel = ok
    elif a.mode == "core":
        sel = [r for r in ok if r["pass_applied"]]
    else:
        sel = ok[: a.k]

    cols = sorted({k for r in rows for k in r})
    with open(out / "all_designs_metrics.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)
    for i, r in enumerate(sel, 1):
        r["select_rank"] = i
        shutil.copy2(r["cif"], out / "designs_in" / f"{r['design']}.cif")
    with open(out / "selection.csv", "w") as fh:
        c2 = ["select_rank", "design", "cif", "length", "n_pass", "base_resid",
              "partner_resid", "wall_resid"] + sorted(
            k for k in sel[0] if k.startswith(("ok_", "rmsd_", "base_", "nitro", "lig_")))
        w = csv.DictWriter(fh, fieldnames=list(dict.fromkeys(c2)), extrasaction="ignore")
        w.writeheader(); w.writerows(sel)

    print(f"\ncriteria applied: {CRITERIA}")
    print("NOT applied (need his custom metrics/TMol): pocket volume, pp clash, TMol dG")
    for c in CRITERIA:
        print(f"  {c:12} {sum(bool(r['ok_' + c]) for r in ok):>4}/{len(ok)}")
    dist = np.bincount([r["n_pass"] for r in ok], minlength=len(CRITERIA) + 1)
    print("designs by number of criteria passed:", dict(enumerate(dist.tolist())))
    print(f"pass ALL applied criteria: {sum(r['pass_applied'] for r in ok)}")
    print(f"mode={a.mode} -> selected {len(sel)} designs  ->  {out / 'selection.csv'}")
    json.dump({"mode": a.mode, "k": a.k, "criteria": CRITERIA, "n_designs": len(ok),
               "n_selected": len(sel)}, open(out / "selection_meta.json", "w"))


if __name__ == "__main__":
    main()
