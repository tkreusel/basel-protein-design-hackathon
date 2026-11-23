"""How do the Boltz poses of the SAME protein differ between ligands?  (TSA 6NT vs substrate H5J vs product 6VP)

  python -m maxi_eval.arm_poses --inputs work/inputs --out work/arm_poses.csv

For each protein with all three arms (same model index) the substrate and product complexes are superposed on the TSA complex (C-alpha fit) and:
  arm.<x>_vs_tsa.ca_rmsd           global backbone difference
  arm.<x>_vs_tsa.pocket_bb_rmsd    backbone of the TSA pocket (6 A shell)       -> does the pocket re-shape for a different ligand?
  arm.<x>_vs_tsa.pocket_allatom_rmsd
  arm.<x>_vs_tsa.base_O_shift      catalytic carboxylate oxygens (A)
  arm.<x>_vs_tsa.donor_tip_shift   oxyanion-donor side-chain tip (A)
  arm.<x>_vs_tsa.lig_rmsd          ligand pose (atoms mapped through ligands.CORRESPONDENCE_TO_TSA) after the fit
  arm.<x>_vs_tsa.lig_centroid_shift
  arm.prod_vs_sub.lig_rmsd         product pose vs substrate pose (both mapped to the TSA frame)
x in {sub, prod}.  Small values = the pocket and ligand position are insensitive to which species is bound (a pre-organised site).
"""
from __future__ import annotations
import argparse, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from . import geometry as G, ligands as L, structures as S
from .apo_holo import transplant_ligand, _rmsd_atoms


def _mapped(lig, code, tsa_xyz):
    xyz = {str(n): c for n, c in zip(lig.atom_name, lig.coord)}
    corr = L._best_corr(L.CORRESPONDENCE_TO_TSA[code], xyz, tsa_xyz)
    return {corr[k]: xyz[k] for k in corr}


def _shift(Pa, Pb, rid, names, R, t):
    out = []
    for n in names:
        a, b = Pa.xyz(rid, n), Pb.xyz(rid, n)
        if a is not None and b is not None:
            out.append(np.linalg.norm((b @ R.T + t) - a))
    return out


def compare(tsa, other, code, site):
    out = {}
    lig_in_tsa, ca, R, t = transplant_ligand(tsa.protein, other.protein, other.ligand)
    out["ca_rmsd"] = ca
    pocket = set(site["pocket_resids"])
    out["pocket_bb_rmsd"] = _rmsd_atoms(tsa.protein, other.protein, pocket, R, t, True)
    out["pocket_allatom_rmsd"] = _rmsd_atoms(tsa.protein, other.protein, pocket, R, t, False)
    Pt, Po = G.Atoms.of(tsa.protein), G.Atoms.of(other.protein)
    b = site.get("base")
    if b:
        s1 = _shift(Pt, Po, b["resid"], [b["O"], b["O_other"]], R, t)
        s2 = _shift(Pt, Po, b["resid"], [b["O_other"], b["O"]], R, t)   # carboxylate O's are interchangeable
        if s1 and s2:
            out["base_O_shift"] = float(min(np.mean(s1), np.mean(s2)))
    sc = next((d for d in site["donors"] if not d["backbone"]), None)
    if sc:
        s = _shift(Pt, Po, sc["resid"], [sc["atom"]], R, t)
        if s:
            out["donor_tip_shift"] = float(s[0])
    tsa_xyz = {str(n): c for n, c in zip(tsa.ligand.atom_name, tsa.ligand.coord)}
    mp = _mapped(lig_in_tsa, code, tsa_xyz)
    keys = [k for k in mp if k in tsa_xyz]
    out["lig_rmsd"] = float(np.sqrt(np.mean([np.sum((mp[k] - tsa_xyz[k]) ** 2) for k in keys])))
    out["lig_centroid_shift"] = float(np.linalg.norm(np.mean([mp[k] for k in keys], 0) - np.mean([tsa_xyz[k] for k in keys], 0)))
    return out, mp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs")
    ap.add_argument("--out", default="work/arm_poses.csv")
    a = ap.parse_args()
    T = L.build_templates()
    man = json.loads((Path(a.inputs) / "manifest.json").read_text())
    rows = []
    for entry, e in man["entries"].items():
        if not entry.endswith("_boltz"):
            continue
        for r in [x for x in e["replicates"] if x["arm"] == "tsa"]:
            k = r.get("boltz_model")
            pick = {arm: next((x for x in e["replicates"] if x["arm"] == arm and x.get("boltz_model") == k), None) for arm in ("sub", "prod")}
            if not all(pick.values()):
                continue
            try:
                tsa = S.extract_replicates(Path(r["pdb"]), entry, "tsa", "boltz", T)[0]
                site = G.detect_site(tsa, e.get("catalytic_overrides"))
                row = {"rep_id": r["rep_id"], "entry": entry}
                mapped = {}
                for arm, code in (("sub", "H5J"), ("prod", "6VP")):
                    other = S.extract_replicates(Path(pick[arm]["pdb"]), entry, arm, "boltz", T)[0]
                    d, mp = compare(tsa, other, code, site)
                    mapped[arm] = mp
                    row.update({f"arm.{arm}_vs_tsa.{k_}": v for k_, v in d.items()})
                common = [k_ for k_ in mapped["sub"] if k_ in mapped["prod"]]
                row["arm.prod_vs_sub.lig_rmsd"] = float(np.sqrt(np.mean([np.sum((mapped["sub"][k_] - mapped["prod"][k_]) ** 2) for k_ in common])))
                rows.append(row); print("arm poses", r["rep_id"])
            except Exception as exc:  # fail loud, continue
                rows.append({"rep_id": r["rep_id"], "entry": entry, "arm.error": f"{type(exc).__name__}: {exc}"}); print("FAILED", r["rep_id"], exc)
    pd.DataFrame(rows).to_csv(a.out, index=False); print("->", a.out)

if __name__ == "__main__":
    main()
