"""Preorganisation from paired Boltz predictions: apo (no ligand) vs TSA-bound (6NT) of the same protein.

  python -m maxi_eval.apo_holo --inputs work/inputs --out work/apo_holo.csv

For every entry that has both an `apo` and a `tsa` Boltz replicate (same sequence, same numbering):
  1. superpose apo onto holo on all matching C-alpha atoms
  2. `ah.ca_rmsd`                      global backbone change on binding
  3. `ah.pocket_bb_rmsd/_allatom_rmsd`  residues within 6 A of the ligand (backbone / all heavy atoms): how much the pocket MOVES
  4. `ah.base_O_shift`, `ah.donor_tip_shift`, `ah.wall_centroid_shift`   displacement of the catalytic groups (A)
  5. `ah.d_chi_base`, `ah.d_chi_donor`   circular change of side-chain chi angles
  6. the holo TSA is TRANSPLANTED into the apo pocket (same rigid transform) and scored with the usual geometry:
     `ah.apo.*` (d_O_X, angles, constraint passes, `pass_core`) and `ah.apo.contact_min` (clash if < ~2.2 A).
A small `ah.*` shift plus a passing transplanted geometry = a pocket already formed before the ligand arrives (preorganised).
LIMIT: both are Boltz predictions; for native-sequence proteins Boltz may have memorised the holo structure, so apo ~ holo can be an artefact.
"""
from __future__ import annotations
import argparse, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from . import cst as C, geometry as G, ligands as L, structures as S

CORE = ["geo.d_O_X", "geo.d_O_H", "geo.ang_O_H_X", "geo.ang_C_O_H"]


def transplant_ligand(apo_prot, holo_prot, holo_lig):
    """Rigid transform (CA fit holo->apo over shared residues) applied to the ligand. Returns (ligand in apo frame, CA rmsd, R, t)."""
    a = {int(r): c for r, n, c in zip(apo_prot.res_id, apo_prot.atom_name, apo_prot.coord) if n == "CA"}
    h = {int(r): c for r, n, c in zip(holo_prot.res_id, holo_prot.atom_name, holo_prot.coord) if n == "CA"}
    ids = sorted(set(a) & set(h))
    if len(ids) < 20:
        raise RuntimeError(f"only {len(ids)} shared residues between apo and holo")
    src, dst = np.array([h[i] for i in ids]), np.array([a[i] for i in ids])
    R, t = L._kabsch(src, dst)
    lig = holo_lig.copy()
    lig.coord = (holo_lig.coord @ R.T + t).astype(np.float32)
    return lig, float(np.sqrt(((src @ R.T + t - dst) ** 2).sum(1).mean())), R, t


def _rmsd_atoms(apo_prot, holo_prot, resids, R, t, backbone_only=False):
    A = {(int(r), str(n)): c for r, n, c in zip(apo_prot.res_id, apo_prot.atom_name, apo_prot.coord)}
    d2 = []
    for r, n, c in zip(holo_prot.res_id, holo_prot.atom_name, holo_prot.coord):
        if int(r) in resids and (not backbone_only or str(n) in ("N", "CA", "C", "O")):
            k = (int(r), str(n))
            if k in A:
                d2.append(((c @ R.T + t) - A[k]) @ ((c @ R.T + t) - A[k]))
    return float(np.sqrt(np.mean(d2))) if d2 else float("nan")


def compare(apo_rep, holo_rep, site, constraints) -> dict:
    out = {}
    lig_in_apo, ca, R, t = transplant_ligand(apo_rep.protein, holo_rep.protein, holo_rep.ligand)
    out["ah.ca_rmsd"] = ca
    pocket = set(site["pocket_resids"])
    out["ah.pocket_bb_rmsd"] = _rmsd_atoms(apo_rep.protein, holo_rep.protein, pocket, R, t, True)
    out["ah.pocket_allatom_rmsd"] = _rmsd_atoms(apo_rep.protein, holo_rep.protein, pocket, R, t, False)
    Ah, Aa = G.Atoms.of(holo_rep.protein), G.Atoms.of(apo_rep.protein)

    def tx(x):
        return None if x is None else x @ R.T + t
    b = site.get("base")
    if b:
        o1, o2 = b["O"], b["O_other"]
        h1, h2, a1, a2 = tx(Ah.xyz(b["resid"], o1)), tx(Ah.xyz(b["resid"], o2)), Aa.xyz(b["resid"], o1), Aa.xyz(b["resid"], o2)
        if all(x is not None for x in (h1, h2, a1, a2)):
            out["ah.base_O_shift"] = float(min((np.linalg.norm(h1 - a1) + np.linalg.norm(h2 - a2)) / 2, (np.linalg.norm(h1 - a2) + np.linalg.norm(h2 - a1)) / 2))
        for k, v in G.chi_angles(apo_rep, b["resid"]).items():
            hv = G.chi_angles(holo_rep, b["resid"]).get(k)
            if hv is not None:
                out[f"ah.d_chi_base_{k}"] = abs(G.wrap180(v - hv))
    sc = next((d for d in site["donors"] if not d["backbone"]), None)
    if sc:
        h, a = tx(Ah.xyz(sc["resid"], sc["atom"])), Aa.xyz(sc["resid"], sc["atom"])
        if h is not None and a is not None:
            out["ah.donor_tip_shift"] = float(np.linalg.norm(h - a))
        for k, v in G.chi_angles(apo_rep, sc["resid"]).items():
            hv = G.chi_angles(holo_rep, sc["resid"]).get(k)
            if hv is not None:
                out[f"ah.d_chi_donor_{k}"] = abs(G.wrap180(v - hv))
    w = site.get("wall")
    if w:
        ph = [tx(Ah.xyz(w["resid"], x)) for x in w["ring"] if Ah.xyz(w["resid"], x) is not None]
        pa = [Aa.xyz(w["resid"], x) for x in w["ring"] if Aa.xyz(w["resid"], x) is not None]
        if len(ph) >= 5 and len(pa) >= 5:
            out["ah.wall_centroid_shift"] = float(np.linalg.norm(np.mean(ph, 0) - np.mean(pa, 0)))
    # holo ligand transplanted into the apo pocket, scored with the standard geometry
    tr = S.Replicate("apo+TSA", apo_rep.entry, "tsa", "boltz", apo_rep.protein, lig_in_apo)
    m, _ = G.all_geometry(tr, site=site)
    ev = C.evaluate(m, constraints)
    for k in ("geo.d_O_X", "geo.d_O_H", "geo.ang_O_H_X", "geo.ang_C_O_H", "geo.tors_syn", "geo.contact_min", "geo.lig_rel_sasa", "geo.enclosure"):
        if k in m:
            out[f"ah.apo.{k[4:]}"] = m[k]
    core = [c.short for c in constraints if c.key in CORE]
    out["ah.apo.pass_core"] = float(all(ev.get(f"cst.{k}.pass", 0.0) == 1.0 for k in core)) if core and "geo.d_O_X" in m else float("nan")
    out["ah.apo.cst_n_pass"] = ev.get("cst.n_pass", float("nan"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs")
    ap.add_argument("--out", default="work/apo_holo.csv")
    ap.add_argument("--reference", default="config/theozyme_reference.yaml")
    ap.add_argument("--group", default="hg3_lineage")
    a = ap.parse_args()
    T = L.build_templates()
    cons = C.load_reference(Path(a.reference), a.group)
    man = json.loads((Path(a.inputs) / "manifest.json").read_text())
    rows = []
    for entry, e in man["entries"].items():
        apos = [r for r in e["replicates"] if r["arm"] == "apo"]
        if not apos:
            continue
        for r in [x for x in e["replicates"] if x["arm"] == "tsa"]:
            try:
                pick = next((a_ for a_ in apos if a_.get("boltz_model") == r.get("boltz_model")), apos[0])
                apo_rep = S.extract_apo(Path(pick["pdb"]), entry)
                holo = S.extract_replicates(Path(r["pdb"]), entry, "tsa", r["source"], T)[0]
                site = G.detect_site(holo, e.get("catalytic_overrides"))
                rows.append({"rep_id": r["rep_id"], "entry": entry, **compare(apo_rep, holo, site, cons)})
                print("apo-holo", r["rep_id"])
            except Exception as exc:  # fail loud, continue
                rows.append({"rep_id": r["rep_id"], "entry": entry, "ah.error": f"{type(exc).__name__}: {exc}"})
                print("FAILED", r["rep_id"], exc)
    pd.DataFrame(rows).to_csv(a.out, index=False)
    print("->", a.out)

if __name__ == "__main__":
    main()
