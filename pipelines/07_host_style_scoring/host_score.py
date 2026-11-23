"""Score candidate sequences the way the hosts' slide describes (2026-11-21).

Four components, each mapped to [0, 1] (1 = best):

  iPDE         Boltz-2 `complex_ipde` (protein-ligand PDE, contact-weighted;
               lower = ligand more confidently placed). Mean over seeds.
               score = clip((4.0 - ipde) / 3.5, 0, 1)  [1 at <=0.5, 0 at >=4.0;
               observed range over 425 of our refolds: 0.42-5.19, median 1.20]
  Geometry     ten distance/angle/dihedral measurements for the base (Asp O ->
               N3) and the H-bond donor (Gln NE2 -> N1), each scored
               exp(-0.5 z^2) against the crystal-calibrated ideal
               (configs/host_reference.yaml); score = mean of the ten. Mean over seeds.
  Burial       CA atoms within 12 A of the catalytic base CG, score = min(1,
               count / crystal mean 45.9). Mean over seeds.
  Consistency  reproducibility across seeds (needs >= 2): mean of
               exp(-ligand-pose RMSD / 1 A), exp(-sd(b_dist)/0.25 A),
               exp(-sd(d_dist)/0.25 A), exp(-sd(burial)/3 CA), exp(-sd(iPDE)/1 A),
               where the pose RMSD is the median pairwise ligand heavy-atom RMSD
               after CA superposition (the pose term counts double).

The hosts' formulas, weights and cut-offs are not public. Everything above is a
documented stand-in calibrated on the HG3-series crystals; component values are
always reported so the ranking can be re-weighted.

composite_abs   = mean of the available component scores
composite_rank  = mean percentile rank within the pool (robust to the fixed maps)

Hard flags (reported, never silently dropped): length > 260, non-standard
residues, catalytic residue identity wrong, ligand missing, a seed missing.

Usage: python -I host_score.py <candidates.csv> <seeds root> <out dir>
       candidates.csv columns: cand_id, sequence, base_pos, donor_pos, source
       seeds root: one dir per seed (seed_*) holding Boltz output; a single
       directory is treated as one seed.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import sys
import warnings
from itertools import combinations
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import yaml
from biotite.structure.io.pdbx import CIFFile, get_structure

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from host_geometry import (ANGULAR, BASE_METRICS, DIHEDRAL, DONOR_METRICS,  # noqa: E402
                           burial_count, circ_diff, geometry)
from kemp import geometry as g  # noqa: E402

REF = yaml.safe_load(open(ROOT / "configs" / "host_reference.yaml"))["metrics"]
AMAP = g.load_atom_map(ROOT / "configs" / "atom_map.yaml")
STD = set("ACDEFGHIKLMNPQRSTVWY")
CANDS: dict = {}
SEED_DIRS: list = []


def geom_scores(m: dict) -> dict:
    out = {}
    for k in BASE_METRICS + DONOR_METRICS:
        x = abs(m[k]) if k in DIHEDRAL else m[k]
        z = (x - REF[k]["mean"]) / REF[k]["sigma"]
        out["sc_" + k] = float(np.exp(-0.5 * z * z))   # "sc_" so the raw measurement m[k] is NOT overwritten
    out["geom_base"] = float(np.mean([out["sc_" + k] for k in BASE_METRICS]))
    out["geom_donor"] = float(np.mean([out["sc_" + k] for k in DONOR_METRICS]))
    out["geometry"] = float(np.mean([out["sc_" + k] for k in BASE_METRICS + DONOR_METRICS]))
    return out


def one_seed(cid: str, sdir: str) -> dict | None:
    c = CANDS[cid]
    hit = glob.glob(f"{sdir}/**/predictions/{cid}/{cid}_model_0.cif", recursive=True)
    if not hit:
        return None
    cif = Path(hit[0])
    s = get_structure(CIFFile.read(cif), model=1)
    r = {"cif": str(cif)}
    conf = cif.parent / f"confidence_{cid}_model_0.json"
    if conf.exists():
        j = json.loads(conf.read_text())
        r.update(ipde=j.get("complex_ipde"), plddt=j.get("complex_plddt"),
                 iptm=j.get("iptm"), ptm=j.get("ptm"))
    base, donor = int(c["base_pos"]), int(c["donor_pos"])
    lc, lr = g.iter_ligand_copies(s, AMAP)[0]
    m = geometry(s, chain="A", base_resid=base, donor_resid=donor, lig_chain=lc, lig_resid=lr)
    r.update({k: v for k, v in m.items()})
    r.update(geom_scores(m))
    b = burial_count(s, chain="A", base_resid=base)
    r.update(b); r["burial"] = min(1.0, b["burial_CG"] / REF["burial_CG"]["mean"])
    # pocket-restricted iPDE (secondary): PDE between ligand tokens and residues <10 A
    pde = cif.parent / f"pde_{cid}_model_0.npz"
    if pde.exists():
        P = np.load(pde)["pde"]
        prot = s[(~s.hetero) & (s.chain_id == "A") & (s.atom_name == "CA")]
        lig = s[s.res_name == "6NT"]
        nL = len(prot)
        d = np.linalg.norm(prot.coord[:, None, :] - lig.coord[None, :, :], axis=-1).min(1)
        near = d < 10.0
        if P.shape[0] == nL + len(lig) and near.any():
            r["ipde_pocket"] = float(0.5 * (P[:nL][near][:, nL:].mean() + P[nL:, :nL][:, near].mean()))
    # arrays for cross-seed pose comparison
    r["_ca"] = prot.coord if pde.exists() or True else None
    r["_lig"] = lig.coord
    tips = []
    for rid, names in ((base, ["CG", "OD1", "OD2"]), (donor, ["CD", "OE1", "NE2"])):
        for n in names:
            mm = (s.chain_id == "A") & (s.res_id == rid) & (s.atom_name == n) & (~s.hetero)
            tips.append(s.coord[mm][0])
    r["_tips"] = np.array(tips)
    return r


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    v, _, wt = np.linalg.svd((P - pc).T @ (Q - qc))
    if np.linalg.det(v @ wt) < 0:
        v[:, -1] *= -1
    return v @ wt, pc, qc


def score_candidate(cid: str) -> dict:
    c = CANDS[cid]
    row = {"cand_id": cid, "source": c.get("source", ""), "length": len(c["sequence"]),
           "base_pos": c["base_pos"], "donor_pos": c["donor_pos"], "flags": []}
    seq = c["sequence"]
    if len(seq) > 260:
        row["flags"].append("too-long")
    if not set(seq) <= STD:
        row["flags"].append("nonstandard-residue")
    if seq[int(c["base_pos"]) - 1] != "D":
        row["flags"].append("base-not-Asp")
    if seq[int(c["donor_pos"]) - 1] != "Q":
        row["flags"].append("donor-not-Gln")
    seeds = []
    for sd in SEED_DIRS:
        try:
            r = one_seed(cid, sd)
        except Exception as exc:  # noqa: BLE001 - record, never drop
            row["flags"].append(f"seed-error:{type(exc).__name__}")
            r = None
        if r is None:
            row["flags"].append("seed-missing")
        else:
            seeds.append(r)
    row["n_seeds"] = len(seeds)
    if not seeds:
        row["flags"].append("NO-PREDICTION"); row["flags"] = ";".join(row["flags"]); return row

    def mean(k):
        v = [s[k] for s in seeds if s.get(k) is not None]
        return float(np.mean(v)) if v else None
    for k in ("ipde", "ipde_pocket", "plddt", "iptm", "geometry", "geom_base", "geom_donor",
              "burial", "burial_CG", "burial_O", "b_dist", "d_dist"):
        row[k] = mean(k)
    row["s_ipde"] = (None if row["ipde"] is None else float(np.clip((4.0 - row["ipde"]) / 3.5, 0, 1)))
    row["s_geometry"], row["s_burial"] = row["geometry"], row["burial"]

    # ---- consistency across seeds --------------------------------------
    row["s_consistency"] = None
    if len(seeds) >= 2:
        pose = []
        for a, b in combinations(seeds, 2):
            R, pc, qc = kabsch(b["_ca"], a["_ca"])
            pose.append(float(np.sqrt((((b["_lig"] - pc) @ R + qc - a["_lig"]) ** 2).sum(1).mean())))
        sd = lambda k: float(np.std([s[k] for s in seeds if s.get(k) is not None]))  # noqa: E731
        parts = {"pose": np.exp(-np.median(pose) / 1.0),
                 "b_dist": np.exp(-sd("b_dist") / 0.25), "d_dist": np.exp(-sd("d_dist") / 0.25),
                 "burial": np.exp(-sd("burial_CG") / 3.0)}
        if all(s.get("ipde") is not None for s in seeds):
            parts["ipde"] = np.exp(-sd("ipde") / 1.0)
        row["lig_pose_rmsd_median"] = float(np.median(pose))
        w = {"pose": 2.0}   # ligand-pose reproducibility is the main preorganisation signal
        row["s_consistency"] = float(np.average(list(parts.values()),
                                                weights=[w.get(k, 1.0) for k in parts]))
        row["consistency_parts"] = json.dumps({k: round(float(v), 3) for k, v in parts.items()})
    comps = [row[k] for k in ("s_ipde", "s_geometry", "s_burial", "s_consistency") if row.get(k) is not None]
    row["composite_abs"] = float(np.mean(comps)) if comps else None
    row["n_components"] = len(comps)
    row["flags"] = ";".join(sorted(set(row["flags"])))
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates"); ap.add_argument("seeds_root"); ap.add_argument("out")
    ap.add_argument("--procs", type=int, default=4); ap.add_argument("--top", type=int, default=20)
    a = ap.parse_args()
    for r in csv.DictReader(open(a.candidates)):
        CANDS[r["cand_id"]] = r
    root = Path(a.seeds_root)
    sd = sorted(p for p in root.glob("seed_*") if p.is_dir())
    SEED_DIRS.extend(str(p) for p in (sd or [root]))
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    print(f"{len(CANDS)} candidates, {len(SEED_DIRS)} seed dir(s): {[Path(p).name for p in SEED_DIRS]}")
    with Pool(a.procs) as p:
        rows = p.map(score_candidate, list(CANDS), chunksize=4)

    scored = [r for r in rows if r.get("composite_abs") is not None]
    # percentile-rank composite (robust to the fixed score maps)
    for key in ("s_ipde", "s_geometry", "s_burial", "s_consistency"):
        vals = np.array([r[key] for r in scored if r.get(key) is not None])
        for r in scored:
            if r.get(key) is None:
                continue
            r["pr_" + key[2:]] = float((vals < r[key]).sum() + 0.5 * (vals == r[key]).sum()) / max(len(vals), 1)
    for r in scored:
        pr = [r[k] for k in ("pr_ipde", "pr_geometry", "pr_burial", "pr_consistency") if r.get(k) is not None]
        r["composite_rank"] = float(np.mean(pr))
    scored.sort(key=lambda r: -r["composite_abs"])
    cols = sorted({k for r in rows for k in r if not k.startswith("_")})
    first = ["cand_id", "source", "composite_abs", "composite_rank", "s_ipde", "s_geometry",
             "s_burial", "s_consistency", "n_seeds", "flags"]
    with open(out / "ranking.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=first + [c for c in cols if c not in first], extrasaction="ignore")
        w.writeheader(); w.writerows(scored + [r for r in rows if r.get("composite_abs") is None])
    f = lambda v: f"{v:6.2f}" if isinstance(v, float) else f"{'-':>6}"  # noqa: E731
    print(f"\n{'#':>3} {'candidate':24}{'comp':>6}{'rank':>6}{'iPDE':>6}{'geom':>6}{'bur':>6}{'cons':>6}  {'ipde_A':>6}{'b_d':>5}{'d_d':>5}{'CA12':>5}  flags")
    for i, r in enumerate(scored[: a.top], 1):
        print(f"{i:>3} {r['cand_id'][:24]:24}{f(r['composite_abs'])}{f(r['composite_rank'])}{f(r['s_ipde'])}"
              f"{f(r['s_geometry'])}{f(r['s_burial'])}{f(r.get('s_consistency'))}  {f(r.get('ipde'))}"
              f"{r['b_dist']:5.2f}{r['d_dist']:5.2f}{r['burial_CG']:5.0f}  {r['flags']}")
    print(f"\nscored {len(scored)}/{len(rows)}; wrote {out / 'ranking.csv'}")


if __name__ == "__main__":
    main()
