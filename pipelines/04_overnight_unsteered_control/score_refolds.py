"""Stage 5: score every Boltz refold against its RFD4 design.

Per prediction (seq_id = <design>__rfd4 | <design>__mpnnK):
  fold          TM-score and full-length CA scRMSD vs the designed backbone
  site          catalytic triad RMSD (9 tip atoms, after superposition), Gln
                donor tip shift, ligand centroid shift
  contacts      base->N3, donor->N1, wall stacking, ligand burial
  base env      polar atoms <=4.5 A, Lys/Arg <=6 A, base SASA, nitro H-bonds
                (Maxi's definitions, our validated reimplementation)
  access        substrate-access bottleneck radius (ligand removed)
  confidence    Boltz pLDDT / pTM / ipTM, and ipSAE protein->ligand if the
                ipsae TSVs are provided

Run on CPU (head node or alongside the GPU job); multiprocessing over
predictions.

Usage: python -I score_refolds.py <selection.csv> <sequences.csv> <boltz_out root>
                                  <out csv> [--ipsae-dir DIR] [--procs 4]
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import subprocess
import sys
import tempfile
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kemp import geometry as g  # noqa: E402
from steering_criteria import measure  # noqa: E402
from substrate_access import access  # noqa: E402

AMAP = g.load_atom_map(ROOT / "configs" / "atom_map.yaml")
USALIGN = ROOT / "bin" / "USalign"
TIPS = {"base": ("ASP", ["CG", "OD1", "OD2"]), "partner": ("GLN", ["CD", "OE1", "NE2"]),
        "wall": ("TRP", ["CD1", "NE1", "CE2"])}
SEL: dict = {}


def _ca(s):
    return s[(~s.hetero) & (s.chain_id == "A") & (s.atom_name == "CA")]


def _tm(design, refold):
    if not USALIGN.exists():
        return float("nan")
    with tempfile.TemporaryDirectory() as t:
        paths = []
        for n, st in (("d", design), ("r", refold)):
            f = PDBFile(); f.set_structure(st[~st.hetero]); p = Path(t) / f"{n}.pdb"; f.write(p)
            paths.append(str(p))
        o = subprocess.run([str(USALIGN), paths[1], paths[0]], capture_output=True, text=True).stdout
    for line in o.splitlines():
        if line.startswith("TM-score="):
            return float(line.split("=")[1].split()[0])
    return float("nan")


def _tips(s, rid, names):
    pts = []
    for n in names:
        m = (~s.hetero) & (s.chain_id == "A") & (s.res_id == rid) & (s.atom_name == n)
        if not m.any():
            return None
        pts.append(s.coord[m][0])
    return np.array(pts)


def score(job):
    seq_id, cif, conf = job
    design = seq_id.split("__")[0]
    sel = SEL[design]
    r = {"seq_id": seq_id, "design": design, "source": seq_id.split("__")[1]}
    try:
        rs = get_structure(CIFFile.read(cif), model=1)
        ds = get_structure(CIFFile.read(sel["cif"]), model=1)
        if conf and Path(conf).exists():
            c = json.loads(Path(conf).read_text())
            r.update(plddt=c.get("complex_plddt"), ptm=c.get("ptm"), iptm=c.get("iptm"))
        dca, rca = _ca(ds), _ca(rs)
        if len(dca) != len(rca):
            r["error"] = "length mismatch"; return r
        r["scrmsd"] = g.kabsch_rmsd(rca.coord, dca.coord)
        r["tm"] = _tm(ds, rs)
        ids = {k: int(float(sel[f"{k}_resid"])) for k in TIPS}
        # backbone superposition refold -> design
        pc, qc = rca.coord.mean(0), dca.coord.mean(0)
        v, _, wt = np.linalg.svd((rca.coord - pc).T @ (dca.coord - qc))
        if np.linalg.det(v @ wt) < 0:
            v[:, -1] *= -1
        R = v @ wt
        tf = lambda x: (x - pc) @ R + qc  # noqa: E731
        d_t, r_t = [], []
        for k, (rn, names) in TIPS.items():
            a, b = _tips(ds, ids[k], names), _tips(rs, ids[k], names)
            if a is None or b is None:
                raise ValueError(f"missing {k} tip atoms")
            d_t.append(a); r_t.append(b)
        r["triad_rmsd"] = g.kabsch_rmsd(np.vstack(r_t), np.vstack(d_t))
        ne2d = _tips(ds, ids["partner"], ["NE2"])[0]
        ne2r = _tips(rs, ids["partner"], ["NE2"])[0]
        r["donor_tip_shift"] = float(np.linalg.norm(tf(ne2r) - ne2d))
        dl, rl = ds[ds.res_name == "6NT"], rs[rs.res_name == "6NT"]
        r["ligand_shift"] = float(np.linalg.norm(tf(rl.coord.mean(0)) - dl.coord.mean(0)))
        lc, lr = g.iter_ligand_copies(rs, AMAP)[0]
        kw = dict(lig_chain=lc, lig_resid=lr, chain="A")
        r["rf_base_N3"] = g.base_to_abstraction(rs, AMAP, base_resid=ids["base"], **kw)
        r["rf_donor_N1"] = g.donor_to_oxyanion(rs, AMAP, donor_resid=ids["partner"], **kw)
        r["rf_wall"] = g.wall_stacking(rs, AMAP, wall_resid=ids["wall"], **kw)
        r.update({f"rf_{k}": v for k, v in g.ligand_burial(rs, AMAP, lig_chain=lc,
                                                            lig_resid=lr).items()})
        m = measure(rs)
        r.update({f"env_{k}": m[k] for k in ("base_sasa_rel", "lig_sasa_frac",
                  "base_sc_polar_4p5", "base_cation_6", "nitro_hb_far", "nitro_hb_near_base",
                  "lig_clash")})
        a = access(rs)
        r["access"] = a["bottleneck"] if a else None
    except Exception as exc:  # noqa: BLE001 - record, never drop
        r["error"] = f"{type(exc).__name__}: {exc}"
    return r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("selection"); ap.add_argument("sequences"); ap.add_argument("boltz_root")
    ap.add_argument("out"); ap.add_argument("--ipsae-dir"); ap.add_argument("--procs", type=int, default=4)
    a = ap.parse_args()
    for r in csv.DictReader(open(a.selection)):
        SEL[r["design"]] = r
    jobs = []
    for r in csv.DictReader(open(a.sequences)):
        sid = r["seq_id"]
        h = glob.glob(f"{a.boltz_root}/**/predictions/{sid}/{sid}_model_0.cif", recursive=True)
        if h:
            conf = glob.glob(f"{a.boltz_root}/**/predictions/{sid}/confidence_{sid}_model_0.json",
                             recursive=True)
            jobs.append((sid, h[0], conf[0] if conf else None))
    print(f"{len(jobs)} predictions to score")
    with Pool(a.procs) as p:
        rows = p.map(score, jobs, chunksize=4)

    if a.ipsae_dir:
        ips = {}
        for f in glob.glob(f"{a.ipsae_dir}/*.tsv"):
            for t in csv.DictReader(open(f), delimiter="\t"):
                if t.get("type1") == "prot" and t.get("type2") == "lig":
                    ips[Path(t["model"]).parent.name] = float(t["ipsae12"])
        for r in rows:
            r["ipsae_p2l"] = ips.get(r["seq_id"])

    def ok(r, k, thr, lt=True):
        v = r.get(k)
        return v is not None and v == v and (v < thr if lt else v >= thr)
    for r in rows:
        r["fold_ok"] = ok(r, "tm", 0.85, False) and ok(r, "scrmsd", 2.0)
        r["buried_ok"] = ok(r, "rf_rel_sasa", 0.15)
        r["const_ok"] = ok(r, "triad_rmsd", 1.5)
        r["donor_ok"] = ok(r, "donor_tip_shift", 1.5)
        r["contacts_ok"] = ok(r, "rf_base_N3", 3.4) and ok(r, "rf_donor_N1", 3.5)
        r["env_ok"] = ok(r, "env_base_sc_polar_4p5", 1.01) and r.get("env_base_cation_6") == 0
        r["strict"] = r["fold_ok"] and r["buried_ok"] and r["const_ok"] and r["donor_ok"]
        r["strict_env"] = r["strict"] and r["env_ok"]
        r["strict_open"] = r["strict"] and ok(r, "access", 1.4, False)
    cols = sorted({k for r in rows for k in r})
    first = ["seq_id", "design", "source"]
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=first + [c for c in cols if c not in first])
        w.writeheader(); w.writerows(rows)
    good = [r for r in rows if "error" not in r]
    print(f"scored {len(good)} ok, {len(rows) - len(good)} errors -> {a.out}")
    print(f"\n{'source':8}{'n':>6}{'fold':>6}{'+bur':>6}{'+const':>7}{'+donor':>7}{'contacts':>9}{'env':>5}{'strict':>7}{'+env':>6}{'+open':>6}")
    for src in ("rfd4", "mpnn"):
        s = [r for r in good if r["source"].startswith(src)]
        f1 = [r for r in s if r["fold_ok"]]; b1 = [r for r in f1 if r["buried_ok"]]
        c1 = [r for r in b1 if r["const_ok"]]; d1 = [r for r in c1 if r["donor_ok"]]
        print(f"{src:8}{len(s):>6}{len(f1):>6}{len(b1):>6}{len(c1):>7}{len(d1):>7}"
              f"{sum(r['contacts_ok'] for r in s):>9}{sum(r['env_ok'] for r in s):>5}"
              f"{sum(r['strict'] for r in s):>7}{sum(r['strict_env'] for r in s):>6}"
              f"{sum(r['strict_open'] for r in s):>6}")


if __name__ == "__main__":
    main()
