"""Gate and run-folder layout for the His-as-base motifs.

Motif-parameterised replacement for the Asp-specific extract/gate path.
Motif residues are found with motif_locate.locate (ligand-anchored), which
also yields the motif RMSD -- so conditioning that was not honoured is caught.

SUCCESS = all of:
  1. motif retained: base and partner tips within 0.5 A, wall (if any) within 1.0 A
  2. ligand buried: rel_sasa <= 0.15, enclosure >= 0.88, >= 35 contacts
  3. no clash: closest protein-ligand heavy-atom pair >= 2.0 A
  4. imidazole tautomer defined: an H-bond ACCEPTOR within 3.3 A of the proton-
     bearing ring N (NE2 when ND1 is the acceptor)
  5. attacking N free: NO H-bond donor within 3.3 A of ND1
  6. length <= 260

Reported but not gated: cations within 6 A of the ring, substrate-access
bottleneck, ND1 -> N3 and partner -> N1 distances.

Usage: python -I his_gate.py <root> <raw rfd4 out dir> <run folder>
"""

from __future__ import annotations

import csv
import glob
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402
from motif_locate import locate  # noqa: E402
from substrate_access import access  # noqa: E402

RING = ["CB", "CG", "ND1", "CD2", "CE1", "NE2"]
MOTIFS = {
    "kemp_his_hg317": dict(
        cif="inputs/kemp_motif_his.cif", lig_chain="L",
        roles={"base": dict(resid=127, resname="HIS", tips=RING),
               "partner": dict(resid=50, resname="GLN", tips=["CD", "OE1", "NE2"]),
               "wall": dict(resid=44, resname="TRP", tips=["CD1", "NE1", "CE2"])},
        partner_atom="NE2", tol={"base": 0.5, "partner": 0.5, "wall": 1.0}),
    "kemp_his_gdm": dict(
        cif="inputs/kemp_theozyme_his_gdm.cif", lig_chain="B",
        roles={"base": dict(resid=1, resname="HIS", tips=RING),
               "partner": dict(resid=2, resname="SER", tips=["CB", "OG"])},
        partner_atom="OG", tol={"base": 0.5, "partner": 0.5}),
}
ACCEPTORS = {("*", "O"), ("*", "OXT"), ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"),
             ("GLU", "OE2"), ("ASN", "OD1"), ("GLN", "OE1"), ("SER", "OG"),
             ("THR", "OG1"), ("TYR", "OH")}
DONORS = {("*", "N"), ("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
          ("ASN", "ND2"), ("GLN", "NE2"), ("SER", "OG"), ("THR", "OG1"),
          ("TYR", "OH"), ("TRP", "NE1")}
CATION = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2")}


def _in(table, rn, an):
    return (rn, an) in table or ("*", an) in table and not (an == "N" and rn == "PRO")


def count_near(prot, centre, own_resid, table, r):
    n = 0
    for rn, an, rid, xyz in zip(prot.res_name, prot.atom_name, prot.res_id, prot.coord):
        if rid == own_resid or not _in(table, rn, an):
            continue
        if np.linalg.norm(xyz - centre) <= r:
            n += 1
    return n


def analyse(cif: Path, key: str, amap) -> dict:
    cfg = MOTIFS[key]
    motif_in = get_structure(CIFFile.read(ROOT / cfg["cif"]), model=1)
    s = get_structure(CIFFile.read(cif), model=1)
    loc = locate(s, motif_in, cfg["roles"], tol=8.0)       # report, gate below
    rec = {f"rmsd_{k}": v["rmsd"] for k, v in loc.items()}
    brid, prid = loc["base"]["resid"], loc["partner"]["resid"]
    prot = s[(~s.hetero) & (s.chain_id == "A")]
    lc = cfg["lig_chain"]
    lig = s[(s.res_name == "6NT")]
    n3 = lig.coord[lig.atom_name == "N3"][0]
    n1 = lig.coord[lig.atom_name == "N1"][0]
    xyz = lambda rid, an: prot.coord[(prot.res_id == rid) & (prot.atom_name == an)][0]  # noqa: E731
    rec.update(base_resid=brid, partner_resid=prid, length=len(set(prot.res_id)))
    rec["nd1_N3"] = float(np.linalg.norm(xyz(brid, "ND1") - n3))
    rec["ne2_N3"] = float(np.linalg.norm(xyz(brid, "NE2") - n3))
    rec["partner_N1"] = float(np.linalg.norm(xyz(prid, cfg["partner_atom"]) - n1))
    # tautomer lock: acceptor near NE2; attacking ND1 must be free of donors
    rec["ne2_acceptors"] = count_near(prot, xyz(brid, "NE2"), brid, ACCEPTORS, 3.3)
    rec["nd1_donors"] = count_near(prot, xyz(brid, "ND1"), brid, DONORS, 3.3)
    ring_c = np.array([xyz(brid, a) for a in ("ND1", "NE2")]).mean(0)
    rec["ring_cations_6A"] = count_near(prot, ring_c, brid, CATION, 6.0)
    rec.update({f"bur_{k}": v for k, v in g.ligand_burial(
        s, amap, lig_chain=lc, lig_resid=int(lig.res_id[0])).items()})
    rec["min_contact"] = g.min_protein_ligand_contact(
        s, amap, lig_chain=lc, lig_resid=int(lig.res_id[0]))
    a = access(s)
    rec["access"] = a["bottleneck"] if a else None
    return rec


def verdict(r: dict, key: str) -> list[str]:
    cfg, why = MOTIFS[key], []
    for role, tol in cfg["tol"].items():
        if r[f"rmsd_{role}"] > tol:
            why.append(f"motif-{role}-moved({r[f'rmsd_{role}']:.2f}A)")
    if r["bur_rel_sasa"] > 0.15 or r["bur_enclosure"] < 0.88 or r["bur_n_contacts"] < 35:
        why.append("ligand-not-buried")
    if r["min_contact"] < 2.0:
        why.append("clash")
    if r["ne2_acceptors"] < 1:
        why.append("no-tautomer-lock")
    if r["nd1_donors"] > 0:
        why.append("ND1-blocked-by-donor")
    if r["length"] > 260:
        why.append("too-long")
    return why


def main(root: Path, raw: Path, run: Path) -> None:
    global ROOT
    ROOT = root
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    for d in ("attempts", "success"):
        if (run / d).exists():
            shutil.rmtree(run / d)
        (run / d).mkdir(parents=True, exist_ok=True)
    rows = []
    for key in MOTIFS:
        for f in sorted(glob.glob(f"{raw}/inference_samples/{key}_*/*_generated.cif")):
            tag = Path(f).stem.replace("_generated", "")
            try:
                r = analyse(Path(f), key, amap)
                why = verdict(r, key)
                v = "SUCCESS" if not why else "FAIL-" + why[0]
            except Exception as exc:  # noqa: BLE001 - record, never drop
                r, why, v = {}, [f"{type(exc).__name__}: {exc}"], "ERROR"
            r.update(motif=key, design=tag, verdict=v, reasons=";".join(why), cif=f)
            rows.append(r)
            shutil.copy2(f, run / "attempts" / f"{tag}_{v}.cif")
    succ = sorted([r for r in rows if r["verdict"] == "SUCCESS"],
                  key=lambda r: (r["bur_rel_sasa"], r["ring_cations_6A"]))
    for i, r in enumerate(succ, 1):
        shutil.copy2(r["cif"], run / "success" / f"rank{i:02d}_{r['design']}.cif")
    cols = sorted({k for r in rows for k in r})
    first = ["motif", "design", "verdict", "reasons"]
    with open(run / "summary.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=first + [c for c in cols if c not in first])
        w.writeheader(); w.writerows(rows)
    for key in MOTIFS:
        s = [r for r in rows if r["motif"] == key]
        n_ok = sum(r["verdict"] == "SUCCESS" for r in s)
        print(f"{key}: {n_ok}/{len(s)} success")
        from collections import Counter
        print("   failures:", dict(Counter(r["verdict"] for r in s if r["verdict"] != "SUCCESS")))
    print(f"run folder: {run}")


ROOT = Path(".")
if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve(), Path(sys.argv[2]), Path(sys.argv[3]))
