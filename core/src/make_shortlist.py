"""Assemble the shortlist of designs for external evaluation (e.g. PLACER).

Selection = union of
  * STRICT  : fold OK (TM>=0.85, scRMSD<2 A) + ligand buried + catalytic
              constellation <1.5 A + donor tip within 1.5 A of the design
  * CONTACTS: fold OK + buried + base->N3 <=3.4 A and donor->N1 <=3.5 A in the
              refold (real contacts, whether or not the shift criteria pass)

Tiers, so a consumer can cut the list however it likes:
  T1  strict AND real contacts AND substrate access >=1.4 A
  T2  real contacts AND access >=1.4 A          (superset of T1)
  T3  strict AND access >=1.4 A
  T4  everything else in the union

Writes (nothing is submitted anywhere):
  structures_cif/<design>.cif     Boltz-2 refold with 6NT (chain A protein, B ligand)
  sequences.fasta
  manifest.csv                    metrics, tiers, catalytic positions (1-indexed, chain A)
  boltz_yamls/<design>_<LIG>.yaml Boltz inputs for 6NT / H5J / 6VP refolds
                                  (NOT run; H5J/6VP CCD acceptance untested)

Usage: python -I make_shortlist.py <project root> <out dir>
"""

from __future__ import annotations

import csv
import glob
import json
import shutil
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from base_environment import RUNS  # noqa: E402
from kemp.record import read_manifest  # noqa: E402

LIGANDS = {"6NT": "transition-state analogue", "H5J": "substrate", "6VP": "product"}


def fl(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


def main(root: Path, out: Path) -> None:
    rows = list(csv.DictReader(open(root / "results_all_arms.csv")))
    acc = {(r["arm"], r["design_id"]): fl(r["bottleneck"])
           for r in csv.DictReader(open(root / "results_substrate_access.csv"))}

    sel = []
    for r in rows:
        tm, sc, sasa = fl(r["tm"]), fl(r["scrmsd"]), fl(r["sasa"])
        fold = tm is not None and tm >= 0.85 and sc is not None and sc < 2.0
        buried = sasa is not None and sasa < 0.15
        strict = r["strict"] == "True"
        rb, rd = fl(r["rf_base"]), fl(r["rf_donor"])
        contacts = (fold and buried and rb is not None and rd is not None
                    and rb <= 3.4 and rd <= 3.5)
        if not (strict or contacts):
            continue
        a = acc.get((r["arm"], r["design_id"]))
        opened = a is not None and a >= 1.4
        tier = ("T1" if strict and contacts and opened else
                "T2" if contacts and opened else
                "T3" if strict and opened else "T4")
        r.update(strict_b=strict, contacts_b=contacts, access=a, tier=tier)
        sel.append(r)
    sel.sort(key=lambda r: (r["tier"], -(fl(r["ipsae"]) or 0)))

    meta, paths, seqs = {}, {}, {}
    for arm, (s1, raw, _) in RUNS.items():
        recs = {x.design_id: x for x in read_manifest(root / s1 / "manifest.jsonl")}
        for name, m in json.loads((root / raw / "designs.json").read_text()).items():
            meta[(arm, name)] = m
            h = glob.glob(f"{root}/{raw}/**/predictions/{name}/{name}_model_0.cif",
                          recursive=True)
            paths[(arm, name)] = h[0] if h else None
            fa = root / raw / "armB_sequences.fasta"
            parent = m.get("parent_design", name)
            if name == parent and parent in recs:
                seqs[(arm, name)] = recs[parent].sequence
    # LigandMPNN sequences live in their yaml inputs
    for arm, (_, raw, _) in RUNS.items():
        for y in glob.glob(f"{root}/{raw}/**/boltz_in/*.yaml", recursive=True):
            name = Path(y).stem
            d = yaml.safe_load(open(y))
            seqs.setdefault((arm, name), d["sequences"][0]["protein"]["sequence"])

    for d in ("structures_cif", "boltz_yamls"):
        (out / d).mkdir(parents=True, exist_ok=True)

    man, fasta = [], []
    for r in sel:
        key = (r["arm"], r["design_id"])
        src, m = paths.get(key), meta[key]
        seq = seqs.get(key)
        if not src or not seq:
            print(f"  skipping {key}: missing structure or sequence")
            continue
        # Design ids are hash(cell, index, seed): the SAME id appears in the
        # full-motif and tip-atom runs for different proteins (32 of 32 shared,
        # 0 identical sequences). Prefix with the backbone family so names are
        # unique; the original id stays in the manifest.
        fam = "full" if r["arm"] in ("A_full_RFD4", "B_full_MPNN") else "tip"
        did = f"{fam}-{r['design_id']}"
        shutil.copy2(src, out / "structures_cif" / f"{did}.cif")
        fasta.append(f">{did} arm={r['arm']} tier={r['tier']} len={len(seq)}\n{seq}")
        for lig in LIGANDS:
            (out / "boltz_yamls" / f"{did}_{lig}.yaml").write_text(yaml.safe_dump({
                "version": 1,
                "sequences": [{"protein": {"id": "A", "sequence": seq, "msa": "empty"}},
                              {"ligand": {"id": "B", "ccd": lig}}]}, sort_keys=False))
        cat = m["catalytic"]
        man.append(dict(
            design_id=did, original_id=r["design_id"], backbone_family=fam, tier=r["tier"], arm=r["arm"], parent_backbone=r["parent"],
            length=len(seq), base_pos=cat["base"], donor_pos=cat["donor"],
            wall_pos=cat["wall"], protein_chain="A", ligand_chain="B", ligand="6NT",
            access_bottleneck_A=r["access"], ipsae_p2l=r["ipsae"], plddt=r["plddt"],
            scrmsd=r["scrmsd"], tm_to_design=r["tm"], rf_base_N3=r["rf_base"],
            rf_donor_N1=r["rf_donor"], rel_sasa=r["sasa"], triad_rmsd=r["triad"],
            donor_tip_shift=r["donor_shift"], base_polar_nbrs=r["n_polar"],
            base_cation_nbrs=r["n_cation"], strict=r["strict_b"],
            real_contacts=r["contacts_b"]))

    (out / "sequences.fasta").write_text("\n".join(fasta) + "\n")
    with open(out / "manifest.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(man[0]))
        w.writeheader()
        w.writerows(man)

    from collections import Counter
    print(f"shortlist: {len(man)} designs   tiers: {dict(sorted(Counter(m['tier'] for m in man).items()))}")
    print(f"           from {len({m['parent_backbone'] for m in man})} backbones, "
          f"arms: {dict(Counter(m['arm'] for m in man))}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
