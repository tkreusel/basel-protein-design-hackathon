"""Extract designs from RFD4 output and prepare the downstream hand-offs.

RFD4-Proteina co-designs sequence with structure, so each output CIF already
carries a sequence -- that is design **arm A**. This script also writes a PDB
per design so LigandMPNN can re-design the same backbone as **arm B**.

The hard part is locating the catalytic residues. Our motif is *unindexed*, so
the model chooses where the base, donor and wall land in sequence; their output
positions are not 44/50/127. We recover them from RFD4's provenance columns
rather than by guessing:

  condition_sequence_1_residue        the conditioned residue name, or <M>
                                      where the model designed the identity
  annotation_src-component_1_atom     which input component an atom came from
  condition_CoordinateGroup_1_atom    our theozyme rigid-body group

Fail loud (plan 0.2): a design whose catalytic residues cannot be located, or
whose conditioned identity does not match what we asked for, is recorded with
status "failed" and a reason. It is not silently dropped and not repaired.

Usage: python -I extract_designs.py <project root> <rfd4 out_dir> <stage dir>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp import geometry as g  # noqa: E402
from kemp.record import DesignRecord, append_manifest, design_id  # noqa: E402

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# What we conditioned, by role. Used to VERIFY the output, not to locate it.
EXPECTED = {"base": "ASP", "donor": "GLN", "wall": "TRP"}
SRC_BY_ROLE = {"base": "A127", "donor": "A50", "wall": "A44"}

MASK_TOKENS = {"<M>", ".", "?", "", "nan", "UNK"}


def _col(atom_site, name):
    if name not in atom_site:
        return None
    return atom_site[name].as_array(str)


def parse_design(cif_path: Path, amap: g.AtomMap) -> dict:
    """Pull sequence, catalytic positions and ligand identity out of one CIF."""
    cif = CIFFile.read(cif_path)
    atom_site = cif.block["atom_site"]

    asym = _col(atom_site, "label_asym_id")
    seq_id = _col(atom_site, "label_seq_id")
    comp = _col(atom_site, "label_comp_id")
    cond_seq = _col(atom_site, "condition_sequence_1_residue")
    src = _col(atom_site, "annotation_src-component_1_atom")

    if asym is None or comp is None:
        raise ValueError(f"{cif_path.name}: atom_site lacks label_asym_id/comp_id")

    # --- the protein chain: the one holding standard amino acids -----------
    prot_mask = np.isin(comp, list(THREE_TO_ONE))
    if not prot_mask.any():
        raise ValueError(f"{cif_path.name}: no standard amino acids found")
    prot_chains = sorted(set(asym[prot_mask]))
    if len(prot_chains) != 1:
        raise ValueError(
            f"{cif_path.name}: expected 1 protein chain, found {prot_chains}"
        )
    pchain = prot_chains[0]

    # --- sequence, in residue order ---------------------------------------
    on_chain = (asym == pchain) & prot_mask
    order, seen = [], {}
    for i in np.flatnonzero(on_chain):
        rid = int(seq_id[i])
        if rid not in seen:
            seen[rid] = comp[i]
            order.append(rid)
    order.sort()
    sequence = "".join(THREE_TO_ONE[seen[r]] for r in order)
    # Output res_id is a fresh contiguous per-chain counter starting at 1, so
    # position in `order` maps to a 1-indexed sequence position.
    pos_of_resid = {rid: n + 1 for n, rid in enumerate(order)}

    # --- catalytic residues, from provenance ------------------------------
    catalytic, cat_detail = {}, {}
    for role, want_src in SRC_BY_ROLE.items():
        hit_resids = set()
        if src is not None:
            for i in np.flatnonzero(on_chain):
                if src[i] == want_src:
                    hit_resids.add(int(seq_id[i]))
        if not hit_resids and cond_seq is not None:
            # Fallback: conditioned residues are the ones not masked.
            for i in np.flatnonzero(on_chain):
                if cond_seq[i] not in MASK_TOKENS and comp[i] == EXPECTED[role]:
                    hit_resids.add(int(seq_id[i]))
        if len(hit_resids) != 1:
            raise ValueError(
                f"{cif_path.name}: could not uniquely locate {role} "
                f"(src={want_src}); candidates={sorted(hit_resids)}"
            )
        rid = hit_resids.pop()
        got = seen[rid]
        if got != EXPECTED[role]:
            raise ValueError(
                f"{cif_path.name}: {role} at res {rid} is {got}, "
                f"expected {EXPECTED[role]} -- conditioning was not honoured"
            )
        catalytic[role] = pos_of_resid[rid]
        cat_detail[role] = rid

    # --- ligand ------------------------------------------------------------
    lig_mask = comp == amap.ligand_ccd_id
    if not lig_mask.any():
        raise ValueError(f"{cif_path.name}: ligand {amap.ligand_ccd_id} absent")
    lig_chain = sorted(set(asym[lig_mask]))[0]
    lig_resid = int(sorted(set(seq_id[lig_mask]))[0])

    return dict(
        sequence=sequence,
        length=len(sequence),
        protein_chain=pchain,
        catalytic=catalytic,
        catalytic_resid=cat_detail,
        lig_chain=lig_chain,
        lig_resid=lig_resid,
    )


def score_geometry(cif_path: Path, info: dict, amap: g.AtomMap) -> dict:
    """Our calibrated metrics on the generated structure."""
    struct = get_structure(CIFFile.read(cif_path), model=1)
    kw = dict(
        lig_chain=info["lig_chain"],
        lig_resid=info["lig_resid"],
        chain=info["protein_chain"],
    )
    out = {}
    out["base_to_abstraction"] = g.base_to_abstraction(
        struct, amap, base_resid=info["catalytic_resid"]["base"], **kw
    )
    out["donor_to_oxyanion"] = g.donor_to_oxyanion(
        struct, amap, donor_resid=info["catalytic_resid"]["donor"], **kw
    )
    out["wall_stacking"] = g.wall_stacking(
        struct, amap, wall_resid=info["catalytic_resid"]["wall"], **kw
    )
    out["min_contact"] = g.min_protein_ligand_contact(
        struct, amap, lig_chain=info["lig_chain"], lig_resid=info["lig_resid"]
    )
    # Burial: necessary because the distance metrics above are satisfiable
    # with the ligand sitting on the protein surface (see thresholds.yaml
    # geometry.ligand_burial).
    out.update(
        g.ligand_burial(
            struct, amap,
            lig_chain=info["lig_chain"], lig_resid=info["lig_resid"],
        )
    )
    return out


def apply_motif_gate(rec: DesignRecord, thr: dict) -> tuple[bool, list[str]]:
    """Gate a design on active-site geometry (plan 6.2 acceptance).

    This exists because RFD4 can produce a structurally fine protein whose
    catalytic residues are nowhere near the ligand, with no warning from the
    model and a fully populated conditioning column (see docs/TOOL_NOTES.md,
    "theozyme conditioning recipe"). Measuring the output is the only thing
    that catches it, so this gate runs on every design before it is used.
    """
    geo = thr["geometry"]
    reasons: list[str] = []
    checks = [
        ("geom.base_to_abstraction", geo["base_to_abstraction"], "base->abstraction"),
        ("geom.donor_to_oxyanion", geo["donor_to_oxyanion"], "donor->oxyanion"),
        ("geom.wall_stacking", geo["wall_stacking"], "wall stacking"),
    ]
    for key, spec, label in checks:
        val = rec.scores.get(key)
        if val is None:
            reasons.append(f"{label}: not measured")
            continue
        if val > spec["pass_max"]:
            reasons.append(f"{label} {val:.2f} > pass_max {spec['pass_max']}")
        lo = spec.get("hard_fail_min")
        if lo is not None and val < lo:
            reasons.append(f"{label} {val:.2f} < hard_fail_min {lo} (clash)")

    contact = rec.scores.get("geom.min_contact")
    min_contact = thr["stage1_backbone"]["min_protein_ligand_contact"]
    if contact is not None and contact < min_contact:
        reasons.append(f"min protein-ligand contact {contact:.2f} < {min_contact}")

    # --- burial: the ligand must be in a pocket, not on the surface --------
    bur = geo.get("ligand_burial")
    if bur:
        rs = rec.scores.get("geom.rel_sasa")
        if rs is not None and rs > bur["rel_sasa_pass_max"]:
            reasons.append(
                f"ligand rel_sasa {rs:.3f} > {bur['rel_sasa_pass_max']} "
                "(solvent-exposed, no pocket)"
            )
        enc = rec.scores.get("geom.enclosure")
        if enc is not None and enc < bur["enclosure_pass_min"]:
            reasons.append(
                f"enclosure {enc:.3f} < {bur['enclosure_pass_min']} "
                "(ligand not surrounded)"
            )
        nc = rec.scores.get("geom.n_contacts")
        if nc is not None and nc < bur["n_contacts_pass_min"]:
            reasons.append(
                f"only {nc:.0f} protein contacts < {bur['n_contacts_pass_min']}"
            )

    if rec.length is not None:
        lmax = thr["stage6_compliance"]["length_max"]
        if rec.length > lmax:
            reasons.append(f"length {rec.length} > {lmax}")

    return (not reasons), reasons


def main(root: Path, rfd4_out: Path, stage_dir: Path) -> None:
    amap = g.load_atom_map(root / "configs" / "atom_map.yaml")
    thr = g.load_thresholds(root / "configs" / "thresholds.yaml")
    stage_dir.mkdir(parents=True, exist_ok=True)
    manifest = stage_dir / "manifest.jsonl"
    pdb_dir = stage_dir / "backbones_pdb"
    pdb_dir.mkdir(exist_ok=True)
    fasta_a = stage_dir / "armA_rfd4_sequences.fasta"

    cifs = sorted(rfd4_out.rglob("*_generated.cif"))
    print(f"found {len(cifs)} generated CIF(s) under {rfd4_out}")
    if not cifs:
        raise SystemExit("no RFD4 outputs to extract")

    fa_lines, n_ok, n_fail, n_gated = [], 0, 0, 0
    run_id = rfd4_out.name

    for n, cif in enumerate(cifs):
        did = design_id("kemp_theozyme_unindexed", n, 0)
        rec = DesignRecord(
            design_id=did,
            run_id=run_id,
            cell={"donor": "a", "wall": "1", "rasa": float("nan")},
            provenance={"stage1_seed": n, "stage2_seed": 0,
                        "parent_id": None, "generation": 0},
            paths={"rfd4_cif": str(cif)},
        )
        try:
            info = parse_design(cif, amap)
            rec.sequence = info["sequence"]
            rec.length = info["length"]
            rec.catalytic_residues = {
                "base": info["catalytic"]["base"],
                "donor": info["catalytic"]["donor"],
                "wall": [info["catalytic"]["wall"]],
            }
            rec.add_scores("rfd4", {"length": info["length"]})
            rec.add_scores("geom", score_geometry(cif, info, amap))

            # PDB hand-off for LigandMPNN (arm B)
            struct = get_structure(CIFFile.read(cif), model=1)
            pdb_path = pdb_dir / f"{did}.pdb"
            pf = PDBFile()
            pf.set_structure(struct)
            pf.write(pdb_path)
            rec.paths["backbone_pdb"] = str(pdb_path)

            passed, reasons = apply_motif_gate(rec, thr)
            rec.set_gate("stage1_motif", passed, reasons)

            # Only designs that actually kept the active site go downstream.
            if passed:
                fa_lines.append(
                    f">{did} armA_rfd4 len={info['length']} "
                    f"base={info['catalytic']['base']} "
                    f"donor={info['catalytic']['donor']} "
                    f"wall={info['catalytic']['wall']}\n{info['sequence']}"
                )
                n_ok += 1
            else:
                n_gated += 1

            tag = "ok  " if passed else "GATE"
            print(
                f"  [{tag}] {did} len={info['length']:3d} "
                f"base={info['catalytic']['base']:3d} "
                f"donor={info['catalytic']['donor']:3d} "
                f"wall={info['catalytic']['wall']:3d}  "
                f"b={rec.scores['geom.base_to_abstraction']:5.2f} "
                f"d={rec.scores['geom.donor_to_oxyanion']:5.2f} "
                f"w={rec.scores['geom.wall_stacking']:5.2f} "
                f"sasa={rec.scores.get('geom.rel_sasa', float('nan')):5.3f} "
                f"encl={rec.scores.get('geom.enclosure', float('nan')):5.3f}"
                + ("" if passed else f"   <- {reasons[0]}")
            )
        except Exception as exc:  # noqa: BLE001 -- record and continue
            rec.fail("extract", exc)
            n_fail += 1
            print(f"  [FAIL] {did} {cif.name}: {type(exc).__name__}: {exc}")

        rec.write(stage_dir / "designs" / did)
        append_manifest(manifest, rec)

    fasta_a.write_text("\n".join(fa_lines) + ("\n" if fa_lines else ""))
    print(f"\npassed_motif_gate={n_ok}  gated_out={n_gated}  failed={n_fail}")
    print(f"arm A sequences: {fasta_a}")
    print(f"backbone PDBs:   {pdb_dir}")
    print(f"manifest:        {manifest}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
