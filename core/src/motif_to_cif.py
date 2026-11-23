"""Convert the Kemp motif PDB to the mmCIF that RFD4 requires.

Two hard facts from docs/reference/yaml-json.md that make this non-trivial:

  1. RFD4's `input` must be .cif, .cif.gz or .bcif. A .pdb is rejected outright,
     so the organisers' HG3_17_kemp_motif.pdb cannot be used directly.

  2. `res_id` in both `contig` and `conditions` resolves against the mmCIF
     **label_seq_id**, NOT auth_seq_id. Our motif uses author numbering
     (Trp44, Gln50, Asp127). If label_seq_id comes out as 1,2,3 then every
     selection written as `sel('A/*/127')` would silently select nothing or the
     wrong residue -- exactly the class of silent error the plan warns about.

This script therefore WRITES the CIF and then PRINTS the numbering that RFD4
will actually see, so the condition spec can be written against observed values
rather than assumed ones. It does not guess.

Usage: python -I motif_to_cif.py <motif.pdb> <out.cif>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, set_structure, get_structure


def main(pdb_path: Path, cif_path: Path, ligand_chain: str = "L") -> None:
    struct = PDBFile.read(pdb_path).get_structure(model=1)
    print(f"read {pdb_path}: {struct.array_length()} atoms")

    # Move the ligand onto its own chain.
    #
    # RFD4's contig puts a chainbreak (/0) before the ligand so it becomes a
    # separate output chain. If the ligand still carries the protein's chain id
    # the contig produces two output chains both called 'A' and apply_contig
    # rejects it:
    #   ValueError: Contig produces two separate output chains both named 'A'.
    # The MCSA benchmark inputs do the same thing -- their ligands sit on their
    # own chains (ORO on C, ZN on B).
    lig = struct.hetero
    if lig.any():
        moved = sorted(set(struct.chain_id[lig]))
        struct.chain_id[lig] = ligand_chain
        print(f"moved {int(lig.sum())} hetero atoms from chain(s) "
              f"{moved} to chain {ligand_chain!r}")

    # Report what we started from, per residue.
    print("\ninput (author) numbering:")
    seen = []
    for cid, rid, rname, het in zip(
        struct.chain_id, struct.res_id, struct.res_name, struct.hetero
    ):
        key = (str(cid), int(rid), str(rname), bool(het))
        if key not in seen:
            seen.append(key)
    for cid, rid, rname, het in seen:
        print(f"  {cid}{rid:>5} {rname:<4} {'HETATM' if het else 'ATOM'}")

    cif_path.parent.mkdir(parents=True, exist_ok=True)
    cif = CIFFile()
    set_structure(cif, struct)
    cif.write(cif_path)
    print(f"\nwrote {cif_path}")

    # --- the part that matters: what numbering is actually in the file? ---
    block = cif.block
    atom_site = block["atom_site"]
    cols = set(atom_site.keys())
    print(f"\natom_site columns present: {sorted(cols)}")

    def col(name):
        return list(atom_site[name].as_array(str)) if name in cols else None

    label_seq = col("label_seq_id")
    auth_seq = col("auth_seq_id")
    label_asym = col("label_asym_id")
    auth_asym = col("auth_asym_id")
    comp = col("label_comp_id") or col("auth_comp_id")

    print("\nnumbering as written to the CIF (unique residues):")
    print(f"  {'label_asym':>10} {'label_seq':>10} {'auth_asym':>10} {'auth_seq':>9}  comp")
    rows = []
    n = len(comp)
    for i in range(n):
        row = (
            label_asym[i] if label_asym else "?",
            label_seq[i] if label_seq else "?",
            auth_asym[i] if auth_asym else "?",
            auth_seq[i] if auth_seq else "?",
            comp[i],
        )
        if row not in rows:
            rows.append(row)
    for r in rows:
        print(f"  {r[0]:>10} {r[1]:>10} {r[2]:>10} {r[3]:>9}  {r[4]}")

    # Round-trip: what does a reader get back?
    back = get_structure(CIFFile.read(cif_path), model=1)
    print("\nround-trip read-back (biotite get_structure):")
    seen2 = []
    for cid, rid, rname in zip(back.chain_id, back.res_id, back.res_name):
        k = (str(cid), int(rid), str(rname))
        if k not in seen2:
            seen2.append(k)
    for cid, rid, rname in seen2:
        print(f"  {cid}{rid:>5} {rname}")

    lab = [r[1] for r in rows]
    aut = [r[3] for r in rows]
    print("\nVERDICT:")
    if lab == aut:
        print("  label_seq_id == auth_seq_id -> selections may use 44/50/127 directly")
    else:
        print("  label_seq_id DIFFERS from auth_seq_id")
        print(f"    label: {lab}")
        print(f"    auth : {aut}")
        print("  -> conditions must be written against the label values above,")
        print("     or the CIF must be rewritten so the two agree.")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]),
         sys.argv[3] if len(sys.argv) > 3 else "L")
