"""Build a His-base variant of the HG3.17 theozyme.

Replaces Asp127 by an imidazole rigid-fitted onto the carboxylate, keeping
Gln50 (oxyanion donor), Trp44 (wall) and 6NT exactly as in the crystal motif.
Same idea as the handoff's GDM Glu->His files, applied to our validated geometry.

Fit: ideal imidazole (taken from the handoff's His theozyme) CB,CG,ND1,CD2 onto
Asp CB,CG,OD2,OD1. OD2 is the carboxylate O proximal to N3 (2.53 A), so ND1
becomes the proton acceptor. The other ring N (NE2) is the one that must carry
the proton, i.e. needs an H-bond acceptor partner (tautomer lock) -- reported.

Usage: python -I make_his_motif.py <ideal His pdb> <motif cif> <out cif>
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from biotite.structure import AtomArray
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure, set_structure


def kabsch(P, Q):
    """Rotation R and translation t so that P @ R + t ~ Q."""
    pc, qc = P.mean(0), Q.mean(0)
    v, _, wt = np.linalg.svd((P - pc).T @ (Q - qc))
    if np.linalg.det(v @ wt) < 0:
        v[:, -1] *= -1
    R = v @ wt
    return R, qc - pc @ R


def main(ideal_pdb: Path, motif_cif: Path, out_cif: Path, base_resid: int = 127) -> None:
    ideal = PDBFile.read(ideal_pdb).get_structure(model=1)
    ideal = ideal[(ideal.res_name == "HIS")]
    motif = get_structure(CIFFile.read(motif_cif), model=1)

    def xyz(a, rid, name):
        m = (a.res_id == rid) & (a.atom_name == name) & (~a.hetero)
        return a.coord[m][0]

    src_names = ["CB", "CG", "ND1", "CD2"]
    tgt_names = ["CB", "CG", "OD2", "OD1"]
    P = np.array([ideal.coord[ideal.atom_name == n][0] for n in src_names])
    Q = np.array([xyz(motif, base_resid, n) for n in tgt_names])
    R, t = kabsch(P, Q)
    rmsd = float(np.sqrt(((P @ R + t - Q) ** 2).sum(1).mean()))

    his_atoms = ["CB", "CG", "ND1", "CD2", "CE1", "NE2"]
    placed = {n: ideal.coord[ideal.atom_name == n][0] @ R + t for n in his_atoms}

    n3 = motif.coord[(motif.res_name == "6NT") & (motif.atom_name == "N3")][0]
    n1 = motif.coord[(motif.res_name == "6NT") & (motif.atom_name == "N1")][0]
    d = lambda a, b: float(np.linalg.norm(a - b))  # noqa: E731
    print(f"fit RMSD {rmsd:.2f} A   (Maxi's Glu->His fits: 0.09 A)")
    print(f"ND1 -> N3 {d(placed['ND1'], n3):.2f} A   (Asp OD2 -> N3 was 2.53; GDM His 3.00)")
    print(f"NE2 -> N3 {d(placed['NE2'], n3):.2f} A   NE2 -> N1 {d(placed['NE2'], n1):.2f} A")

    # steric check of the placed ring against everything else in the motif
    others = motif[~((motif.res_id == base_resid) & (~motif.hetero))]
    for n, c in placed.items():
        dm = np.linalg.norm(others.coord - c, axis=1)
        i = dm.argmin()
        flag = "  <-- CLOSE" if dm[i] < 2.8 and n not in ("CB",) else ""
        print(f"   {n:3} nearest other atom {others.atom_name[i]:4}{others.res_name[i]} {dm[i]:.2f} A{flag}")

    # write the new motif: HIS at the base position, keep Asp backbone atoms
    keep = ~((motif.res_id == base_resid) & (~motif.hetero) &
             np.isin(motif.atom_name, ["CB", "CG", "OD1", "OD2"]))
    new = motif[keep]
    res_mask = (new.res_id == base_resid) & (~new.hetero)
    new.res_name[res_mask] = "HIS"
    add = AtomArray(len(his_atoms))
    for k, n in enumerate(his_atoms):
        add.coord[k] = placed[n]
        add.atom_name[k] = n
        add.element[k] = n[0]
    for attr in ("chain_id", "res_id", "res_name", "hetero"):
        add.set_annotation(attr, np.array([motif.get_annotation(attr)[
            (motif.res_id == base_resid) & (~motif.hetero)][0]] * len(his_atoms)))
    add.res_name[:] = "HIS"
    add.hetero[:] = False
    combined = new + add
    # order: residues by id then keep atoms grouped
    order = np.argsort(combined.res_id, kind="stable")
    combined = combined[order]
    cif = CIFFile()
    set_structure(cif, combined)
    cif.write(out_cif)
    print(f"wrote {out_cif}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
