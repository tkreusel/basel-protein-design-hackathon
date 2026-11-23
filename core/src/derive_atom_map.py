"""Derive the 6NT atom map from deposited coordinates (plan 6.1).

The plan is explicit: do not hard-code atom names into geometry.py, and do not
trust a chemical argument about which 6NT atom plays which catalytic role.
Derive it from the crystal structures, which show where the catalytic base and
the oxyanion donor actually point.

Chemistry under test (to be confirmed or refuted by the numbers):
  6NT is 6-nitrobenzotriazole, a transition-state ANALOGUE of
  5-nitrobenzisoxazole. Mapping the triazole onto the isoxazole ring:
      6NT N3 (bonded C3A, carries H3)  <->  substrate C3  (H abstracted by base)
      6NT N1 (bonded C7A and N2)       <->  substrate O1  (becomes phenolate O)
  Prediction: Asp127 carboxylate sits near N3; the oxyanion donor (Gln50 NE2
  when present) sits near N1 -- NOT near the nitro group.

The nitro group is the alternative oxyanion hypothesis (plan cell family 'b':
Tyr/Ser/Thr hydroxyl -> nitro group), so we measure it too and let the
distances decide.

Usage: python -I derive_atom_map.py <benchmark/pdb dir>
"""

import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

LIGAND = "6NT"

# Which variant each structure is, for reading the table. Position 50 identity
# is read from the coordinates, not assumed.
VARIANT = {
    "7K4Q": "HG3            (kcat/KM 146)",
    "7K4X": "HG3.7          (27,000)",
    "7K4Z": "HG3.17         (126,000)",
    "8RD5": "HG3.R5         (n/a)",
    "6SS1": "HG3.17 Q50A    (oxyanion ablated)",
    "6SS3": "HG3.17 Q50K    (oxyanion ablated)",
}

BASE_RESID = 127          # Asp catalytic base
DONOR_RESID = 50          # oxyanion donor position (Gln/Lys/Ala across series)
WALL_RESID = 44           # Trp pocket wall

# Candidate ligand acceptor/donor atoms, by hypothesis.
TRIAZOLE = ["N1", "N2", "N3"]
NITRO = ["NO1", "O11", "O21"]

# Side-chain polar tips per residue type at position 50.
DONOR_TIPS = {
    "GLN": ["NE2", "OE1"],
    "ASN": ["ND2", "OD1"],
    "LYS": ["NZ"],
    "TYR": ["OH"],
    "SER": ["OG"],
    "THR": ["OG1"],
    "ALA": [],
}
BASE_TIPS = ["OD1", "OD2"]


def tip_coords(struct, chain, resid, names):
    out = {}
    for n in names:
        m = (
            (struct.chain_id == chain)
            & (struct.res_id == resid)
            & (struct.atom_name == n)
            & (~struct.hetero)
        )
        if m.any():
            out[n] = struct.coord[m][0]
    return out


def main(pdb_dir: Path) -> None:
    print(
        "Distances in Angstrom. '*' marks the closest ligand atom for that "
        "protein tip.\nAn H-bond is ~2.6-3.2 A between heavy atoms.\n"
    )
    for cif_path in sorted(pdb_dir.glob("*.cif")):
        pdb_id = cif_path.stem
        cif = CIFFile.read(cif_path)
        struct = get_structure(cif, model=1)

        lig_mask = struct.res_name == LIGAND
        copies = sorted(set(zip(struct.chain_id[lig_mask], struct.res_id[lig_mask])))

        print(f"{'=' * 78}")
        print(f"{pdb_id}   {VARIANT.get(pdb_id, '?')}")
        print(f"{'=' * 78}")

        for chain, resid in copies:
            lig = struct[
                lig_mask & (struct.chain_id == chain) & (struct.res_id == resid)
            ]
            lig_xyz = {n: c for n, c in zip(lig.atom_name, lig.coord)}

            # identity at the donor position, read from coordinates
            dm = (
                (struct.chain_id == chain)
                & (struct.res_id == DONOR_RESID)
                & (~struct.hetero)
            )
            donor_res = sorted(set(struct.res_name[dm]))
            donor_res = donor_res[0] if donor_res else "ABSENT"

            print(f"\n  chain {chain}, {LIGAND} {resid}   position 50 = {donor_res}")

            # --- base -> abstraction centre -------------------------------
            base = tip_coords(struct, chain, BASE_RESID, BASE_TIPS)
            if base:
                print(f"    Asp{BASE_RESID} carboxylate -> triazole N:")
                for tip, bxyz in base.items():
                    ds = {
                        n: float(np.linalg.norm(bxyz - lig_xyz[n]))
                        for n in TRIAZOLE
                        if n in lig_xyz
                    }
                    best = min(ds, key=ds.get)
                    row = "  ".join(
                        f"{n}={ds[n]:5.2f}{'*' if n == best else ' '}" for n in TRIAZOLE
                    )
                    print(f"      {tip}:  {row}")

            # --- donor -> oxyanion centre vs nitro ------------------------
            tips = DONOR_TIPS.get(donor_res, [])
            donor = tip_coords(struct, chain, DONOR_RESID, tips)
            if not donor:
                print(f"    position 50 ({donor_res}): no polar tip -- ablated")
            else:
                print(f"    {donor_res}{DONOR_RESID} tip -> triazole N  |  nitro:")
                for tip, dxyz in donor.items():
                    tri = {
                        n: float(np.linalg.norm(dxyz - lig_xyz[n]))
                        for n in TRIAZOLE
                        if n in lig_xyz
                    }
                    nit = {
                        n: float(np.linalg.norm(dxyz - lig_xyz[n]))
                        for n in NITRO
                        if n in lig_xyz
                    }
                    allds = {**tri, **nit}
                    best = min(allds, key=allds.get)
                    trirow = "  ".join(
                        f"{n}={tri[n]:5.2f}{'*' if n == best else ' '}"
                        for n in TRIAZOLE
                        if n in tri
                    )
                    nitrow = "  ".join(
                        f"{n}={nit[n]:5.2f}{'*' if n == best else ' '}"
                        for n in NITRO
                        if n in nit
                    )
                    print(f"      {tip}:  {trirow}   |  {nitrow}")

            # --- wall stacking -------------------------------------------
            wm = (
                (struct.chain_id == chain)
                & (struct.res_id == WALL_RESID)
                & (~struct.hetero)
                & np.isin(
                    struct.atom_name,
                    ["CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
                )
            )
            ring = [n for n in ("C3A", "C4", "C5", "C6", "C7", "C7A") if n in lig_xyz]
            if wm.any() and ring:
                wc = struct.coord[wm].mean(axis=0)
                lc = np.array([lig_xyz[n] for n in ring]).mean(axis=0)
                print(
                    f"    Trp{WALL_RESID} ring centroid -> benzo centroid: "
                    f"{float(np.linalg.norm(wc - lc)):5.2f}"
                )
        print()


if __name__ == "__main__":
    main(Path(sys.argv[1]))
