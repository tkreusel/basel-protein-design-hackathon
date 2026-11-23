"""Inspect the Milestone 1 benchmark structures.

Read-only reconnaissance, run before any geometry code is written. Answers:
what chains/ligands does each structure hold, is 6NT actually present, and
which residues surround it? Plan 6.1 depends on the real answers, not on
assumed residue numbering.

Usage: python -I explore_benchmark.py <benchmark/pdb dir>
"""

import sys
from pathlib import Path

import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure

# The catalytic set reported for the HG3 series, used only as a hypothesis to
# check against the deposited coordinates -- never as an assumption.
HG3_CATALYTIC_HINTS = {127: "Asp base", 50: "Gln oxyanion donor", 44: "Trp wall"}

SOLVENT = {"HOH", "DOD", "WAT"}


def main(pdb_dir: Path) -> None:
    for cif_path in sorted(pdb_dir.glob("*.cif")):
        pdb_id = cif_path.stem
        print(f"\n{'=' * 70}\n{pdb_id}\n{'=' * 70}")
        try:
            cif = CIFFile.read(cif_path)
            # model 1 only; keep every hetero atom so ligands survive
            struct = get_structure(cif, model=1)
        except Exception as exc:  # noqa: BLE001 - report and continue (plan 0.2)
            print(f"  FAILED to parse: {type(exc).__name__}: {exc}")
            continue

        chains = sorted(set(struct.chain_id))
        print(f"  atoms: {struct.array_length()}   chains: {chains}")

        hetero = struct[struct.hetero]
        het_names = sorted(
            {r for r in hetero.res_name if r not in SOLVENT}
        )
        print(f"  non-solvent hetero residues: {het_names}")

        # Is the transition-state analogue there?
        for lig in ("6NT", "XGU"):
            mask = struct.res_name == lig
            if mask.any():
                lig_atoms = struct[mask]
                copies = sorted(
                    set(zip(lig_atoms.chain_id, lig_atoms.res_id))
                )
                print(f"  >> {lig} present: {len(copies)} copy/copies {copies}")
                print(
                    f"     atom names: {list(lig_atoms.atom_name[: len(lig_atoms) // max(len(copies), 1)])}"
                )

                # What protein residues line the pocket around copy 1?
                first_chain, first_res = copies[0]
                centre = struct[
                    mask
                    & (struct.chain_id == first_chain)
                    & (struct.res_id == first_res)
                ]
                prot = struct[~struct.hetero]
                if prot.array_length():
                    d = np.linalg.norm(
                        prot.coord[:, None, :] - centre.coord[None, :, :], axis=-1
                    ).min(axis=1)
                    near = sorted(
                        {
                            (c, int(r), n)
                            for c, r, n in zip(
                                prot.chain_id[d < 4.5],
                                prot.res_id[d < 4.5],
                                prot.res_name[d < 4.5],
                            )
                        },
                        key=lambda t: (t[0], t[1]),
                    )
                    print(f"     pocket residues within 4.5 A ({len(near)}):")
                    for c, r, n in near:
                        hint = HG3_CATALYTIC_HINTS.get(r, "")
                        flag = f"   <-- {hint}" if hint else ""
                        print(f"       {c}{r:>4} {n}{flag}")

        # Report the hypothesised catalytic positions as actually deposited.
        print("  catalytic-position check (chain A):")
        for resid, role in sorted(HG3_CATALYTIC_HINTS.items()):
            m = (struct.chain_id == "A") & (struct.res_id == resid) & (~struct.hetero)
            if m.any():
                names = sorted(set(struct.res_name[m]))
                atoms = list(struct.atom_name[m])
                print(f"    A{resid}: {names} ({role}) atoms={atoms}")
            else:
                print(f"    A{resid}: ABSENT  (expected {role})")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
