"""Rosetta .params files for 6NT / H5J / 6VP (PyRosetta env; molfile_to_params.py is not shipped in the wheel).

Path: SDF -> sdf.MolFileIOReader -> convert_to_ResidueType -> rename atoms to the canonical names -> set MMFF94 partial
charges from the mol2 (identical to the tmol parameterisation) -> write_topology_file.
Run:  /mnt/work/envs/pyrosetta/bin/python -m maxi_eval.build_rosetta_params [ligand_dir]
"""
import sys
from pathlib import Path


def mol2_names_charges(p):
    names, q, on = [], [], False
    for line in Path(p).read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            on = line.startswith("@<TRIPOS>ATOM"); continue
        if on and line.strip():
            f = line.split(); names.append(f[1]); q.append(float(f[8]))
    return names, q


def main(d="ligands"):
    import pyrosetta
    from pyrosetta import rosetta
    pyrosetta.init("-mute all")
    sdf, chem = rosetta.core.chemical.sdf, rosetta.core.chemical
    d = Path(d)
    for k in ("6NT", "H5J", "6VP"):
        names, charges = mol2_names_charges(d / f"{k}.mol2")
        reader = sdf.MolFileIOReader()
        txt = (d / f"{k}.sdf").read_text()
        tmp = d / f"tmp_{k}.sdf"
        tmp.write_text(txt if "$$$$" in txt else txt.rstrip("\n") + "\n$$$$\n")   # RDKit .mol files lack the SDF terminator
        mols = reader.parse_file(str(tmp))
        tmp.unlink()
        rt = sdf.convert_to_ResidueType(mols)
        vds = list(rt.all_atoms())
        assert len(vds) == len(names), (k, len(vds), len(names))
        # SDF atom order == mol2 order == template order. Rename through unique temporaries (C5->C6 collisions).
        for i, v in enumerate(vds):
            rt.rename_atom(v, f"T{i:02d}")
        for v, n, q in zip(vds, names, charges):
            rt.rename_atom(v, n)
            rt.atom(v).charge(q)
        rt.name(k); rt.name3(k); rt.name1("X")
        out = d / f"{k}.params"
        chem.write_topology_file(chem.ResidueType.make(rt), str(out))
        print("wrote", out, "natoms", len(vds), "net charge %.3f" % sum(charges))

if __name__ == "__main__":
    main(*sys.argv[1:])
