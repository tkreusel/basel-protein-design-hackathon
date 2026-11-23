"""Build ligand files once: ligands/{6NT,H5J,6VP}.{sdf,mol2,tmol}  (CCD codes, CCD atom names, ideal coordinates).

6NT = transition-state analogue, H5J = substrate, 6VP = product (the CCD entry is the neutral phenol).
tmol params are written with mode="keep" from a mol2 carrying MMFF94 partial charges (the same charges the Rosetta params use).
Run:  python -m maxi_eval.build_ligands   (RFD4 venv),  then  python -m maxi_eval.build_rosetta_params   (pyrosetta env)
"""
from pathlib import Path
import sys
from . import ligands as L


def main(d="ligands"):
    d = Path(d)
    T = L.build_templates()
    for k, t in T.items():
        L.write_mol2(t, d / f"{k}.mol2"); L.write_sdf(t, d / f"{k}.sdf")
    from tmol.ligand import write_params_from_mol2
    for k in T:
        write_params_from_mol2(d / f"{k}.mol2", d / f"{k}.tmol", res_name=k, mode="keep")
        print("wrote", d / f"{k}.tmol")

if __name__ == "__main__":
    main(*sys.argv[1:])
