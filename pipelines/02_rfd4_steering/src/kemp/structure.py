"""Structure helpers shared by the stage scripts (biotite only; no RFD4 imports)."""
import numpy as np
import biotite.structure as bs
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

from . import config

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
       "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def residue_tips():
    """{'ASP': ('CB','CG','OD1','OD2'), ...} from config/motifs.yaml."""
    return {k: tuple(v["tips"]) for k, v in config.load("motifs")["residues"].items()}


def load_cif(path, extra_fields=None):
    """Heavy atoms with finite coordinates, model 1, label (not author) identifiers."""
    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1, use_author_fields=False, extra_fields=extra_fields or [])
    return atoms[(atoms.element != "H") & ~np.isnan(atoms.coord).any(axis=1)]


def load_reference():
    """Theozyme (HG3.17 motif + 6NT), heavy atoms."""
    ref = pdb.PDBFile.read(str(config.PIPELINE / config.load("motifs")["reference_pdb"])).get_structure(model=1)
    return ref[ref.element != "H"]


def atom_xyz(atoms, names):
    """Coordinates of the first atom with each name, shape (len(names), 3). IndexError if one is missing."""
    return np.array([atoms.coord[atoms.atom_name == n][0] for n in names])


def rmsd_xyz(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def kabsch(P, Q):
    """Rotation R and translation t that move P onto Q (least squares): P @ R.T + t."""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def locate_motif(protein, ligand, motif, reference=None):
    """Find the unindexed motif residues of an RFD4 design.

    The reference ligand is superposed on the design's ligand; for each motif residue type the design
    residue whose tip atoms lie closest to the transformed reference tips is taken.
    Returns [(res_name, res_id or None, tip_rmsd)] in motif order.
    """
    reference = load_reference() if reference is None else reference
    ref_lig, ref_prot = reference[reference.hetero], reference[~reference.hetero]
    shared = [n for n in ref_lig.atom_name if n in set(ligand.atom_name)]
    _, transform = bs.superimpose(atom_xyz(ligand, shared), atom_xyz(ref_lig, shared))
    tips = residue_tips()
    found = []
    for res in config.load("motifs")["motifs"][motif]["residues"]:
        names = tips[res]
        ref_res = ref_prot[ref_prot.res_name == res]
        target = transform.apply(atom_xyz(ref_res, names))
        best = (9.0, None)
        for rid in np.unique(protein.res_id[protein.res_name == res]):
            cand = protein[(protein.res_id == rid) & (protein.res_name == res)]
            try:
                d = rmsd_xyz(target, atom_xyz(cand, names))
            except IndexError:
                continue
            if d < best[0]:
                best = (d, int(rid))
        found.append((res, best[1], best[0]))
    return found


def sequence(protein):
    return "".join(AA3.get(n, "X") for n in protein[protein.atom_name == "CA"].res_name)
