"""Ligand definitions for the Kemp reaction, using the real CCD codes and the CCD's own atom names.

    arm    CCD code  compound                                   role in the reaction
    tsa    6NT       6-nitrobenzotriazole                        transition-state analogue
    sub    H5J       5-nitro-1,2-benzoxazole                     substrate
    prod   6VP       5-nitro-2-hydroxybenzonitrile (neutral!)    product (the CCD entry is the PHENOL, not the phenolate anion)

Because the codes are genuine CCD components and we keep the CCD atom names, tmol / atomworks / PLACER / Rosetta all see
one consistent definition (an earlier version used stand-in codes whose names clashed with real CCD entries and were silently
re-protonated). Incoming coordinates (Boltz, crystals, PLACER) are renamed to the CCD names by graph matching, so any naming works.

Role atoms (per ligand) are the only place atom names appear; geometry code reads ROLES, never literal names.

    role                  6NT     H5J     6VP
    abstraction centre    N3(-H3) CAE(-H2) CAE (nitrile C; no H)
    leaving-group atom    N1      OAH     OAB (phenol O)
    ring partner of C3    N2      NAG     NAA (nitrile N)
    C3a / C7a             C3A/C7A CAJ/CAK CAJ/CAI
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

LIGAND_CODES = ("6NT", "H5J", "6VP")

# source-ligand atom -> 6NT atom, used to superpose substrate/product onto a TSA pose (same-pocket comparisons).
# 6VP's nitrile N is not in a ring, so it is left out of the product superposition.
CORRESPONDENCE_TO_TSA = {
    "H5J": {"CAJ": "C3A", "CAF": "C4", "CAI": "C5", "CAC": "C6", "CAD": "C7", "CAK": "C7A",
            "NAL": "NO1", "OAA": "O11", "OAB": "O21", "OAH": "N1", "NAG": "N2", "CAE": "N3"},
    "6VP": {"CAJ": "C3A", "CAH": "C4", "CAK": "C5", "CAG": "C6", "CAF": "C7", "CAI": "C7A",
            "NAL": "NO1", "OAD": "O11", "OAC": "O21", "OAB": "N1", "CAE": "N3"},
}

ROLES = {
    "6NT": dict(abstraction="N3", abstraction_H="H3", oxyanion="N1", oxyanion_partner="C7A",
                ring_partner="N2", ring_attached="C3A", nitro_N="NO1", nitro_O=["O11", "O21"],
                benzo=["C3A", "C4", "C5", "C6", "C7", "C7A"]),
    "H5J": dict(abstraction="CAE", abstraction_H="H2", oxyanion="OAH", oxyanion_partner="CAK",
                ring_partner="NAG", ring_attached="CAJ", nitro_N="NAL", nitro_O=["OAA", "OAB"],
                benzo=["CAJ", "CAF", "CAI", "CAC", "CAD", "CAK"]),
    "6VP": dict(abstraction="CAE", abstraction_H=None, oxyanion="OAB", oxyanion_partner="CAI",
                ring_partner="NAA", ring_attached="CAJ", nitro_N="NAL", nitro_O=["OAD", "OAC"],
                benzo=["CAJ", "CAH", "CAK", "CAG", "CAF", "CAI"]),
}
ARM_TO_RESNAME = {"tsa": "6NT", "sub": "H5J", "prod": "6VP"}
RESNAME_TO_ARM = {v: k for k, v in ARM_TO_RESNAME.items()}


@dataclass
class LigandTemplate:
    resname: str
    mol: object  # rdkit Mol with explicit H, 3D conformer, atom property 'name'
    names: list[str]

    def coords(self) -> np.ndarray:
        return np.array(self.mol.GetConformer().GetPositions())


def _from_ccd(code: str) -> LigandTemplate:
    """Template from the CCD entry (biotite's bundled copy): CCD atom names, ideal coordinates, explicit H, formal charges."""
    import biotite.structure.info as info
    from rdkit import Chem

    at = info.get_from_ccd("chem_comp_atom", code)
    bd = info.get_from_ccd("chem_comp_bond", code)
    m = Chem.RWMol()
    idx, pos = {}, []
    for name, el, ch, x, y, z in zip(at["atom_id"].as_array(), at["type_symbol"].as_array(), at["charge"].as_array(),
                                     at["pdbx_model_Cartn_x_ideal"].as_array(), at["pdbx_model_Cartn_y_ideal"].as_array(),
                                     at["pdbx_model_Cartn_z_ideal"].as_array()):
        a = Chem.Atom(str(el).capitalize())
        a.SetFormalCharge(int(ch))
        a.SetProp("name", str(name))
        idx[str(name)] = m.AddAtom(a)
        pos.append((float(x), float(y), float(z)))
    order = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE}
    for a1, a2, o in zip(bd["atom_id_1"].as_array(), bd["atom_id_2"].as_array(), bd["value_order"].as_array()):
        m.AddBond(idx[str(a1)], idx[str(a2)], order[str(o)])
    conf = Chem.Conformer(m.GetNumAtoms())
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, p)
    m.AddConformer(conf)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return LigandTemplate(code, mol, [a.GetProp("name") for a in mol.GetAtoms()])


def build_templates(ligand_dir=None) -> dict[str, LigandTemplate]:
    """{'6NT','H5J','6VP'} -> template. `ligand_dir` kept for call-compatibility (CCD data comes from biotite)."""
    return {c: _from_ccd(c) for c in LIGAND_CODES}


# ---- file writers (mol2 with MMFF94 partial charges -> tmol params) ----------------------
def write_mol2(t: LigandTemplate, path: Path) -> None:
    from rdkit.Chem import AllChem

    mol = t.mol
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        raise RuntimeError(f"MMFF typing failed for {t.resname}")
    conf = mol.GetConformer()
    L = ["@<TRIPOS>MOLECULE", t.resname, f" {mol.GetNumAtoms()} {mol.GetNumBonds()} 1 0 0", "SMALL",
         "USER_CHARGES", "", "@<TRIPOS>ATOM"]
    for i, a in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        sym = a.GetSymbol()
        typ = {"C": "C.ar" if a.GetIsAromatic() else ("C.1" if a.GetDegree() == 2 else "C.3"),
               "N": "N.ar" if a.GetIsAromatic() else ("N.1" if a.GetDegree() == 1 else "N.2"),
               "O": "O.2", "H": "H"}[sym]
        L.append(f"{i+1:>5} {t.names[i]:<6} {p.x:9.4f} {p.y:9.4f} {p.z:9.4f} {typ:<6} 1 {t.resname} "
                 f"{props.GetMMFFPartialCharge(i):8.4f}")
    L.append("@<TRIPOS>BOND")
    for j, b in enumerate(mol.GetBonds()):
        bt = "ar" if b.GetIsAromatic() else {1.0: "1", 2.0: "2", 3.0: "3"}[b.GetBondTypeAsDouble()]
        L.append(f"{j+1:>5} {b.GetBeginAtomIdx()+1:>4} {b.GetEndAtomIdx()+1:>4} {bt}")
    Path(path).write_text("\n".join(L) + "\n")


def write_sdf(t: LigandTemplate, path: Path) -> None:
    from rdkit import Chem

    Chem.MolToMolFile(t.mol, str(path))


# ---- naming incoming coordinates by graph matching --------------------------------------
_COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66}


def _graph_from_coords(elements, xyz):
    n = len(elements)
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    r = np.array([_COV[e] for e in elements])
    adj = (d < (r[:, None] + r[None] + 0.45)) & ~np.eye(n, dtype=bool)
    return adj


def match_names(template: LigandTemplate, elements: list[str], xyz: np.ndarray) -> list[str | None]:
    """Return canonical atom names for each input heavy atom (None for unmatched / H).

    Connectivity is inferred from interatomic distances (covalent radii + 0.45 A) and
    matched to the template heavy-atom graph by element-labelled subgraph isomorphism, with
    ties (nitro O's) resolved by lowest RMSD after superposition.
    """
    import networkx as nx
    from networkx.algorithms import isomorphism as iso

    heavy = [i for i, e in enumerate(elements) if e != "H"]
    el = [elements[i] for i in heavy]
    adj = _graph_from_coords(el, xyz[heavy])
    G = nx.Graph()
    for i, e in enumerate(el):
        G.add_node(i, e=e)
    for i, j in zip(*np.nonzero(np.triu(adj))):
        G.add_edge(int(i), int(j))
    T = nx.Graph()
    tm = template.mol
    t_heavy = [a.GetIdx() for a in tm.GetAtoms() if a.GetSymbol() != "H"]
    for k, i in enumerate(t_heavy):
        T.add_node(k, e=tm.GetAtomWithIdx(i).GetSymbol())
    pos = {i: k for k, i in enumerate(t_heavy)}
    for b in tm.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a in pos and c in pos:
            T.add_edge(pos[a], pos[c])
    gm = iso.GraphMatcher(G, T, node_match=lambda a, b: a["e"] == b["e"])
    best, best_rmsd = None, 1e9
    t_xyz = template.coords()[t_heavy]
    for m in gm.isomorphisms_iter():  # G node -> T node
        if len(m) != len(el):
            continue
        src = np.array([xyz[heavy[g]] for g in m])
        dst = np.array([t_xyz[t] for t in m.values()])
        rmsd = _kabsch_rmsd(src, dst)
        if rmsd < best_rmsd:
            best, best_rmsd = m, rmsd
    if best is None:
        return [None] * len(elements)
    out: list[str | None] = [None] * len(elements)
    for g, t in best.items():
        out[heavy[g]] = template.names[t_heavy[t]]
    return out


def _kabsch(src: np.ndarray, dst: np.ndarray):
    """R, t minimising |R src + t - dst|."""
    sc, dc = src.mean(0), dst.mean(0)
    H = (src - sc).T @ (dst - dc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, dc - R @ sc


def _kabsch_rmsd(src, dst) -> float:
    R, t = _kabsch(src, dst)
    return float(np.sqrt(((src @ R.T + t - dst) ** 2).sum(1).mean()))


def _best_corr(corr: dict, src_xyz: dict, dst_xyz: dict) -> dict:
    """The two nitro oxygens are chemically equivalent: pick the pairing (O11<->O21 or swapped) with the lower RMSD."""
    o_src = [k for k, v in corr.items() if v in ("O11", "O21")]
    if len(o_src) != 2:
        return corr
    alt = dict(corr)
    alt[o_src[0]], alt[o_src[1]] = corr[o_src[1]], corr[o_src[0]]
    def rmsd(c):
        A = np.array([src_xyz[k] for k in c]); B = np.array([dst_xyz[v] for v in c.values()])
        return _kabsch_rmsd(A, B)
    return corr if rmsd(corr) <= rmsd(alt) else alt


def place_onto_tsa(template: LigandTemplate, tsa_xyz_by_name: dict[str, np.ndarray]):
    """Superpose a H5J/6VP template onto the coordinates of a 6NT pose.

    Returns (names, xyz[n,3], rmsd_of_mapped_atoms). Hydrogens come along rigidly.
    """
    names = template.names
    t_xyz = template.coords()
    src_xyz = dict(zip(names, t_xyz))
    corr = _best_corr(CORRESPONDENCE_TO_TSA[template.resname], src_xyz, tsa_xyz_by_name)
    src = np.array([src_xyz[k] for k in corr])
    dst = np.array([tsa_xyz_by_name[v] for v in corr.values()])
    R, t = _kabsch(src, dst)
    return list(names), t_xyz @ R.T + t, _kabsch_rmsd(src, dst)
