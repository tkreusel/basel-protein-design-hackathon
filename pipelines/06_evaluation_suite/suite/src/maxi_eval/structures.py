"""Read crystal / NMR / Boltz / PLACER structures into one standard 'Replicate' form.

A Replicate = one protein chain + one ligand copy, heavy atoms (and H if present), ligand atoms
renamed to the canonical scheme in ligands.py. Everything downstream (geometry, PLACER input,
tmol input) consumes Replicates, so the source of a structure never matters after this module.

Sources:
  crystal  : several ligand copies per file (chains A/B) -> one Replicate per copy; waters kept
             (coordinates only) for validating the water-site proxy.
  nmr      : one Replicate per model.
  boltz    : one Replicate per model_k.cif (protein chain A + ligand chain B).
  placer   : multi-MODEL PDB written by PLACER.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import ligands as L

ADDITIVES = {"SO4", "GOL", "ACT", "EDO", "PEG", "PGE", "PG4", "DMS", "FMT", "CL", "NA", "MG", "ZN", "CA", "K", "TRS", "MPD", "BME", "IMD"}


@dataclass
class Replicate:
    rep_id: str                 # e.g. "4BS0:m1:A1305"
    entry: str
    arm: str                    # tsa | sub | prod
    source: str                 # crystal | nmr | boltz | placer
    protein: object             # biotite AtomArray (one chain)
    ligand: object              # biotite AtomArray, canonical names
    waters: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    meta: dict = field(default_factory=dict)

    @property
    def lig_resname(self) -> str | None:
        return L.ARM_TO_RESNAME.get(self.arm)      # None for the ligand-free 'apo' arm


def _read_models(path: Path):
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure.io.pdb as pdb

    path = Path(path)
    if path.suffix in (".cif", ".mmcif"):
        f = pdbx.CIFFile.read(path)
        try:
            stack = pdbx.get_structure(f, altloc="first", use_author_fields=True, include_bonds=False)
        except Exception:  # noqa: BLE001 - e.g. no auth_* columns (some predictors)
            stack = pdbx.get_structure(f, altloc="first", use_author_fields=False, include_bonds=False)
    else:
        stack = pdb.PDBFile.read(path).get_structure(altloc="first", include_bonds=False)
    return [stack] if not hasattr(stack, "stack_depth") else [stack[i] for i in range(stack.stack_depth())]


def extract_replicates(path: Path, entry: str, arm: str, source: str,
                       templates: dict[str, L.LigandTemplate], model_ids: list[int] | None = None,
                       keep_h: bool = False) -> list[Replicate]:
    import biotite.structure as struc

    tpl = templates[L.ARM_TO_RESNAME[arm]]
    t_el = sorted(a.GetSymbol() for a in tpl.mol.GetAtoms() if a.GetSymbol() != "H")
    reps: list[Replicate] = []
    models = _read_models(path)
    for mi, st in enumerate(models, start=1):
        if model_ids is not None and mi not in model_ids:
            continue
        if not keep_h:
            st = st[st.element != "H"]
        is_prot = struc.filter_amino_acids(st)
        prot = st[is_prot]
        non = st[~is_prot & ~struc.filter_solvent(st)]
        waters = st[struc.filter_solvent(st) & (st.element == "O")].coord
        # candidate ligand instances
        keys = []
        for ch, rid, rn in dict.fromkeys(zip(non.chain_id, non.res_id, non.res_name)):
            if rn in ADDITIVES:
                continue
            m = (non.chain_id == ch) & (non.res_id == rid) & (non.res_name == rn)
            heavy = non[m]
            heavy = heavy[heavy.element != "H"]
            if sorted(heavy.element) == t_el:
                keys.append((ch, rid, rn))
        for ch, rid, rn in keys:
            m = (non.chain_id == ch) & (non.res_id == rid) & (non.res_name == rn)
            lig = non[m]
            lig_heavy = lig[lig.element != "H"]
            names = L.match_names(tpl, list(lig_heavy.element), lig_heavy.coord)
            if any(n is None for n in names):
                continue
            lig_c = lig_heavy.copy()
            lig_c.atom_name = np.array(names)
            lig_c.res_name = np.full(len(lig_c), tpl.resname)
            lig_c.chain_id = np.full(len(lig_c), "B")
            lig_c.res_id = np.ones(len(lig_c), dtype=int)
            lig_c.hetero = np.ones(len(lig_c), dtype=bool)
            # protein chain = chain closest to the ligand
            d = np.linalg.norm(prot.coord[:, None] - lig_heavy.coord[None], axis=-1).min(1)
            chain = prot.chain_id[int(np.argmin(d))]
            pc = prot[prot.chain_id == chain].copy()
            pc.chain_id = np.full(len(pc), "A")
            near_w = waters[np.linalg.norm(waters[:, None] - lig_heavy.coord.mean(0), axis=-1).ravel() < 15.0] \
                if len(waters) else np.zeros((0, 3))
            rid_s = f"{chain}{rid}"
            reps.append(Replicate(
                rep_id=f"{entry}:{arm}:m{mi}:{rid_s}" if len(models) > 1 else f"{entry}:{arm}:{rid_s}",
                entry=entry, arm=arm, source=source, protein=pc, ligand=lig_c, waters=near_w,
                meta=dict(file=str(path), model=mi, src_chain=str(chain), src_lig_resid=int(rid))))
    return reps


def write_complex_pdb(rep: Replicate, path: Path, with_h: bool = False) -> None:
    import biotite.structure.io.pdb as pdb
    import biotite.structure as struc

    p = rep.protein if with_h else rep.protein[rep.protein.element != "H"]
    l = rep.ligand if with_h else rep.ligand[rep.ligand.element != "H"]
    both = p + l if len(l) else p
    f = pdb.PDBFile()
    f.set_structure(both)
    f.write(str(path))


def split_placer_models(pdb_path: Path, entry: str, arm: str, templates) -> list[Replicate]:
    """PLACER writes one multi-MODEL PDB; ligand keeps its input name (canonical, chain B)."""
    return extract_replicates(pdb_path, entry, arm, "placer", templates)


def make_ligand_array(names, xyz, resname, elements):
    import biotite.structure as struc

    a = struc.AtomArray(len(names))
    a.coord = np.asarray(xyz, dtype=np.float32)
    a.atom_name = np.array(names)
    a.element = np.array(elements)
    a.res_name = np.full(len(names), resname)
    a.chain_id = np.full(len(names), "B")
    a.res_id = np.ones(len(names), dtype=int)
    a.hetero = np.ones(len(names), dtype=bool)
    return a


def extract_apo(path: Path, entry: str, source: str = "boltz") -> "Replicate":
    """Ligand-free (apo) prediction: protein chain only, empty ligand."""
    import biotite.structure as struc
    st = _read_models(Path(path))[0]
    st = st[st.element != "H"]
    prot = st[struc.filter_amino_acids(st)]
    chain = prot.chain_id[0]
    prot = prot[prot.chain_id == chain].copy()
    prot.chain_id = np.full(len(prot), "A")
    return Replicate(rep_id=f"{entry}:apo", entry=entry, arm="apo", source=source, protein=prot,
                     ligand=struc.AtomArray(0), meta=dict(file=str(path), model=1))
