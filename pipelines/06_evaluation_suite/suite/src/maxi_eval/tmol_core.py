"""tmol (Rosetta beta2016 in PyTorch) helpers: pose building, relax protocols, energy decomposition.

Everything here runs in the RFD4 cu12 venv (tmol 0.2.0). GPU strongly preferred for relax;
scoring-only paths work on CPU (used by the smoke tests).

What tmol can / cannot do is listed in PIPELINE.md section 4. Short version: all energy terms,
block-pair (residue x ligand) energies, relax/repack/mutate, coordinate and geometric constraints
are native; PackStat, holes, shape complementarity, BUNS, enzdes CST files and the cartesian_ddg
*protocol* are not (proxies live in stability.py / cst.py).
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import attrs
import numpy as np
import torch

warnings.filterwarnings("ignore")

from tmol.io import atom_array_from_pose_stack, build_context_from_atom_array, pose_stack_from_atom_array  # noqa: E402
from tmol.kinematics import CartesianMoveMap, FoldForest  # noqa: E402
from tmol.pack import PackerPalette, PackerTask  # noqa: E402
from tmol.pack.rotamer import FixedAAChiSampler, IncludeCurrentSampler  # noqa: E402
from tmol.pack.rotamer.dunbrack import create_dunbrack_sampler_from_database  # noqa: E402
from tmol.database import ParameterDatabase  # noqa: E402
from tmol.relax import fast_relax  # noqa: E402
from tmol.score import ScoreType, beta2016_score_function  # noqa: E402
from tmol.score.constraint._utility import constrain_all_ca  # noqa: E402

LIGAND_PARAMS = ["6NT.tmol", "H5J.tmol", "6VP.tmol"]


def default_device() -> torch.device:
    return torch.device("cuda", torch.cuda.current_device()) if torch.cuda.is_available() else torch.device("cpu")


@dataclass
class Session:
    """One PoseBuildContext + score function, reused across structures (ligands injected once)."""
    device: torch.device
    ctx: object
    sfxn: object
    cst_weight: float

    @classmethod
    def create(cls, ligand_dir: Path, example_protein, templates, device=None, cst_weight: float = 1.0):
        """Context from `example_protein` + ALL three ligands (6NT, H5J, 6VP).

        tmol only registers ligand types that are present in the structure used to build the context, so a
        context built from a 6NT complex cannot later build a substrate complex (verified: tests/smoke_ctx.py).
        """
        from .structures import make_ligand_array

        device = device or default_device()
        prot = example_protein[example_protein.element != "H"]
        arrays = [prot]
        for k, (resname, dx) in enumerate((("6NT", 300.0), ("H5J", 320.0), ("6VP", 340.0)), start=1):
            t = templates[resname]
            els = [a.GetSymbol() for a in t.mol.GetAtoms()]
            keep = [i for i, e in enumerate(els) if e != "H"]
            a = make_ligand_array([t.names[i] for i in keep], t.coords()[keep] + np.array([dx, 0, 0]), resname,
                                  [els[i] for i in keep])
            a.res_id[:] = k
            arrays.append(a)
        dummy = arrays[0]
        for a in arrays[1:]:
            dummy = dummy + a
        files = [str(Path(ligand_dir) / f) for f in LIGAND_PARAMS]
        ctx = build_context_from_atom_array(dummy, device, prepare_ligands=True, strict_ligands=True,
                                            ligand_params_files=files)
        sf = beta2016_score_function(device, param_db=ctx.parameter_database)
        sf.set_weight(ScoreType.constraint, cst_weight)
        return cls(device, ctx, sf, cst_weight)


def build_pose(sess: Session, atom_array, no_optH: bool = False):
    return pose_stack_from_atom_array(atom_array, sess.device, context=sess.ctx, no_optH=no_optH)


# ------------------------------------------------------------------ block bookkeeping
def block_map(pose) -> dict:
    """{'res': {resid: block}, 'lig': block or None, 'n_blocks': n}. Protein chain A, ligand chain B."""
    pi = pose.pdb_info
    res, lig = {}, None
    for b in range(pose.max_n_blocks):
        if int(pose.block_type_ind[0, b]) < 0:
            continue
        chain = str(pi.chain_labels[0, b])
        rid = int(pi.residue_labels[0, b])
        if chain == "B":
            lig = b
        else:
            res[rid] = b
    return {"res": res, "lig": lig, "n_blocks": pose.max_n_blocks}


def atoms_of_block(pose, b: int) -> slice:
    off = int(pose.block_coord_offset[0, b])
    return slice(off, off + pose.block_type(0, b).n_atoms)


def atom_index(pose, b: int, name: str) -> int:
    return int(pose.block_coord_offset[0, b]) + pose.block_type(0, b).atom_to_idx[name]


def shell_blocks(pose, lig_b: int, radius: float) -> torch.Tensor:
    """bool [1, max_blocks] - blocks with any atom within `radius` of any ligand atom (ligand included)."""
    xyz = pose.coords[0]
    lig_xyz = xyz[atoms_of_block(pose, lig_b)]
    mask = torch.zeros((1, pose.max_n_blocks), dtype=torch.bool, device=xyz.device)
    for b in range(pose.max_n_blocks):
        if int(pose.block_type_ind[0, b]) < 0:
            continue
        x = xyz[atoms_of_block(pose, b)]
        if (torch.cdist(x, lig_xyz).min() <= radius):
            mask[0, b] = True
    return mask


def coord_mask_for_blocks(pose, block_mask: torch.Tensor, lig_b: int | None, move_ligand: bool) -> torch.Tensor:
    m = torch.zeros(pose.coords.shape[:2], dtype=torch.bool, device=pose.coords.device)
    for b in range(pose.max_n_blocks):
        if int(pose.block_type_ind[0, b]) < 0:
            continue
        if lig_b is not None and b == lig_b:
            m[0, atoms_of_block(pose, b)] = move_ligand
        elif bool(block_mask[0, b]):
            m[0, atoms_of_block(pose, b)] = True
    return m


# ------------------------------------------------------------------ constraints
def with_constraints(pose, fn, atom_idx: torch.Tensor, params: torch.Tensor):
    """Return a pose whose ConstraintSet additionally holds `fn` on `atom_idx` [n,k,3] (pose,block,atom)."""
    from tmol.pose._constraint_set import ConstraintSet

    cs = pose.constraint_set or ConstraintSet.create_empty(pose.device, pose.n_poses)
    cs = cs.add_constraints(fn, atom_idx.to(torch.int32), params.to(torch.float32))
    return attrs.evolve(pose, constraint_set=cs)


def add_ca_restraints(pose):
    return constrain_all_ca(pose)


def strip_constraints(pose):
    return attrs.evolve(pose, constraint_set=None)


# ------------------------------------------------------------------ relax
def _shell_task_op(shell: torch.Tensor):
    def op(task: PackerTask):
        task.restrict_to_repacking()
        task.disable_packing_by_block_mask(~shell)
        task.add_conformer_sampler(create_dunbrack_sampler_from_database(ParameterDatabase.get_default(), task.device))
        task.add_conformer_sampler(FixedAAChiSampler())
        task.add_conformer_sampler(IncludeCurrentSampler())
    return op


def relax(sess: Session, pose, *, pocket_radius: float | None = None, move_ligand: bool = True,
          repeats: int = 2, keep_constraints: bool = True):
    """FastRelax (Cartesian). Full protein if pocket_radius is None, else only residues within
    `pocket_radius` of the ligand repack/move (rest frozen). Constraints already on `pose` stay at
    constant weight (no ramp) so the result *is held by* them.
    """
    bm = block_map(pose)
    ff = FoldForest.reasonable_fold_forest(pose)
    ops = None
    if pocket_radius is None:
        mm = CartesianMoveMap()
        if not move_ligand and bm["lig"] is not None:
            allm = torch.ones(pose.coords.shape[:2], dtype=torch.bool, device=pose.coords.device)
            allm[0, atoms_of_block(pose, bm["lig"])] = False
            mm = CartesianMoveMap(coord_mask=allm)
    else:
        shell = shell_blocks(pose, bm["lig"], pocket_radius)
        mm = CartesianMoveMap(coord_mask=coord_mask_for_blocks(pose, shell, bm["lig"], move_ligand))
        ops = [_shell_task_op(shell)]
    has_cst = pose.constraint_set is not None and keep_constraints
    try:
        return fast_relax(pose, sess.sfxn, PackerPalette(), mm, ff, task_operations=ops, num_repeats=repeats,
                          ramp_constraints=False if has_cst else None)
    finally:
        # tmol 0.2.0 leaves the constraint weight at 0 after fast_relax (verified: tests/smoke_cst2.py);
        # without this every later energy_report silently shows E_constraint = 0.
        sess.sfxn.set_weight(ScoreType.constraint, sess.cst_weight)


# ------------------------------------------------------------------ energies
def term_names(sess: Session) -> list[str]:
    return [st.name for st in sess.sfxn.all_score_types()]


def block_pair_tensor(sess: Session, pose) -> torch.Tensor:
    """[terms, B, B] weighted energies, upper triangle (i<j) = pair, diagonal = one-body/intra."""
    scorer = sess.sfxn.render_block_pair_scoring_module(pose)
    return scorer(pose.coords, sum_terms=False)[:, 0]


def energy_report(sess: Session, pose, pocket_resids: list[int] | None = None) -> dict:
    """Totals, ligand-protein interface energy per term, repulsion split, per-residue energies."""
    bm = block_map(pose)
    lb = bm["lig"]
    M = block_pair_tensor(sess, pose).detach().cpu()
    tn = term_names(sess)
    ti = {n: i for i, n in enumerate(tn)}
    nb = M.shape[-1]
    valid = torch.tensor([int(pose.block_type_ind[0, b]) >= 0 for b in range(nb)])
    Mv = M * (valid[:, None] & valid[None, :]).float()[None]
    total = float(Mv.sum())
    cst = float(Mv[ti["constraint"]].sum()) if "constraint" in ti else 0.0
    out: dict = {"E_total": total, "E_constraint": cst, "E_total_noncst": total - cst}
    n_prot = len(bm["res"])
    out["E_total_per_res"] = (total - cst) / max(n_prot, 1)
    # per-residue energy (one-body + half of every pair)
    S = Mv.sum(0)  # [B,B] all terms
    dg = torch.diagonal(S)
    per_block = dg + 0.5 * (S.sum(1) + S.sum(0) - 2 * dg)
    inv = {b: r for r, b in bm["res"].items()}
    pr = {inv[b]: float(per_block[b]) for b in inv}
    out["per_res_energy"] = pr
    vals = np.array(list(pr.values()))
    out["res_E_max"] = float(vals.max())
    out["n_res_E_gt5"] = float((vals > 5.0).sum())
    if pocket_resids:
        pv = np.array([pr[r] for r in pocket_resids if r in pr])
        out["pocket_E_per_res"] = float(pv.mean()) if len(pv) else float("nan")
        out["pocket_E_max"] = float(pv.max()) if len(pv) else float("nan")
    # repulsion split (fa_ljrep): intra-residue = diagonal, inter = upper triangle
    rep = Mv[ti["fa_ljrep"]]
    diag = torch.diagonal(rep)
    prot_blocks = torch.tensor([b in inv for b in range(nb)])
    out["rep_intra_total"] = float(diag[prot_blocks].sum())
    out["rep_inter_total"] = float(rep.sum() - diag.sum())
    if pocket_resids:
        pb = torch.tensor([bm["res"][r] for r in pocket_resids if r in bm["res"]])
        out["rep_intra_pocket"] = float(diag[pb].sum())
        sub = rep[pb][:, pb]
        out["rep_inter_pocket"] = float(sub.sum() - torch.diagonal(sub).sum())
    if lb is not None:
        cross = Mv[:, :lb, lb]  # [terms, protein blocks]; ligand is last block so pairs are (i<lb, lb)
        for n, i in ti.items():
            v = float(cross[i].sum())
            if abs(v) > 1e-9:
                out[f"int.{n}"] = v
        out["int.total"] = float(cross.sum())
        out["int.rep_over_atr"] = float(cross[ti["fa_ljrep"]].sum() / max(-float(cross[ti["fa_ljatr"]].sum()), 1e-6))
        out["E_ligand_intra"] = float(Mv[:, lb, lb].sum())
        # restraint energy (CA + theozyme) is bookkeeping, not strain: exclude it
        out["E_protein_in_complex"] = total - cst - float(cross.sum()) - float(Mv[:, lb, lb].sum())
        # per-residue interaction with the ligand (for hot-spot tables)
        per_res_int = cross.sum(0)
        out["int_per_res"] = {inv[b]: float(per_res_int[b]) for b in inv}
        out["int_hbond_per_res"] = {inv[b]: float(cross[ti["hbond"]][b]) for b in inv if abs(float(cross[ti["hbond"]][b])) > 1e-4}
    return out


# ------------------------------------------------------------------ pose <-> atoms
def pose_to_replicate_arrays(pose):
    """(protein AtomArray chain A incl. H, ligand AtomArray chain B incl. H)."""
    aa = atom_array_from_pose_stack(pose)
    aa = aa[0] if hasattr(aa, "stack_depth") else aa
    return aa[aa.chain_id == "A"], aa[aa.chain_id == "B"]


def pose_charges(sess: Session, pose):
    """(xyz [N,3], q [N], block_of_atom [N]) with Rosetta/MMFF partial charges used by fa_elec."""
    pbt = pose.packed_block_types
    sess.sfxn.render_whole_pose_scoring_module(pose)  # ensures elec annotations exist
    q_tab = pbt.elec_partial_charge.cpu().numpy()  # [n_types, max_atoms]
    xyz, q, blk = [], [], []
    for b in range(pose.max_n_blocks):
        t = int(pose.block_type_ind[0, b])
        if t < 0:
            continue
        sl = atoms_of_block(pose, b)
        n = sl.stop - sl.start
        xyz.append(pose.coords[0, sl].detach().cpu().numpy())
        q.append(q_tab[t, :n])
        blk.append(np.full(n, b))
    return np.vstack(xyz), np.concatenate(q), np.concatenate(blk)
