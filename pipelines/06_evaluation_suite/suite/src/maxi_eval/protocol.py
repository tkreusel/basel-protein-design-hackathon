"""The tmol evaluation protocol for one TSA-bound replicate (see PIPELINE.md section 4).

    raw        as given (tmol protonation only)
    rlx        full Cartesian FastRelax, CA restrained (0.5 A), ligand free      -> "does the site hold by itself?"
    rlxcst     same + flat-bottomed theozyme constraints held (weight `cst_scale`)
               strain = E_noncst(rlxcst) - E_noncst(rlx)                          -> cost of imposing the theozyme
    tsa_pk     pocket-only relax (8 A) of the rlx complex  \\
    sub        substrate superposed on the TSA, pocket relax, ligand free          } same protein, same protocol,
    prod       product superposed on the TSA, pocket relax, ligand free            } so differences are the ligand
    sub_ts / prod_ts   ligand held at the TSA pose, only side chains repack       -> charge/shape at the TS geometry
    delta.*    ddE_tsa_sub, ddE_prod_sub, ddE_prod_tsa (+ per-term)
    apo        (optional) protein-only relax -> reorganisation energy
    abl.*      (optional) X->Ala of base / donor / wall: change in stability and in ligand interaction

HONEST LIMITS (also in PIPELINE.md): 6NT is a neutral TSA, not the delocalised-anion TS; Rosetta energies are not
free energies; relax is a local minimisation, so "sub free" is "substrate refined near the TS pose", not a
Michaelis-complex search (use the Boltz substrate arm for an independent pose). Use deltas to RANK, not to predict kcat/KM.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import torch

from . import cst as C
from . import field as F
from . import geometry as G
from . import ligands as L
from . import stability as ST
from . import structures as S
from . import tmol_core as TC


@dataclass
class Config:
    ligand_dir: Path
    reference_yaml: Path
    reference_group: str = "hg3_lineage"
    relax_repeats: int = 2
    pocket_radius: float = 8.0
    cst_scale: float = 10.0
    do_apo: bool = False
    do_ablation: bool = False
    do_field: bool = True
    do_stability: bool = True
    cst_file: Path | None = None
    cst_alias: dict = field(default_factory=dict)


def _scalars(prefix: str, d: dict) -> dict:
    return {f"{prefix}.{k}": v for k, v in d.items() if isinstance(v, (int, float, np.floating)) and not isinstance(v, bool)}


def swap_complex(prot_noH, tsa_lig, target: str, templates):
    """Protein (heavy) + target ligand (H5J/6VP) superposed on the TSA pose."""
    tpl = templates[target]
    tsa_xyz = {str(n): c for n, c in zip(tsa_lig.atom_name, tsa_lig.coord)}
    names, xyz, rmsd = L.place_onto_tsa(tpl, tsa_xyz)
    heavy = [i for i, n in enumerate(names) if not n.startswith("H")]
    els = {n: a.GetSymbol() for n, a in zip(tpl.names, tpl.mol.GetAtoms())}
    arr = S.make_ligand_array([names[i] for i in heavy], xyz[heavy], target, [els[names[i]] for i in heavy])
    return prot_noH + arr, rmsd


def _strip_h(a):
    return a[a.element != "H"]


class Evaluator:
    def __init__(self, sess: TC.Session, templates, cfg: Config):
        self.sess, self.T, self.cfg = sess, templates, cfg
        self.constraints = C.load_reference(cfg.reference_yaml, cfg.reference_group)
        self.cst_blocks = None
        if cfg.cst_file:
            self.cst_blocks = C.parse_enzdes_cst(Path(cfg.cst_file).read_text())

    # ---------------------------------------------------------------- one pose -> metrics
    def eval_pose(self, pose, site, tag, arm="tsa", detail=None, geometry=True) -> dict:
        sess = self.sess
        rn = L.ARM_TO_RESNAME[arm]
        site_a = {**site, "roles": L.ROLES[rn]}
        en = TC.energy_report(sess, pose, site["pocket_resids"])
        out = _scalars(tag, en)
        if detail is not None:
            detail[tag] = {k: en[k] for k in ("per_res_energy", "int_per_res", "int_hbond_per_res") if k in en}
        prot, lig = TC.pose_to_replicate_arrays(pose)
        rep = S.Replicate(f"{tag}", "x", arm, "tmol", prot, lig)
        if geometry:
            if arm == "prod":      # no abstraction proton / base geometry for the product
                m = {**G.oxyanion_metrics(rep, site_a), **G.burial_metrics(rep)}
            else:
                m, _ = G.all_geometry(rep, site=site_a)
            out.update({f"{tag}.{k}": v for k, v in m.items()})
            if arm in ("tsa", "sub"):
                out.update({f"{tag}.{k}": v for k, v in C.evaluate(m, self.constraints).items()})
                if self.cst_blocks is not None:
                    out.update({f"{tag}.{k}": v for k, v in C.evaluate_enzdes_cst(rep, self.cst_blocks, self.cfg.cst_alias).items()})
        bm = TC.block_map(pose)
        if self.cfg.do_field and arm in ("tsa", "sub") and site.get("base"):
            xyz, q, blk = TC.pose_charges(sess, pose)
            ly = G.lig_xyz(rep)
            out.update({f"{tag}.{k}": v for k, v in F.project(xyz, q, blk, ly, L.ROLES[rn], bm["lig"], bm["res"].get(site["base"]["resid"])).items()})
        if self.cfg.do_stability:
            b = ST.buried_unsatisfied_polars(prot, lig, site["pocket_resids"])
            out.update({f"{tag}.{k}": v for k, v in b.items() if k != "buns.atoms"})
            if detail is not None:
                detail.setdefault(tag, {})["buns_atoms"] = b["buns.atoms"]
            out.update({f"{tag}.{k}": v for k, v in ST.pocket_void_volume(prot, lig).items()})
        return out

    # ---------------------------------------------------------------- the protocol
    def run_tsa_replicate(self, rep: S.Replicate, overrides: dict | None = None) -> dict:
        sess, cfg = self.sess, self.cfg
        t0 = time.time()
        site = G.detect_site(rep, overrides)
        res = {"rep_id": rep.rep_id, "entry": rep.entry, "source": rep.source, "n_res": int(len(set(rep.protein.res_id)))}
        detail: dict = {"site": {k: v for k, v in site.items() if k != "roles"}}
        aa = _strip_h(rep.protein) + _strip_h(rep.ligand)
        pose0 = TC.build_pose(sess, aa)
        bm = TC.block_map(pose0)
        res.update(self.eval_pose(pose0, site, "raw", detail=detail))

        # --- rlx: held only by CA restraints
        p_ca = TC.add_ca_restraints(pose0)
        pose1 = TC.relax(sess, p_ca, repeats=cfg.relax_repeats)
        res.update(self.eval_pose(pose1, site, "rlx", detail=detail))
        res["time_rlx_s"] = time.time() - t0

        # --- rlxcst: theozyme held
        pose_c, applied, skipped = C.add_to_pose(p_ca, site, self.constraints, bm, scale=cfg.cst_scale)
        res["cst.applied"], res["cst.skipped"] = ";".join(applied), ";".join(skipped)
        pose2 = TC.relax(sess, pose_c, repeats=cfg.relax_repeats)
        res.update(self.eval_pose(pose2, site, "rlxcst", detail=detail))
        res["strain.E_noncst"] = res["rlxcst.E_total_noncst"] - res["rlx.E_total_noncst"]
        res["strain.int_total"] = res["rlxcst.int.total"] - res["rlx.int.total"]

        # --- same-pocket ligand comparison from the relaxed (free) TSA complex
        prot1, lig1 = TC.pose_to_replicate_arrays(pose1)
        prot1h, lig1h = _strip_h(prot1), _strip_h(lig1)
        arms = {}
        arms["tsa_pk"] = ("tsa", prot1h + lig1h, True)
        for arm, tgt in (("sub", "H5J"), ("prod", "6VP")):
            cplx, rmsd = swap_complex(prot1h, lig1h, tgt, self.T)
            res[f"{arm}.place_rmsd"] = rmsd
            arms[arm] = (arm, cplx, True)
            arms[f"{arm}_ts"] = (arm, cplx, False)
        for tag, (arm, cplx, move) in arms.items():
            p = TC.build_pose(sess, cplx)
            p = TC.relax(sess, p, pocket_radius=cfg.pocket_radius, move_ligand=move, repeats=1)
            res.update(self.eval_pose(p, site, tag, arm=arm, detail=detail, geometry=move or tag == "tsa_pk"))
        for term in ("total", "hbond", "fa_elec", "fa_lk", "fa_ljatr", "fa_ljrep", "lk_bridge", "lk_ball"):
            k = f"int.{term}"
            g = lambda t: res.get(f"{t}.{k}", 0.0)  # noqa: E731
            res[f"delta.tsa_minus_sub.{term}"] = g("tsa_pk") - g("sub")
            res[f"delta.prod_minus_sub.{term}"] = g("prod") - g("sub")
            res[f"delta.prod_minus_tsa.{term}"] = g("prod") - g("tsa_pk")
            res[f"delta.prod_ts_minus_sub.{term}"] = g("prod_ts") - g("sub")
            res[f"delta.tsa_minus_sub_ts.{term}"] = g("tsa_pk") - g("sub_ts")

        # --- optional: apo reorganisation
        if cfg.do_apo:
            apo = TC.build_pose(sess, _strip_h(prot1))
            apo = TC.relax(sess, TC.add_ca_restraints(apo), repeats=cfg.relax_repeats)
            e_apo = TC.energy_report(sess, apo)["E_total_noncst"]
            res["apo.E_total_noncst"] = e_apo
            res["apo.reorg_E"] = res["rlx.E_protein_in_complex"] - e_apo  # includes CA-restraint-free energy of holo protein

        # --- optional: ablation scan (X -> Ala)
        if cfg.do_ablation:
            res.update(self.ablate(pose1, site, res))
        res["time_total_s"] = time.time() - t0
        return {"scalars": res, "detail": detail}

    def run_native_arm(self, rep: S.Replicate, site: dict) -> dict:
        """Independent Boltz refold of the substrate/product arm: same full CA-restrained relax as `rlx`."""
        aa = _strip_h(rep.protein) + _strip_h(rep.ligand)
        pose = TC.relax(self.sess, TC.add_ca_restraints(TC.build_pose(self.sess, aa)), repeats=self.cfg.relax_repeats)
        return self.eval_pose(pose, site, f"native_{rep.arm}", arm=rep.arm, geometry=True)

    # ---------------------------------------------------------------- ablations
    def ablate(self, pose1, site, res) -> dict:
        out = {}
        bm = TC.block_map(pose1)
        roles = {"base": site["base"]["resid"] if site.get("base") else None,
                 "donor": next((d["resid"] for d in site["donors"] if not d["backbone"]), None),
                 "wall": site["wall"]["resid"] if site.get("wall") else None}
        for role, rid in roles.items():
            if rid is None or rid not in bm["res"]:
                continue
            p = mutate_pack(self.sess, pose1, bm["res"][rid], "ALA", self.cfg.pocket_radius)
            p = TC.relax(self.sess, p, pocket_radius=self.cfg.pocket_radius, repeats=1)
            en = TC.energy_report(self.sess, p, site["pocket_resids"])
            out[f"abl.{role}{rid}.dE_total"] = en["E_total_noncst"] - res["rlx.E_total_noncst"]
            out[f"abl.{role}{rid}.d_int"] = en["int.total"] - res["rlx.int.total"]
            out[f"abl.{role}{rid}.d_hbond"] = en.get("int.hbond", 0.0) - res.get("rlx.int.hbond", 0.0)
        return out


def mutate_pack(sess: TC.Session, pose, block: int, name3: str, radius: float):
    """Mutate one block to `name3` and repack its neighbourhood (adapted from the RFD4-Proteina tmol ddG helpers)."""
    from tmol.pack import PackerPalette, PackerTask, pack_rotamers
    from tmol.pack.rotamer import FixedAAChiSampler, IncludeCurrentSampler
    from tmol.pack.rotamer.dunbrack import create_dunbrack_sampler_from_database
    from tmol.database import ParameterDatabase

    bm = TC.block_map(pose)
    shell = TC.shell_blocks(pose, bm["lig"], radius)
    mut = torch.zeros_like(shell)
    mut[0, block] = True
    task = PackerTask(pose, PackerPalette())
    task.restrict_to_repacking()
    considered = task.per_block_considered_block_types
    p_i, b_i, slot = (mut.unsqueeze(-1) & (considered != -1)).nonzero(as_tuple=True)
    match = torch.tensor([bt.name3 == name3 for bt in task.pbt.active_block_types], device=task.device)
    task.per_block_is_block_type_allowed[p_i, b_i, slot] = match[considered[p_i, b_i, slot]]
    task.disable_packing_by_block_mask(~(shell | mut))
    task.add_conformer_sampler(create_dunbrack_sampler_from_database(ParameterDatabase.get_default(), task.device))
    task.add_conformer_sampler(FixedAAChiSampler())
    task.add_conformer_sampler(IncludeCurrentSampler())
    return pack_rotamers(pose, sess.sfxn, task)
