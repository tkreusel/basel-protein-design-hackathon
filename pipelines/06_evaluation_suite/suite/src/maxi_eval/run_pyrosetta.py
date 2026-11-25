"""PyRosetta stage: the Rosetta filters tmol does not have, plus an independent ligand ddG as a cross-check of tmol.

Runs in /mnt/work/envs/pyrosetta (CPU only; PyRosetta has no GPU path). Reads the same manifest as every other stage.

Per TSA replicate (all keys prefixed `ros.`):
  relax       Rosetta FastRelax (beta_nov16), coordinate-restrained to the start, ligand free
  interface   InterfaceAnalyzer (jump 1 = ligand): dG_separated (repack-separated), delta SASA, packstat,
              delta unsatisfied H-bonds, shape complementarity (Sc)
  filters     BuriedUnsatHbonds (pocket + ligand), PackStat, Sc.  NOT Holes: it needs the external DAlphaBall binary,
              which is not in the PyRosetta wheel, and a missing binary kills the whole process (stability.pocket_void_volume stays as the proxy).
  ddG         substrate / product superposed on the relaxed TSA pose, pocket-relaxed: dG_sep(TSA) - dG_sep(substrate) etc.
              (independent of tmol; agreement between the two is the sanity check, disagreement is informative)
  cart_ddg    (opt.) X->Ala of base/donor/wall, cartesian pocket relax, mean of `nstruct` runs:
              ddG_fold-like (E_mut - E_wt) and the change in ligand dG. APPROXIMATION of the cartesian_ddg app.
  enzdes      (opt., --cst) match constraints added from your .cst via REMARK 666 templates; constraint energies reported.
Failures are recorded per metric (`ros.error.<name>`), never swallowed.
"""
from __future__ import annotations
import argparse, json, os, sys, time, traceback, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")

LIG_RESNAME = {"tsa": "6NT", "sub": "H5J", "prod": "6VP"}


def init(ligand_dir: Path, extra: str = ""):
    import pyrosetta
    params = " ".join(str(ligand_dir / f"{k}.params") for k in ("6NT", "H5J", "6VP"))
    pyrosetta.init(f"-mute all -corrections::beta_nov16 true -ex1 -ex2 -use_input_sc -extra_res_fa {params} -no_optH false "
                   f"-in:ignore_unrecognized_res false -relax:constrain_relax_to_start_coords {extra}", silent=True)
    return pyrosetta


def read_mol2(path):
    names, xyz, on = [], [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            on = line.startswith("@<TRIPOS>ATOM"); continue
        if on and line.strip():
            f = line.split(); names.append(f[1]); xyz.append([float(f[2]), float(f[3]), float(f[4])])
    return names, np.array(xyz)


def pocket_resids(pose, lig, radius=8.0):
    import pyrosetta.rosetta as R
    lx = [pose.residue(lig).xyz(i) for i in range(1, pose.residue(lig).natoms() + 1)]
    out = []
    for r in range(1, pose.total_residue()):
        res = pose.residue(r)
        if any((res.xyz(i) - x).norm() <= radius for i in range(1, res.natoms() + 1) for x in lx):
            out.append(r)
    return out


def xml_obj(xml: str):
    from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
    return XmlObjects.create_from_string(xml)


def relax(pose, sfxn, repeats=2, cartesian=False, movable=None, ligand_free=True):
    import pyrosetta.rosetta as R
    from pyrosetta.rosetta.protocols.relax import FastRelax
    mm = R.core.kinematics.MoveMap()
    n = pose.total_residue()
    for r in range(1, n + 1):
        mv = movable is None or r in movable
        mm.set_bb(r, bool(mv) and r != n)
        mm.set_chi(r, bool(mv))
    mm.set_jump(1, ligand_free)
    # FastRelax(constrain_relax_to_start_coords) leaves coordinate_constraint weight 0.5 on the score function it is given
    # (verified: +5214 REU in later scores). Always hand it a clone so the caller's scorefxn stays clean.
    fr = FastRelax(sfxn.clone(), repeats)
    fr.constrain_relax_to_start_coords(True)
    fr.set_movemap(mm)
    if cartesian:
        fr.cartesian(True)
        fr.min_type("lbfgs_armijo_nonmonotone")
    fr.apply(pose)
    return pose


def sc_ligand(pose, lig: int) -> float:
    """Lawrence-Colman Sc between protein and the ligand.

    Rosetta's ShapeComplementarityCalculator/filter silently ignores non-polymer residues (AddResidue adds 0 atoms and
    the filter returns -1/0), so ligand heavy atoms are added by hand with element vdW radii.
    """
    import pyrosetta.rosetta as R
    A = R.core.scoring.sc.Atom
    calc = R.core.scoring.sc.ShapeComplementarityCalculator(); calc.Init()
    for i in range(1, lig):
        calc.AddResidue(0, pose.residue(i))
    rsd = pose.residue(lig)
    rad = {"C": 1.7, "N": 1.55, "O": 1.52}
    n = 0
    for i in range(1, rsd.natoms() + 1):
        if rsd.atom_is_hydrogen(i):
            continue
        at = A(); x = rsd.xyz(i)
        at.x, at.y, at.z = x.x, x.y, x.z
        at.radius = rad.get(rsd.atom_type(i).element(), 1.7)
        n += calc.AddAtom(1, at)
    if n == 0 or calc.Calc() != 1:
        raise RuntimeError("sc calculation failed")
    return float(calc.GetResults().sc)


def interface(pose, sfxn):
    import pyrosetta.rosetta as R
    sfxn(pose)  # IAM reads cached energies: without this the cross-term energy is 0.0 for freshly modified poses
    iam = R.protocols.analysis.InterfaceAnalyzerMover(1, False, sfxn, True, True, True, False)
    iam.apply(pose)
    D = iam.get_all_data()

    def f(v):  # some InterfaceData fields are vector1_double (index 1 = total)
        try:
            return float(v)
        except TypeError:
            return float(list(v)[0])
    d = {k: f(v) for k, v in {"dG_sep": D.dG, "dSASA": D.dSASA,
         "packstat_int": D.packstat, "delta_unsat_hb": D.delta_unsat_hbonds, "n_hb": D.interface_hbonds,
         "hbond_E": D.total_hb_E, "n_int_res": D.interface_nres}.items()}
    try:
        d["sc"] = sc_ligand(pose, pose.total_residue())       # Sc is known to be unreliable: its failure must not discard dG_sep / packstat / hbonds
    except Exception:
        d["sc"] = float("nan")
    return d


def filters(pose, lig, pocket):
    out, err = {}, {}
    res_str = ",".join(str(r) for r in pocket)
    xml = f"""<ROSETTASCRIPTS>
  <SCOREFXNS><ScoreFunction name="sfxn" weights="beta_nov16"/></SCOREFXNS>
  <RESIDUE_SELECTORS><Index name="pocket" resnums="{res_str}"/><Index name="ligsel" resnums="{lig}"/></RESIDUE_SELECTORS>
  <FILTERS>
    <PackStat name="packstat" threshold="0.0" chain="1" repeats="1" confidence="0"/>
    <BuriedUnsatHbonds name="buns_all" scorefxn="sfxn" jump_number="0" cutoff="100" residue_selector="pocket" confidence="0"/>
    <BuriedUnsatHbonds name="buns_lig" scorefxn="sfxn" jump_number="0" cutoff="100" residue_selector="ligsel" confidence="0"/>
  </FILTERS>
  <PROTOCOLS/>
</ROSETTASCRIPTS>"""
    xo = xml_obj(xml)
    for name in ("packstat", "buns_all", "buns_lig"):
        try:
            out[name] = float(xo.get_filter(name).report_sm(pose))
        except Exception as exc:  # noqa: BLE001
            err[name] = f"{type(exc).__name__}: {str(exc)[:120]}"
    return out, err


def place_ligand(pose, lig, target: str, ligand_dir: Path):
    """Replace the TSA residue by substrate/product superposed on it (mapped atoms from ligands.CORRESPONDENCE_TO_TSA)."""
    import pyrosetta.rosetta as R
    from maxi_eval import ligands as L
    ts = R.core.chemical.ChemicalManager.get_instance().residue_type_set("fa_standard")
    rsd_tsa = pose.residue(lig)
    tsa_xyz = {rsd_tsa.atom_name(i).strip(): np.array(list(rsd_tsa.xyz(i))) for i in range(1, rsd_tsa.natoms() + 1)}
    names, xyz = read_mol2(ligand_dir / f"{target}.mol2")
    corr = L._best_corr(L.CORRESPONDENCE_TO_TSA[target], dict(zip(names, xyz)), tsa_xyz)
    src = np.array([xyz[names.index(k)] for k in corr]); dst = np.array([tsa_xyz[v] for v in corr.values()])
    Rm, t = L._kabsch(src, dst)
    new = R.core.conformation.ResidueFactory.create_residue(ts.name_map(target))
    placed = xyz @ Rm.T + t
    for n, p in zip(names, placed):
        if new.has(n):
            new.set_xyz(n, R.numeric.xyzVector_double_t(*p))
    pose.replace_residue(lig, new, False)
    return L._kabsch_rmsd(src, dst)


def enzdes_remarks(pdb_text_lines, cst_blocks, rep_site):
    """REMARK 666 lines (one per CST block) so enzdes knows which residue plays which block."""
    out = []
    for i, b in enumerate(cst_blocks, start=1):
        pick = rep_site["base"] if i == 1 else None  # first block = base unless the CST says otherwise; see PIPELINE.md
        if pick:
            out.append(f"REMARK 666 MATCH TEMPLATE B 6NT    1 MATCH MOTIF A {pick['resname']} {pick['resid']:>4}  {i}  1")
    return out


def run_replicate(r: dict, e: dict, a, pr) -> dict:
    import pyrosetta.rosetta as R
    t0 = time.time()
    res: dict = {"rep_id": r["rep_id"], "entry": e["_entry"], "status": "ok"}
    err: dict = {}
    sfxn = pr.create_score_function("beta_nov16")
    pose = pr.pose_from_file(r["pdb"])
    lig = pose.total_residue()
    assert pose.residue(lig).name3().strip() == "6NT", pose.residue(lig).name3()
    site = r.get("site") or {}
    base, donors, wall = site.get("base"), site.get("donors") or [], site.get("wall")

    # raw
    res["ros.raw.score_per_res"] = sfxn(pose) / (lig - 1)
    if a.cst:
        pass  # handled below, after relax
    # relax (restrained, ligand free)
    if a.relax_repeats > 0:
        relax(pose, sfxn, a.relax_repeats)
    res["ros.rlx.score_per_res"] = sfxn(pose) / (lig - 1)
    pocket = pocket_resids(pose, lig, 8.0)
    res["ros.n_pocket_res"] = len(pocket)

    def block(tag, p):
        try:
            res.update({f"ros.{tag}.{k}": v for k, v in interface(p, sfxn).items()})
        except Exception as exc:  # noqa: BLE001
            err[f"{tag}.interface"] = f"{type(exc).__name__}: {str(exc)[:150]}"
    block("tsa", pose)
    f, fe = filters(pose, lig, pocket)
    res.update({f"ros.tsa.filter.{k}": v for k, v in f.items()}); err.update({f"filter.{k}": v for k, v in fe.items()})

    # same-pocket substrate / product (independent ddG vs tmol)
    for arm, tgt in (("sub", "H5J"), ("prod", "6VP")):
        try:
            p2 = pose.clone()
            res[f"ros.{arm}.place_rmsd"] = place_ligand(p2, lig, tgt, Path(a.ligand_dir))
            relax(p2, sfxn, 1, movable=set(pocket) | {lig})
            block(arm, p2)
        except Exception as exc:  # noqa: BLE001
            err[f"{arm}.swap"] = f"{type(exc).__name__}: {str(exc)[:150]}"
    for k in ("dG_sep", "sc", "packstat_int", "delta_unsat_hb", "hbond_E"):
        t_, s_, p_ = (res.get(f"ros.{x}.{k}") for x in ("tsa", "sub", "prod"))
        if t_ is not None and s_ is not None:
            res[f"ros.delta.tsa_minus_sub.{k}"] = t_ - s_
        if p_ is not None and s_ is not None:
            res[f"ros.delta.prod_minus_sub.{k}"] = p_ - s_

    # independent native refolds of the substrate (H5J) and product (6VP): same full relax + interface analysis as the TSA
    if r["source"] == "boltz":
        for x in e["replicates"]:
            if x["arm"] in ("sub", "prod") and x["source"] == "boltz" and x.get("boltz_model") == r.get("boltz_model"):
                try:
                    pn = pr.pose_from_file(x["pdb"])
                    assert pn.residue(pn.total_residue()).name3().strip() == LIG_RESNAME[x["arm"]], pn.residue(pn.total_residue()).name3()
                    if a.relax_repeats > 0:
                        relax(pn, sfxn, a.relax_repeats)
                    res.update({f"ros.native_{x['arm']}.{k}": v for k, v in interface(pn, sfxn).items()})
                    res[f"ros.native_{x['arm']}.score_per_res"] = sfxn(pn) / (pn.total_residue() - 1)
                except Exception as exc:  # noqa: BLE001
                    err[f"native_{x['arm']}"] = f"{type(exc).__name__}: {str(exc)[:150]}"
        for k in ("dG_sep", "sc", "packstat_int", "hbond_E"):
            t_, s_, p_ = res.get(f"ros.tsa.{k}"), res.get(f"ros.native_sub.{k}"), res.get(f"ros.native_prod.{k}")
            if t_ is not None and s_ is not None:
                res[f"ros.delta_native.tsa_minus_sub.{k}"] = t_ - s_
            if p_ is not None and s_ is not None:
                res[f"ros.delta_native.prod_minus_sub.{k}"] = p_ - s_

    # enzdes constraints (only with a CST file)
    if a.cst:
        try:
            from maxi_eval import cst as C
            blocks = C.parse_enzdes_cst(Path(a.cst).read_text())
            lines = [l for l in Path(r["pdb"]).read_text().splitlines() if not l.startswith("REMARK 666")]
            rem = enzdes_remarks(lines, blocks, {"base": base})
            tmp = Path(a.tmpdir) / (r["rep_id"].replace(":", "_") + "_enz.pdb")
            tmp.write_text("\n".join(rem + lines) + "\n")
            pe = pr.pose_from_file(str(tmp))
            from pyrosetta.rosetta.protocols.enzdes import AddOrRemoveMatchCsts
            pr.rosetta.basic.options.set_file_option("enzdes:cstfile", str(a.cst))
            mv = AddOrRemoveMatchCsts(); mv.set_cst_action(R.protocols.enzdes.CstAction.ADD_NEW); mv.apply(pe)
            sf2 = sfxn.clone()
            for st in (R.core.scoring.atom_pair_constraint, R.core.scoring.angle_constraint, R.core.scoring.dihedral_constraint):
                sf2.set_weight(st, 1.0)
            sf2(pe)
            en = pe.energies().total_energies()
            for st, nm in ((R.core.scoring.atom_pair_constraint, "atom_pair"), (R.core.scoring.angle_constraint, "angle"),
                           (R.core.scoring.dihedral_constraint, "dihedral")):
                res[f"ros.enzdes.cst_{nm}"] = float(en[st])
        except Exception as exc:  # noqa: BLE001
            err["enzdes"] = f"{type(exc).__name__}: {str(exc)[:200]}"

    # cartesian ddG-like ablations
    if a.ablate:
        scart = pr.create_score_function("beta_nov16_cart")
        movable = set(pocket) | {lig}
        roles = {"base": base["resid"] if base else None,
                 "donor": next((d["resid"] for d in donors if not d.get("backbone")), None),
                 "wall": wall["resid"] if wall else None}
        # pose numbering == PDB numbering after pose_from_file? use pdb_info to map
        pinfo = pose.pdb_info()
        for role, rid in roles.items():
            if rid is None:
                continue
            try:
                seq = pinfo.pdb2pose("A", int(rid))
                es_wt, es_mut, dints = [], [], []
                for k in range(a.nstruct):
                    pw = pose.clone(); relax(pw, scart, 1, cartesian=True, movable=movable); es_wt.append(scart(pw))
                    pm = pose.clone()
                    R.protocols.simple_moves.MutateResidue(seq, "ALA").apply(pm)
                    relax(pm, scart, 1, cartesian=True, movable=movable); es_mut.append(scart(pm))
                    dints.append(interface(pm, sfxn)["dG_sep"] - interface(pw, sfxn)["dG_sep"])
                res[f"ros.cartddg.{role}{rid}.ddG_fold_like"] = float(np.mean(es_mut) - np.mean(es_wt))
                res[f"ros.cartddg.{role}{rid}.d_dG_ligand"] = float(np.mean(dints))
            except Exception as exc:  # noqa: BLE001
                err[f"cartddg.{role}"] = f"{type(exc).__name__}: {str(exc)[:150]}"
    for k, v in err.items():
        res[f"ros.error.{k}"] = v
    res["ros.time_s"] = time.time() - t0
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs")
    ap.add_argument("--out", default="work/pyrosetta")
    ap.add_argument("--ligand-dir", default="ligands")
    ap.add_argument("--entries", default=None)
    ap.add_argument("--max-replicates", type=int, default=None)
    ap.add_argument("--relax-repeats", type=int, default=2)
    ap.add_argument("--reverse", action="store_true", help="walk entries from the END; lock files stop several runners colliding")
    ap.add_argument("--ablate", action="store_true")
    ap.add_argument("--nstruct", type=int, default=2)
    ap.add_argument("--cst", default=None)
    ap.add_argument("--tmpdir", default=os.environ.get("TMPDIR", "/tmp"))
    a = ap.parse_args()
    pr = init(Path(a.ligand_dir))
    man = json.loads((Path(a.inputs) / "manifest.json").read_text())
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    items = list(man["entries"].items())
    if a.reverse:
        items.reverse()
    for entry, e in items:
        if a.entries and entry not in a.entries.split(","):
            continue
        e["_entry"] = entry
        reps = [r for r in e["replicates"] if r["arm"] == "tsa"]
        if a.max_replicates:
            reps = reps[: a.max_replicates]
        for r in reps:
            jf = out / (r["rep_id"].replace(":", "_") + ".json")
            lock = jf.with_suffix(".lock")
            if jf.exists() or lock.exists():
                print("skip (done/locked)", r["rep_id"]); continue
            lock.write_text("running")
            try:
                res = run_replicate(r, e, a, pr)
            except Exception as exc:  # fail loud, continue
                res = {"rep_id": r["rep_id"], "entry": entry, "status": "failed", "error": f"{type(exc).__name__}: {exc}",
                       "traceback": traceback.format_exc()}
                print("FAILED", r["rep_id"], exc, file=sys.stderr)
            jf.write_text(json.dumps(res, indent=1, default=float))
            lock.unlink(missing_ok=True)
            print(f"done {r['rep_id']} {res.get('ros.time_s', 0):.0f}s status={res['status']}", flush=True)

if __name__ == "__main__":
    main()
