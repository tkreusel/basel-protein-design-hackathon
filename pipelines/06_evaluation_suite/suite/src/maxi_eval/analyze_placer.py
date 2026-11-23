"""PLACER ensemble -> preorganisation statistics (runs in the RFD4/tmol venv; reads run_placer.py output).

Preorganisation = how much of the PLACER ensemble sits in the catalytic geometry, and how tightly.
  holo  ensemble (ligand + side chains predicted)
     frac_pass_core / frac_pass_all   fraction of models satisfying the theozyme constraints (core = distances + linear H-bond)
     frac_donor                       fraction with the oxyanion donor within 3.5 A of the leaving-group atom
     lig_rmsd_spread                  mean pairwise ligand heavy-atom RMSD across models   (small = defined pose)
     chi_sd_{base,donor,wall}         circular SD of chi angles of the catalytic residues  (small = rigid)
     syn_frac, mean/sd of key distances/angles
     prmsd_mean / top10 / plddt_mean  PLACER's own confidence
  apo   ensemble (ligand removed, side chains predicted around the base): the *reference* TSA pose from the input is
        transplanted into each model, so frac_pass_* = "fraction of the unbound pocket already pre-formed".
  field (optional, --field-models N): tmol charges -> field on the C-H and leaving-group axes per model; the SD across the
        ensemble is the quantitative 'electrostatic preorganisation'.
"""
from __future__ import annotations
import argparse, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from . import cst as C, geometry as G, ligands as L, structures as S

CORE = ["geo.d_O_X", "geo.d_O_H", "geo.ang_O_H_X", "geo.ang_C_O_H"]


def load_models(pdb: Path):
    """Read models one by one: PLACER models can differ in atom count (stacked reading would fail)."""
    import biotite.structure.io.pdb as pdb_io
    f = pdb_io.PDBFile.read(str(pdb))
    return [f.get_structure(model=i, altloc="first", include_bonds=False) for i in range(1, f.get_model_count() + 1)]


def ligand_rmsd_matrix(ligs: list[np.ndarray]) -> np.ndarray:
    X = np.stack(ligs)
    n = len(X)
    D = np.zeros((n, n))
    for i in range(n):
        D[i] = np.sqrt(((X - X[i]) ** 2).sum(-1).mean(-1))
    return D


def summarize(models, ref_rep, site, constraints, csv: Path, mode: str, field_ctx=None, n_field=0) -> tuple[dict, pd.DataFrame]:
    sc = pd.read_csv(csv)
    rows, ligs = [], []
    chis: dict[str, list] = {"base": [], "donor": [], "wall": []}
    roles_res = {"base": site["base"]["resid"] if site.get("base") else None,
                 "donor": next((d["resid"] for d in site["donors"] if not d["backbone"]), None),
                 "wall": site["wall"]["resid"] if site.get("wall") else None}
    for i, m in enumerate(models):
        prot = m[(m.chain_id == "A") & (m.element != "H")]
        lig = m[(m.chain_id == "B") & (m.element != "H")]
        if mode == "apo" or len(lig) == 0:
            lig = ref_rep.ligand[ref_rep.ligand.element != "H"]
        rep = S.Replicate(f"m{i+1}", ref_rep.entry, ref_rep.arm, "placer", prot, lig)
        met, _ = G.all_geometry(rep, site=site)
        row = {"model": i + 1, **{k: v for k, v in met.items() if not k.startswith("geo.wsite")}}
        ev = C.evaluate(met, constraints)
        row.update({k: v for k, v in ev.items() if k.endswith(".pass") or k in ("cst.n_pass", "cst.frac_pass", "cst.rms_z")})
        core_keys = [c.short for c in constraints if c.key in CORE]
        has_base = "geo.d_O_X" in met      # the product (6VP) has no abstractable proton -> no catalytic-geometry pass/fail
        row["pass_core"] = float(all(ev.get(f"cst.{k}.pass", 0.0) == 1.0 for k in core_keys)) if (core_keys and has_base) else float("nan")
        row["pass_all"] = float(ev["cst.frac_pass"] == 1.0) if has_base else float("nan")
        row["donor_3p5"] = float(met.get("geo.donor_best_d", 99) <= 3.5) if site["donors"] else float("nan")
        for role, rid in roles_res.items():
            if rid is not None:
                chis[role].append(G.chi_angles(rep, rid))
        ligs.append(lig.coord)
        rows.append(row)
    df = pd.DataFrame(rows)
    df = df.merge(sc[["model_idx", "prmsd", "plddt", "plddt_pde"]].rename(columns={"model_idx": "model"}), on="model", how="left")
    out = {"n_models": float(len(df)), "frac_pass_core": float(df.pass_core.mean()), "frac_pass_all": float(df.pass_all.mean()),
           "frac_donor": float(df.donor_3p5.mean()) if df.donor_3p5.notna().any() else float("nan"),
           "syn_frac": float(df["geo.syn_flag"].mean()) if "geo.syn_flag" in df else float("nan"),
           "prmsd_mean": float(df.prmsd.mean()), "prmsd_top10": float(df.prmsd.nsmallest(max(1, len(df) // 10)).mean()),
           "plddt_mean": float(df.plddt.mean())}
    for k in ("geo.d_O_X", "geo.ang_O_H_X", "geo.ang_C_O_H"):
        if k in df:
            out[f"{k[4:]}_mean"], out[f"{k[4:]}_sd"] = float(df[k].mean()), float(df[k].std())
    if "geo.tors_syn" in df:
        out["tors_syn_circ_sd"] = G.circ_std_deg(df["geo.tors_syn"])
    if mode == "holo" and len(ligs) > 1:
        D = ligand_rmsd_matrix(ligs)
        out["lig_rmsd_spread"] = float(D[np.triu_indices(len(D), 1)].mean())
        ref = ref_rep.ligand[ref_rep.ligand.element != "H"].coord
        out["lig_rmsd_to_input"] = float(np.mean([np.sqrt(((x - ref) ** 2).sum(-1).mean()) for x in ligs]))
    for role, lst in chis.items():
        if lst:
            sds = [G.circ_std_deg([c[k] for c in lst if k in c]) for k in lst[0] if sum(k in c for c in lst) == len(lst)]
            if sds:
                out[f"chi_sd_{role}"] = float(np.mean(sds))
    # optional tmol field per model
    if field_ctx is not None and n_field > 0:
        fl = field_per_model(field_ctx, models, ref_rep, site, mode, n_field, df)
        for k in fl.columns:
            df.loc[fl.index, k] = fl[k]
        for k in fl.columns:
            out[f"{k}_mean"], out[f"{k}_sd"] = float(fl[k].mean()), float(fl[k].std())
    return out, df


def field_per_model(ctx, models, ref_rep, site, mode, n, df) -> pd.DataFrame:
    from . import field as F, tmol_core as TC
    sess, templates = ctx
    order = df.sort_values("prmsd").index[:n]
    res = {}
    for i in order:
        m = models[i]
        prot = m[(m.chain_id == "A") & (m.element != "H")]
        lig = m[(m.chain_id == "B") & (m.element != "H")]
        if mode == "apo" or len(lig) == 0:
            lig = ref_rep.ligand[ref_rep.ligand.element != "H"]
        try:
            pose = TC.build_pose(sess, prot + lig)
            bm = TC.block_map(pose)
            xyz, q, blk = TC.pose_charges(sess, pose)
            p2, l2 = TC.pose_to_replicate_arrays(pose)
            ly = G.lig_xyz(S.Replicate("x", "x", ref_rep.arm, "p", p2, l2))
            f = F.project(xyz, q, blk, ly, L.ROLES[ref_rep.lig_resname], bm["lig"], bm["res"].get(site["base"]["resid"]) if site.get("base") else None, modes=("vac",))
            res[i] = {"field_CH": f["field.CH.env_nobase.vac"], "field_LG": f["field.LG.env.vac"]}
        except Exception as exc:  # noqa: BLE001
            print("  field failed for model", i, exc)
    return pd.DataFrame.from_dict(res, orient="index")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs")
    ap.add_argument("--placer", default="work/placer")
    ap.add_argument("--out", default="work/placer_summary.csv")
    ap.add_argument("--reference", default="config/theozyme_reference.yaml")
    ap.add_argument("--group", default="hg3_lineage")
    ap.add_argument("--ligand-dir", default="ligands")
    ap.add_argument("--field-models", type=int, default=0, help="per replicate/mode, top-N models by prmsd get a tmol field calculation")
    a = ap.parse_args()
    T = L.build_templates(Path(a.ligand_dir))
    constraints = C.load_reference(Path(a.reference), a.group)
    man = json.loads((Path(a.inputs) / "manifest.json").read_text())
    field_ctx = None
    summary = []
    for entry, e in man["entries"].items():
        for r in e["replicates"]:
            base = Path(a.placer) / entry / r["rep_id"].replace(":", "_")
            if not Path(f"{base}_holo.csv").exists() and not Path(f"{base}_apo.csv").exists():
                continue
            if r["arm"] == "apo":      # ligand-free input: transplant the entry's TSA pose into the apo frame as the reference
                from .apo_holo import transplant_ligand
                trep = next((x for x in e["replicates"] if x["arm"] == "tsa"), None)
                if trep is None:
                    continue
                holo = S.extract_replicates(Path(trep["pdb"]), entry, "tsa", trep["source"], T)[0]
                apo_p = S.extract_apo(Path(r["pdb"]), entry).protein
                lig, _, _, _ = transplant_ligand(apo_p, holo.protein, holo.ligand)
                ref_rep = S.Replicate(r["rep_id"], entry, "tsa", r["source"], apo_p, lig)
                site = G.detect_site(holo, e.get("catalytic_overrides"))
            else:
                ref_rep = S.extract_replicates(Path(r["pdb"]), entry, r["arm"], r["source"], T)[0]
                site = G.detect_site(ref_rep, e.get("catalytic_overrides"))
            row = {"rep_id": r["rep_id"], "entry": entry, "arm": r["arm"]}
            pre = "placer" if r["arm"] == "tsa" else f"placer_{r['arm']}arm" if r["arm"] == "apo" else f"placer_{r['arm']}"   # tsa keeps the original column names
            for mode in ("holo", "apo"):
                csv, pdb = Path(f"{base}_{mode}.csv"), Path(f"{base}_{mode}_model.pdb")
                if not csv.exists():
                    continue
                if a.field_models and field_ctx is None:
                    from . import tmol_core as TC
                    field_ctx = (TC.Session.create(Path(a.ligand_dir), ref_rep.protein, T), T)
                try:
                    out, df = summarize(load_models(pdb), ref_rep, site, constraints, csv, mode, field_ctx, a.field_models)
                except Exception as exc:     # one bad model must not cost the whole run its PLACER summary
                    print("ANALYSIS FAILED", r["rep_id"], mode, f"{type(exc).__name__}: {exc}"[:200], flush=True); continue
                df.to_csv(f"{base}_{mode}_models.csv", index=False)
                row.update({f"{pre}.{mode}.{k}": v for k, v in out.items()})
            summary.append(row)
            print("analysed", r["rep_id"])
    pd.DataFrame(summary).to_csv(a.out, index=False)
    print("->", a.out)

if __name__ == "__main__":
    main()
