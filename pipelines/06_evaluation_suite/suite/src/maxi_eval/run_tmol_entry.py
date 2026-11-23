"""CLI: tmol protocol on all TSA replicates of the given entries. Resumable (one JSON per replicate).

  python -m maxi_eval.run_tmol_entry --inputs work/inputs --out work/tmol [--entries 7K4Q,4BS0] [--apo] [--ablate]
"""
from __future__ import annotations
import argparse, json, sys, time, traceback, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import torch
from . import ligands as L, structures as S, tmol_core as TC, protocol as P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs")
    ap.add_argument("--out", default="work/tmol")
    ap.add_argument("--entries", default=None, help="comma list; default all in manifest")
    ap.add_argument("--max-replicates", type=int, default=None, help="cap replicates per entry (e.g. 3 of 10 NMR models)")
    ap.add_argument("--ligand-dir", default="ligands")
    ap.add_argument("--reference", default="config/theozyme_reference.yaml")
    ap.add_argument("--group", default="hg3_lineage")
    ap.add_argument("--cst", default=None, help="Rosetta enzdes .cst file (optional)")
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--cst-scale", type=float, default=10.0)
    ap.add_argument("--reverse", action="store_true", help="walk the entry list from the END (a second runner meets the forward runners in the middle; finished/locked replicates are skipped)")
    ap.add_argument("--apo", action="store_true")
    ap.add_argument("--ablate", action="store_true")
    a = ap.parse_args()

    inputs, out = Path(a.inputs), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    man = json.loads((inputs / "manifest.json").read_text())
    T = L.build_templates(Path(a.ligand_dir))
    cfg = P.Config(Path(a.ligand_dir), Path(a.reference), a.group, a.repeats, 8.0, a.cst_scale, a.apo, a.ablate, cst_file=a.cst)
    first = True
    items = list(man["entries"].items())
    if a.reverse:
        items.reverse()
    for entry, e in items:
        if a.entries and entry not in a.entries.split(","):
            continue
        reps = [r for r in e["replicates"] if r["arm"] == "tsa"]
        if a.max_replicates:
            reps = reps[: a.max_replicates]
        for r in reps:
            jf = out / (r["rep_id"].replace(":", "_") + ".json")
            lock = jf.with_suffix(".lock")
            if jf.exists() or lock.exists():
                print("skip (done/locked)", r["rep_id"]); continue
            lock.write_text("running")
            rep = S.extract_replicates(Path(r["pdb"]), entry, "tsa", r["source"], T)
            if not rep:
                print("!! no replicate parsed from", r["pdb"]); continue
            rep = rep[0]; rep.rep_id = r["rep_id"]
            if first:
                sess = TC.Session.create(Path(a.ligand_dir), rep.protein, T, cst_weight=1.0)
                ev = P.Evaluator(sess, T, cfg); first = False
            t = time.time()
            try:
                res = ev.run_tsa_replicate(rep, e.get("catalytic_overrides"))
                sc = res["scalars"]
                if r["source"] == "boltz":      # independent substrate / product refolds of the same model index
                    from . import geometry as G
                    site = G.detect_site(rep, e.get("catalytic_overrides"))
                    for arm in ("sub", "prod"):
                        cand = [x for x in e["replicates"] if x["arm"] == arm and x.get("boltz_model") == r.get("boltz_model")] \
                               or [x for x in e["replicates"] if x["arm"] == arm]
                        if not cand:
                            continue
                        arep = S.extract_replicates(Path(cand[0]["pdb"]), entry, arm, "boltz", T)[0]
                        sc.update(ev.run_native_arm(arep, site))
                    for term in ("total", "hbond", "fa_elec", "fa_lk", "fa_ljatr", "fa_ljrep", "lk_bridge", "lk_ball"):
                        tv = sc.get(f"rlx.int.{term}")
                        if tv is not None and f"native_sub.int.{term}" in sc:
                            sc[f"delta_native.tsa_minus_sub.{term}"] = tv - sc[f"native_sub.int.{term}"]
                        if tv is not None and f"native_prod.int.{term}" in sc and f"native_sub.int.{term}" in sc:
                            sc[f"delta_native.prod_minus_sub.{term}"] = sc[f"native_prod.int.{term}"] - sc[f"native_sub.int.{term}"]
                sc["status"] = "ok"
            except Exception as exc:  # fail loud per replicate, continue the batch
                res = {"scalars": {"rep_id": r["rep_id"], "entry": entry, "status": "failed", "error": f"{type(exc).__name__}: {exc}"},
                       "detail": {"traceback": traceback.format_exc()}}
                print("FAILED", r["rep_id"], exc, file=sys.stderr)
            jf.write_text(json.dumps(res, indent=1, default=float))
            lock.unlink(missing_ok=True)
            print(f"done {r['rep_id']} {time.time()-t:.0f}s status={res['scalars']['status']}", flush=True)

if __name__ == "__main__":
    main()
