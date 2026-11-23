"""Collect everything that has finished so far into OUT/{funnel.csv, candidates.csv, arm_summary.csv}.

usage: collect.py [--camp NAME] [--out DIR]       (default OUT = <runs>/_collect)

funnel.csv       one row per run (motif, arm, batch): designs, lineages, gate tiers, sequences, Boltz 1x strict/fallback, stage flags
candidates.csv   one row per sequence that went through the Boltz 1x gate: design metrics, gate flags, stage-5 report (means over the samples,
                 TSA-substrate and TSA-product differences) and, when stage B is done, the evaluation aggregate (ipSAE2, PLACER, tmol, PyRosetta)
arm_summary.csv  per (motif, arm family, sheet): designs, lineages, pass rates, Boltz strict passes, designs with a strict pass
"""
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import config


def _read_ids(path):
    return set(open(path).read().split()) if os.path.exists(path) else set()


def _run_candidates(run):
    """Per-sequence table of one run, or None if the Boltz 1x gate has not run."""
    if not os.path.exists(f"{run}/boltz1_eval.csv"):
        return None
    c = pd.read_csv(f"{run}/boltz1_eval.csv").assign(batch=run.name)
    reports = [pd.read_csv(p).drop(columns=["design", "variant"], errors="ignore")
               for p in (f"{run}/eval_forest/stage5_report.csv", f"{run}/eval_forest_x/stage5_report.csv") if os.path.exists(p)]
    if reports:
        c = c.merge(pd.concat(reports, ignore_index=True), on="seq_id", how="left")
    evals = []
    for p in (f"{run}/eval/work/report/summary_entries.csv", f"{run}/eval_x/work/report/summary_entries.csv"):
        if os.path.exists(p):
            e = pd.read_csv(p)
            e["seq_id"] = e.iloc[:, 0].astype(str).str.replace("_boltz$", "", regex=True)
            evals.append(e.drop(columns=[e.columns[0]]))
    if evals:
        c = c.merge(pd.concat(evals, ignore_index=True), on="seq_id", how="left", suffixes=("", "_eval"))
    main, extended = _read_ids(f"{run}/survivors.txt"), _read_ids(f"{run}/survivors_x.txt")
    c["eval_set"] = np.where(c.seq_id.isin(main), "main", np.where(c.seq_id.isin(extended), "extended", "none"))
    return c


def arm_summary(gates, cands):
    g = gates.copy()
    g["family"], g["sheet"] = g.arm.map(config.arm_family), g.arm.str.endswith("_sheet")
    key = ["motif", "family", "sheet"]
    s = g.groupby(key).agg(designs=("design", "size"), lineages=("lineage", "nunique"), pass_geom=("pass_geom", "mean"), pass_seq=("pass_seq", "mean"),
                           pass_core=("pass_core", "mean"), ch18_ok=("channel18_ok", "mean"), tmol_median=("tmol_ddg", "median"),
                           bottleneck_median=("bottleneck_radius", "median")).round(3)
    if cands is not None and len(cands):
        c = cands.copy()
        c["family"], c["sheet"] = c.arm.map(config.arm_family), c.arm.str.endswith("_sheet")
        s = s.join(c.groupby(key).agg(sequences_boltz1=("seq_id", "size"), boltz_strict=("boltz_gate_strict", "sum"),
                                      designs_with_strict=("design", lambda x: x[c.loc[x.index, "boltz_gate_strict"]].nunique())))
    return s


def collect(camp=None, out=None):
    root = config.runs_root(camp)
    out = Path(out) if out else root / "_collect"
    out.mkdir(parents=True, exist_ok=True)
    funnel, cands, gates = [], [], []
    for run in sorted(root.glob("*/*/b*")):
        motif, arm = run.parent.parent.name, run.parent.name
        f = dict(motif=motif, arm=arm, batch=run.name, stageA=(run / "STAGE_A_DONE").exists(), stageB=(run / "STAGE_B_DONE").exists(),
                 generated=(run / "gen/metrics.csv").exists())
        if (run / "batch_summary.csv").exists():
            f.update(pd.read_csv(run / "batch_summary.csv").iloc[0].drop("run").to_dict())
        funnel.append(f)
        if (run / "gate.csv").exists():
            gates.append(pd.read_csv(run / "gate.csv").assign(batch=run.name))
        c = _run_candidates(run)
        if c is not None:
            cands.append(c)
    funnel = pd.DataFrame(funnel)
    funnel.to_csv(out / "funnel.csv", index=False)
    cands = pd.concat(cands, ignore_index=True) if cands else None
    if cands is not None:
        cands.to_csv(out / "candidates.csv", index=False)
    if gates:
        s = arm_summary(pd.concat(gates, ignore_index=True), cands)
        s.to_csv(out / "arm_summary.csv")
        pd.set_option("display.width", 250)
        print(s.to_string())
    print(len(funnel), "runs;", int(funnel.stageA.sum()), "stage A done;", int(funnel.stageB.sum()), "stage B done ->", out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--camp")
    ap.add_argument("--out")
    a = ap.parse_args()
    collect(a.camp, a.out)


if __name__ == "__main__":
    main()
