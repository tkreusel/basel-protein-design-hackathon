"""Status report from the collected tables: OUT/REPORT.md and OUT/top_candidates.csv.

usage: report.py [--camp NAME] [--top 40]

The panel rank is descriptive only. No activity labels exist, so the panel is not calibrated: it is the mean percentile over a
pre-declared set of evaluated metrics (direction in PANEL). TMol is also a steering reward in the steered arms and is weighted like the others.
"""
import argparse

import numpy as np
import pandas as pd

from . import collect as collect_mod
from . import config

PANEL = {"ipsae.tsa": +1, "boltz_ipde.tsa": -1, "6NT.pocket_rmsd.mean": -1, "placer.apo.frac_pass_core": +1, "placer.holo.frac_pass_core": +1,
         "ros.tsa.dG_sep": -1, "ros.delta_native.tsa_minus_sub.dG_sep": -1, "delta.tsa_minus_sub.total": -1}


def panel_rank(c):
    """Evaluated sequences with panel_rank (mean percentile, higher = better on average) and panel_n (metrics available)."""
    e = c[(c.eval_set != "none") & c["ipsae.tsa"].notna()].copy()
    panel = {k: v for k, v in PANEL.items() if k in e}
    ranks = pd.DataFrame({k: e[k].rank(ascending=(v > 0), pct=True) for k, v in panel.items()})
    e["panel_rank"], e["panel_n"] = ranks.mean(axis=1, skipna=True), ranks.notna().sum(axis=1)
    return e.sort_values("panel_rank", ascending=False), panel


def build_report(camp=None, top=40):
    """Refresh the collected tables and write REPORT.md and top_candidates.csv next to them. Returns the output directory."""
    out = collect_mod.collect(camp)
    f = pd.read_csv(out / "funnel.csv")
    c = pd.read_csv(out / "candidates.csv") if (out / "candidates.csv").exists() else pd.DataFrame()
    lines = [f"# Kemp campaign {camp or config.campaign_name()}: status\n", f"{len(f)} run directories; stage A done {int(f.stageA.sum())}; stage B done {int(f.stageB.sum())}\n"]
    done = f[f.stageA == True]
    if len(done):
        cols = [x for x in ("designs", "lineages", "pass_geom", "pass_core", "gate_tier1", "gate_tier2", "gate_tier3", "ch18_ok", "sequences", "boltz_strict", "boltz_fallback") if x in done]
        lines += ["## Funnel per finished run\n", done[["motif", "arm", "batch"] + cols].to_markdown(index=False), "\n"]
        summed = done.groupby(["motif", "arm"])[[x for x in ("designs", "lineages", "pass_geom", "pass_core", "sequences", "boltz_strict") if x in done]].sum()
        lines += ["## Summed per motif and arm\n", summed.to_markdown(), "\n"]
    if len(c):
        c["arm_family"] = c.arm.map(config.arm_family) + np.where(c.arm.str.endswith("_sheet"), " (sheet)", "")
        e, panel = panel_rank(c)
        if len(e):
            keep = ["seq_id", "motif", "arm", "variant", "eval_set", "boltz_tier", "panel_rank", "panel_n", "ca_rmsd", "ligand_rmsd_ca_aligned",
                    "motif_tip_rmsd_ca_aligned"] + list(panel) + ["delta_iptm.TSA_minus_sub", "delta_iptm.TSA_minus_prod"]
            keep = [k for k in keep if k in e]
            e[keep].to_csv(out / "top_candidates.csv", index=False)
            lines += [f"## Top {top} evaluated sequences by panel rank\n", e[keep].head(top).round(2).to_markdown(index=False), "\n"]
            per_arm = e.groupby(["motif", "arm_family"]).agg(n=("seq_id", "size"), strict=("boltz_tier", lambda x: int((x == 1).sum())),
                                                             panel_rank_mean=("panel_rank", "mean"), **{k: (k, "median") for k in panel}).round(2)
            lines += ["## Evaluated sequences per arm (medians)\n", per_arm.to_markdown(), "\n"]
    lines.append("The panel is pre-declared but uncalibrated. Unsteered designs are gated on geometry only; sheet runs are a small diversity add-on; "
                 "the best Boltz sample is the one with the lowest pocket RMSD to the design.\n")
    (out / "REPORT.md").write_text("\n".join(lines))
    print("wrote", out / "REPORT.md")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--camp")
    ap.add_argument("--top", type=int, default=40)
    a = ap.parse_args()
    build_report(a.camp, a.top)


if __name__ == "__main__":
    main()
