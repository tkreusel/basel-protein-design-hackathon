"""Build the results folder from the collected tables: tables, sequences, and one folder of structures per selection.

usage: export_results.py [--camp NAME] [--top 50]       writes to $KEMP_RESULTS (re-runnable)

Selections (folders with sequences.fasta, summary.csv, metrics_all_columns.csv and one subfolder per sequence):
  passing_strict_boltz_gate      every evaluated sequence that passed the strict Boltz-1x gate
  tsa_preferred_over_substrate   z_tsa_vs_sub > 2
  shortlist_all_criteria         strict gate + PLACER signal + ipSAE2(TSA) >= 0.7 + TSA preferred over substrate + PyRosetta dG_sep <= -15
  top50_by_panel_rank            descriptive panel rank (report.py)
  final_top25_ipsae0.8_composite / final_strict_all_gates_ipsae0.8   from final_scoring.py, if present
Per-sequence subfolder: design_RFD4.cif, boltz_TSA_6NT_best.cif, boltz_substrate_H5J_best.cif, boltz_product_6VP_best.cif,
placer_holo/apo_ensemble20.pdb, sequence.fasta, metrics.csv.
"""
import argparse
import glob
import os
import shutil

import numpy as np
import pandas as pd

from . import config, flags
from . import report as report_mod

SUMMARY_COLUMNS = ["final_rank", "composite_z", "composite_rank", "seq_id", "motif", "arm", "variant", "eval_set", "boltz_tier", "z_tsa_vs_sub", "z_tsa_vs_prod",
                   "ipsae.tsa", "ipsae.sub", "ipsae.prod", "boltz_ipde.tsa", "placer.apo.frac_pass_core", "placer.holo.frac_pass_core", "ros.tsa.dG_sep",
                   "delta.tsa_minus_sub.total", "6NT.pocket_rmsd.mean", "ligand_rmsd_ca_aligned", "motif_tip_rmsd_ca_aligned"]


def all_sequences(root):
    seqs = {}
    for p in glob.glob(f"{root}/*/*/b*/mpnn/sequences.csv"):
        s = pd.read_csv(p)
        seqs.update(dict(zip(s.seq_id, s.sequence)))
    return seqs


def export_group(out, root, seqs, name, df):
    """Write one selection folder. Structures that no longer exist on disk are skipped."""
    group = f"{out}/{name}"
    shutil.rmtree(group, ignore_errors=True)
    os.makedirs(group)
    with open(group + "/sequences.fasta", "w") as fh:
        for x in df.itertuples():
            fh.write(f">{x.seq_id} {x.motif} {x.arm} {x.variant}\n{seqs.get(x.seq_id, '')}\n")
    df.drop(columns=["family"], errors="ignore").to_csv(group + "/metrics_all_columns.csv", index=False)
    df[[k for k in SUMMARY_COLUMNS if k in df]].round(3).to_csv(group + "/summary.csv", index=False)
    for i, (_, x) in enumerate(df.iterrows(), 1):
        run = f"{root}/{x.motif}/{x.arm}/{x.batch}"
        forest, ev_dir = (f"{run}/eval_forest_x", f"{run}/eval_x") if x.eval_set == "extended" else (f"{run}/eval_forest", f"{run}/eval")
        d = f"{group}/{i:02d}_{x.motif}_{x.arm}_{x.seq_id[-34:]}"
        os.makedirs(d, exist_ok=True)
        files = ((f"{run}/gen/designs/{x.design}.cif", "design_RFD4.cif"),
                 (f"{forest}/{x.seq_id}__6NT/{x.seq_id}__6NT_model_0.cif", "boltz_TSA_6NT_best.cif"),
                 (f"{forest}/{x.seq_id}__H5J/{x.seq_id}__H5J_model_0.cif", "boltz_substrate_H5J_best.cif"),
                 (f"{forest}/{x.seq_id}__6VP/{x.seq_id}__6VP_model_0.cif", "boltz_product_6VP_best.cif"))
        for src, dst in files:
            if os.path.exists(src):
                shutil.copy(src, f"{d}/{dst}")
        for mode in ("holo", "apo"):
            for p in glob.glob(f"{ev_dir}/work/placer/{x.seq_id}_boltz/*_tsa_boltz_m0_{mode}_model.pdb"):
                shutil.copy(p, f"{d}/placer_{mode}_ensemble20.pdb")
        open(d + "/sequence.fasta", "w").write(f">{x.seq_id}\n{seqs.get(x.seq_id, '')}\n")
        x.dropna().to_csv(d + "/metrics.csv")
    return len(df)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--camp")
    ap.add_argument("--top", type=int, default=50)
    a = ap.parse_args()
    root, out = config.runs_root(a.camp), str(config.results_dir())
    collected = root / "_collect"
    for d in ("tables", "sequences", "report"):
        os.makedirs(f"{out}/{d}", exist_ok=True)
    report_mod.build_report(a.camp, 40)       # refreshes _collect/ (funnel, candidates, REPORT.md, top_candidates)
    for f in ("funnel.csv", "arm_summary.csv", "candidates.csv", "top_candidates.csv"):
        shutil.copy(collected / f, f"{out}/tables/{f}")
    shutil.copy(collected / "REPORT.md", out + "/report/REPORT.md")
    gates = [pd.read_csv(p).drop(columns=["cif"], errors="ignore") for p in sorted(glob.glob(f"{root}/*/*/b*/gate.csv"))]
    designs = pd.concat(gates, ignore_index=True)
    designs.to_csv(out + "/tables/all_designs_metrics.csv", index=False)

    cands = pd.read_csv(collected / "candidates.csv")
    ev = flags.add_flags(cands[(cands.eval_set != "none") & cands["ipsae.tsa"].notna()])
    ev.to_csv(out + "/tables/evaluated_sequences_full_suite.csv", index=False)
    seqs = all_sequences(root)

    sel = {"passing_strict_boltz_gate": ev[ev.flag_strict_boltz].sort_values(["motif", "arm", "z_tsa_vs_sub"], ascending=[True, True, False]),
           "tsa_preferred_over_substrate": ev[ev.flag_tsa_over_sub].sort_values("z_tsa_vs_sub", ascending=False),
           "shortlist_all_criteria": ev[ev.flag_shortlist].sort_values("ipsae.tsa", ascending=False)}
    top = pd.read_csv(collected / "top_candidates.csv").head(a.top)
    sel[f"top{a.top}_by_panel_rank"] = ev.set_index("seq_id").loc[top.seq_id].reset_index()
    final = f"{out}/tables/final_scoring.csv"
    if os.path.exists(final):
        f = pd.read_csv(final)
        f = f[f.composite_z.notna()]
        sel["final_top25_ipsae0.8_composite"] = f.head(25)
        tier = designs.drop_duplicates("design").set_index("design")["gate_tier"]
        strict = f.assign(design_gate_tier=f.design.map(tier))
        strict = strict[(strict.boltz_tier == 1) & (strict.design_gate_tier == 1)].copy()
        strict["final_rank_strict"] = np.arange(1, len(strict) + 1)
        strict.to_csv(out + "/tables/final_strict_all_gates_ipsae0.8.csv", index=False)
        sel["final_strict_all_gates_ipsae0.8"] = strict
    counts = {name: export_group(out, root, seqs, name, df) for name, df in sel.items()}
    with open(out + "/sequences/evaluated_sequences.fasta", "w") as fh:
        for x in ev.itertuples():
            fh.write(f">{x.seq_id} {x.motif} {x.arm} {x.variant} tier{x.boltz_tier}\n{seqs.get(x.seq_id, '')}\n")
    print("results:", out, "| evaluated:", len(ev), "|", counts)


if __name__ == "__main__":
    main()
