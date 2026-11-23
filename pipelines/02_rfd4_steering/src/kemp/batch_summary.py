"""Funnel numbers of one run, written to RUN/batch_summary.csv.   usage: batch_summary.py RUN"""
import sys

import pandas as pd


def summarize(run):
    r = {"run": run}
    try:
        g = pd.read_csv(f"{run}/gate.csv")
        r.update(designs=len(g), lineages=g.lineage.nunique(), pass_geom=int(g.pass_geom.sum()), pass_core=int(g.pass_core.sum()),
                 gate_tier1=int((g.gate_tier == 1).sum()), gate_tier2=int((g.gate_tier == 2).sum()), gate_tier3=int((g.gate_tier == 3).sum()),
                 ch18_ok=int(g.channel18_ok.sum()), tmol_median=float(g.tmol_ddg.median()), strand_mean=float(g.strand.mean()), helix_mean=float(g.helix.mean()))
        r["sequences"] = len(pd.read_csv(f"{run}/mpnn/sequences.csv"))
        b = pd.read_csv(f"{run}/boltz1_eval.csv")
        r.update(boltz1_predictions=len(b), boltz_strict=int(b.boltz_gate_strict.sum()), boltz_fallback=int((b.boltz_tier == 2).sum()),
                 boltz_designs_with_pass=int(b[b.boltz_gate].design.nunique()))
        for c in ("ca_rmsd", "ligand_rmsd_ca_aligned", "ligand_iptm", "complex_ipde"):
            if c in b:
                r["boltz1_median_" + c] = float(b[c].median())
    except Exception as exc:                 # partial runs still get a summary of what exists
        r["error"] = repr(exc)[:100]
    pd.DataFrame([r]).to_csv(f"{run}/batch_summary.csv", index=False)
    return r


if __name__ == "__main__":
    print({k: v for k, v in summarize(sys.argv[1]).items() if k != "run"})
