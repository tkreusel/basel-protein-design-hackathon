"""Merge tmol / PLACER / ipSAE2 results into one per-entry table + a readable report.

  python -m maxi_eval.aggregate --work work --out work/report

Inputs (any may be missing; absent stages just leave their columns out):
  work/tmol/*.json            per-replicate tmol protocol (run_tmol_entry)
  work/placer_summary.csv     per-replicate PLACER ensemble statistics (analyze_placer)
  work/ipsae/ipsae_by_entry_arm.csv   (run_ipsae)
  config/entries.yaml         labels (kcat_KM)
Outputs: summary_replicates.csv, summary_entries.csv, report.md, and rank correlations vs log10(kcat/KM) when >= 4
distinct labelled entries exist (otherwise it says so instead of printing a meaningless number).
"""
from __future__ import annotations
import argparse, glob, json, math
from pathlib import Path
import numpy as np, pandas as pd, yaml

HEADLINE = [
    # (column, label, what 'good' looks like)
    ("placer.holo.frac_pass_core", "PLACER holo: frac ensemble in catalytic geometry", "high"),
    ("placer.apo.frac_pass_core", "PLACER apo: frac ensemble pre-formed", "high"),
    ("placer.holo.lig_rmsd_spread", "PLACER holo: ligand pose spread (A)", "low"),
    ("placer.holo.chi_sd_base", "PLACER: base chi circular SD (deg)", "low"),
    ("placer.holo.prmsd_mean", "PLACER: mean prmsd (A)", "low"),
    ("rlx.cst.n_pass", "tmol rlx: theozyme constraints satisfied", "high"),
    ("rlx.cst.rms_z", "tmol rlx: rms deviation (tol units)", "low"),
    ("strain.E_noncst", "tmol: strain to impose theozyme (REU)", "low"),
    ("rlx.geo.d_O_X", "base O...X (A)", "ref ~2.5"),
    ("rlx.geo.ang_O_H_X", "O...H-X angle (deg)", "ref ~165"),
    ("rlx.geo.tors_syn", "carboxylate O'-C-O...H (0 syn / 180 anti)", "see report"),
    ("rlx.geo.donor_best_d", "oxyanion donor distance (A)", "ref ~3.0"),
    ("delta.tsa_minus_sub.total", "ddE TSA-substrate (REU)", "negative"),
    ("delta.prod_minus_sub.total", "dE product-substrate (REU)", "positive (no product inhibition)"),
    ("rlx.rep_intra_pocket", "pocket intra-residue repulsion", "low"),
    ("rlx.E_total_per_res", "tmol E/res (REU)", "low"),
    ("rlx.buns.n_unsat", "buried unsatisfied polars (approx.)", "low"),
    ("rlx.void.volume_A3", "pocket void volume (A^3, proxy)", "low"),
    ("rlx.field.CH.env_nobase.vac", "field on C-H, env w/o base (MV/cm)", "positive"),
    ("rlx.field.LG.env.vac", "field on leaving group (MV/cm)", "positive"),
    ("ros.tsa.dG_sep", "Rosetta dG_separated TSA (REU)", "negative"),
    ("ros.delta.tsa_minus_sub.dG_sep", "Rosetta ddG TSA-substrate (REU)", "negative"),
    ("ros.delta.prod_minus_sub.dG_sep", "Rosetta dG product-substrate (REU)", "positive"),
    ("ros.tsa.sc", "Rosetta shape complementarity (ligand)", "high"),
    ("ros.tsa.packstat_int", "Rosetta interface packstat", "high"),
    ("ros.tsa.filter.buns_all", "Rosetta BUNS (pocket)", "low"),
    ("delta_native.tsa_minus_sub.total", "tmol ddE TSA-substrate, independent refolds (REU)", "negative"),
    ("delta_native.prod_minus_sub.total", "tmol dE product-substrate, independent refolds", "positive"),
    ("ros.delta_native.tsa_minus_sub.dG_sep", "Rosetta ddG TSA-substrate, independent refolds", "negative"),
    ("ipsae.prod", "ipSAE2 protein->product (6VP)", "low (product should release)"),
    ("ipsae.tsa", "ipSAE2 protein->TSA", "high"),
    ("ipsae.sub", "ipSAE2 protein->substrate", "high"),
    ("ipsae.delta_tsa_sub", "ipSAE2 TSA - substrate", "positive"),
    ("site.pae_cat_to_lig.tsa", "Boltz PAE catalytic residues->ligand (A)", "low"),
    ("boltz_ipde.tsa", "Boltz ipDE TSA, mean over samples (A)", "low"), ("boltz_ipde.sub", "Boltz ipDE substrate (A)", "low"), ("boltz_ipde.prod", "Boltz ipDE product (A)", "ref"),
]


def load_tmol(work: Path) -> pd.DataFrame:
    rows = []
    for f in sorted(glob.glob(str(work / "tmol" / "*.json"))):
        d = json.load(open(f))["scalars"]
        rows.append({k: v for k, v in d.items() if not isinstance(v, (list, dict))})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default="work")
    ap.add_argument("--out", default="work/report")
    ap.add_argument("--config", default="config/entries.yaml")
    a = ap.parse_args()
    work, out = Path(a.work), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(Path(a.config).read_text())
    labels = {k: v.get("kcat_KM") for k, v in {**cfg["entries"], **{k: e for k, e in (cfg.get("extra_entries") or {}).items() if e.get("enabled")}}.items()}

    rep = load_tmol(work)
    if rep.empty:
        print("no tmol results yet"); rep = pd.DataFrame(columns=["entry", "rep_id"])
    ok = rep[rep.get("status", "ok") == "ok"] if "status" in rep else rep
    if (rep.get("status") == "failed").any():
        print("FAILED replicates:\n", rep[rep.status == "failed"][["rep_id", "error"]].to_string())
    ros = []
    for f in sorted(glob.glob(str(work / "pyrosetta" / "*.json"))):
        d = json.load(open(f))
        if d.get("status") == "ok":
            ros.append({k: v for k, v in d.items() if not isinstance(v, (list, dict)) and not k.startswith("ros.error")})
        else:
            print("PyRosetta FAILED:", d.get("rep_id"), d.get("error"))
    if ros:
        r = pd.DataFrame(ros).drop(columns=["entry", "status"], errors="ignore")
        ok = ok.merge(r, on="rep_id", how="outer") if len(ok) else r
    ah = work / "apo_holo.csv"
    if ah.exists():
        a_ = pd.read_csv(ah).drop(columns=["entry"], errors="ignore")
        ok = ok.merge(a_, on="rep_id", how="outer") if len(ok) else a_
    pl = work / "placer_summary.csv"
    if pl.exists():
        p = pd.read_csv(pl).drop(columns=["entry"], errors="ignore")
        # per-arm PLACER rows (placer_sub.*, placer_prod.*) fold onto the entry's TSA row for the same model
        p["rep_id"] = p.rep_id.str.replace(":sub:", ":tsa:", regex=False).str.replace(":prod:", ":tsa:", regex=False).str.replace(":apo:", ":tsa:", regex=False)
        p = p.drop(columns=["arm"], errors="ignore").groupby("rep_id", as_index=False).first()
        ok = ok.merge(p, on="rep_id", how="outer") if len(ok) else p
        ok["entry"] = ok["entry"].fillna(ok["rep_id"].str.split(":").str[0])
    ok.to_csv(out / "summary_replicates.csv", index=False)

    num = ok.select_dtypes("number")
    ent = num.groupby(ok["entry"]).mean()
    sd = num.groupby(ok["entry"]).std().add_suffix("|sd")
    ent = ent.join(sd)
    ip = work / "ipsae" / "ipsae_by_entry_arm.csv"
    if ip.exists():
        d = pd.read_csv(ip)
        for arm in ("tsa", "sub", "prod"):
            s = d[d.arm == arm].set_index("entry")
            for col, name in (("ipsae.ipsae|mean", f"ipsae.{arm}"), ("ipsae.ipsae|std", f"ipsae.{arm}|sd"),
                              ("site.pae_cat_to_lig|mean", f"site.pae_cat_to_lig.{arm}"), ("site.plddt_cat|mean", f"site.plddt_cat.{arm}"),
                              ("boltz.complex_ipde|mean", f"boltz_ipde.{arm}"), ("boltz.complex_ipde|std", f"boltz_ipde.{arm}|sd"), ("boltz.complex_ipde|min", f"boltz_ipde.{arm}|min"), ("boltz.complex_ipde|max", f"boltz_ipde.{arm}|max"),
                              ("boltz.ligand_iptm|mean", f"boltz_ligand_iptm.{arm}"), ("boltz.ligand_iptm|std", f"boltz_ligand_iptm.{arm}|sd"), ("boltz.complex_iplddt|mean", f"boltz_iplddt.{arm}")):
                if col in s:
                    ent[name] = s[col]
        if "ipsae.tsa" in ent and "ipsae.sub" in ent:
            ent["ipsae.delta_tsa_sub"] = ent["ipsae.tsa"] - ent["ipsae.sub"]
    from .labels import load_labels, label_for
    lab = load_labels(a.config)
    ent["kcat_KM"] = pd.Series({e: (label_for(e, lab) or {}).get("kcat_KM") for e in ent.index}, dtype=float)
    ent["class"] = pd.Series({e: (label_for(e, lab) or {}).get("class") for e in ent.index})
    ent["log10_kcat_KM"] = np.log10(ent["kcat_KM"].astype(float))
    ent.to_csv(out / "summary_entries.csv")

    L = ["# Kemp eliminase in-silico evaluation - summary\n", f"{len(ok)} replicate(s), {len(ent)} entries.\n",
         "Heuristic features for RANKING; none is calibrated to predict kcat/KM (see PIPELINE.md).\n", "## Headline table (entry means)\n"]
    cols = [(c, l, g) for c, l, g in HEADLINE if c in ent.columns]
    if cols:
        t = ent[[c for c, _, _ in cols] + ["kcat_KM"]].copy()
        t.columns = [l for _, l, _ in cols] + ["kcat/KM"]
        try:
            L.append(t.round(2).to_markdown())
        except ImportError:   # tabulate not installed
            L.append(t.round(2).to_string())
        L.append("\n\nDirection of 'good': " + "; ".join(f"{l}: {g}" for _, l, g in cols))
    lab = ent.dropna(subset=["log10_kcat_KM"])
    L.append("\n\n## Rank correlation with log10(kcat/KM)\n")
    if lab["log10_kcat_KM"].nunique() < 4:
        L.append(f"Not computed: only {lab['log10_kcat_KM'].nunique()} distinct labelled activity value(s) (>= 4 needed for a meaningful rank test; "
                 "enable the validation ladder in config/entries.yaml and fill labels).")
    else:
        from scipy.stats import spearmanr
        rows = []
        for c, l, g in cols:
            v = lab[[c, "log10_kcat_KM"]].dropna()
            if len(v) >= 4 and v[c].nunique() > 1:
                r, pv = spearmanr(v[c], v["log10_kcat_KM"]); rows.append((l, len(v), r, pv))
        L.append(pd.DataFrame(rows, columns=["metric", "n", "spearman", "p"]).round(3).to_string(index=False))
    (out / "report.md").write_text("\n".join(L) + "\n")
    print("\n".join(L[:6])); print("->", out)

if __name__ == "__main__":
    main()
