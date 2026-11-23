"""Compute the 7 selected ranking scores (S01,S03,S04,S05,S06,S07,S12) for any run of the pipeline.

  python -m maxi_eval.selected_scores --report work/report --extra work/features_extra.csv --out work/selected_scores.csv

Needs: <report>/summary_entries.csv (aggregate) and <extra> (python -m maxi_eval.features_extra --inputs work/inputs --out work/features_extra.csv;
reads work/tmol/*.json + the Boltz PAE npz + PLACER outputs listed in work/inputs/manifest.json).
All seven: LOWER (more negative) = more active in the 19-enzyme reference set (see kemp_selected_scores_list.csv).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

SCORES = {  # name -> (source, a, b)   value = a - b  (b None: a alone);  source 'ent' = summary_entries.csv, 'ex' = features_extra.csv
    "S01_tmol_lkballiso_TSA-minus-PROD": ("ent", "tsa_pk.int.lk_ball_iso", "prod.int.lk_ball_iso"),
    "S03_tmol_faelec_SUB-minus-PROD": ("ent", "sub.int.fa_elec", "prod.int.fa_elec"),
    "S04_tmol_interaction_TSA": ("ent", "tsa_pk.int.total", None),
    "S05_tmol_Etotal_TSA-minus-PROD": ("ent", "tsa_pk.E_total", "prod.E_total"),
    "S06_tmol_unsat_ligand_polars_TSA-minus-PROD": ("ent", "tsa_pk.buns.n_unsat_ligand", "prod.buns.n_unsat_ligand"),
    "S07_tmol_polar_residue_int_TSA-minus-SUB": ("ex", "pr.tsa_pk-sub.int_polar", None),
    "S12_boltz_catalytic_to_pocket_PAE_SUB-minus-PROD": ("ex", "bz.sub-prod.cat_to_pocket_pae", None),
}


def compute(report: Path, extra: Path) -> pd.DataFrame:
    ent = pd.read_csv(report / "summary_entries.csv", index_col=0); ent = ent[ent.index.str.endswith("_boltz")]; ent.index = ent.index.str.replace("_boltz", "", regex=False)
    ex = pd.read_csv(extra, index_col=0); ex.index = ex.index.astype(str).str.replace("_boltz", "", regex=False)
    out = pd.DataFrame(index=ent.index)
    for name, (src, a, b) in SCORES.items():
        T = ent if src == "ent" else ex
        if a not in T.columns or (b and b not in T.columns):
            out[name] = float("nan"); print(f"[warn] {name}: missing column(s) {a} {b or ''}"); continue
        v = T[a].astype(float) - (T[b].astype(float) if b else 0.0)
        out[name] = v.reindex(out.index)
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--report", default="work/report"); ap.add_argument("--extra", default="work/features_extra.csv"); ap.add_argument("--out", default="work/selected_scores.csv")
    a = ap.parse_args(); df = compute(Path(a.report), Path(a.extra)); df.insert(0, "entry", df.index); df.to_csv(a.out, index=False, float_format="%.5g"); print(df.round(3).to_string(index=False)); print("->", a.out)


if __name__ == "__main__":
    main()
