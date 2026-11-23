"""Export the reference dataset with the user-selected scores + plots.

  PYTHONPATH=src python -m maxi_eval.export_selected --pick 1,3,4,5,6,7,12
Outputs (work/report/deliverables/): kemp_reference_dataset_selected.csv, kemp_selected_scores_list.csv, figures/selected_scores_all19.png, figures/selected_scores_HG3.png
"""
from __future__ import annotations
import argparse, shutil
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from .metric_proposal import PICKS, load, values
from .make_deliverables import load_ipae21

SHORT = {1: "S01_tmol_lkballiso_TSA-minus-PROD", 3: "S03_tmol_faelec_SUB-minus-PROD", 4: "S04_tmol_interaction_TSA", 5: "S05_tmol_Etotal_TSA-minus-PROD",
         6: "S06_tmol_unsat_ligand_polars_TSA-minus-PROD", 7: "S07_tmol_polar_residue_int_TSA-minus-SUB", 12: "S12_boltz_catalytic_to_pocket_PAE_SUB-minus-PROD"}
ARM = {"tsa": "TSA", "sub": "SUB", "prod": "PROD"}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pick", default="1,3,4,5,6,7,12"); ap.add_argument("--out", default="work/report/deliverables")
    a = ap.parse_args(); out = Path(a.out); ids = [int(x) for x in a.pick.split(",")]
    M, R, ex, ent = load(); work = Path("work")
    lab = ent[ent.log10_kcat_KM.notna()]
    ds = pd.read_csv("/mnt/work/manually_created_dataset.csv", sep="\t"); ds.columns = [c.strip() for c in ds.columns]
    ds = ds.rename(columns={"PDB ID": "pdb_ids", "Name": "name", "source": "kcat_KM_source"}).set_index("name")
    ipae = load_ipae21(work)
    df = pd.DataFrame(index=lab.index); df.insert(0, "enzyme", lab.index)
    df["group"] = lab.family
    df["pdb_ids_literature"] = ds.pdb_ids.reindex(lab.index); df["kcat_over_KM_M-1s-1"] = lab.kcat_KM; df["log10_kcat_over_KM"] = lab.log10_kcat_KM
    df["kcat_KM_source"] = ds.kcat_KM_source.reindex(lab.index); df["n_mutations_vs_HG3"] = ex.get("cmp.n_mut_vs_HG3")
    for arm, A in ARM.items():
        for key, nm in (("boltz.iptm", "iptm"), ("boltz.ligand_iptm", "ligand_iptm"), ("boltz.complex_iplddt", "complex_iplddt"), ("site.plddt_cat", "plddt_catalytic_residues"), ("site.plddt_lig", "plddt_ligand")):
            df[f"boltz_{nm}_{A}"] = ent[f"ip.{arm}.{key}"]
    for arm, A in ARM.items():
        for key, nm in (("ipsae.p2l", "ipSAE_protein_to_ligand"), ("ipsae.l2p", "ipSAE_ligand_to_protein"), ("ipsae.iPAE12", "iPAE_protein_to_ligand"), ("site.pae_cat_to_lig", "PAE_catalytic_to_ligand"), ("site.pae_lig_to_cat", "PAE_ligand_to_catalytic")):
            df[f"ipsae2_{nm}_{A}"] = ent[f"ip.{arm}.{key}"]
        if not ipae.empty:
            df[f"ipsae2_iPAE_ligand_to_protein_{A}"] = ipae[f"ip.{arm}.ipsae.iPAE21"]
    for arm, A in {**ARM, "apo": "APO"}.items():
        c = f"bz.{arm}.cat_to_pocket_pae"
        if c in ex.columns:
            df[f"boltz_PAE_catalytic_to_pocket_{A}"] = ex[c]
    scores = []
    for i in ids:
        cid, rec, tool, name, expl, cav = PICKS[i - 1][:6]; assert True
        v = values(cid, M, ex, ent).reindex(df.index); df[SHORT[i]] = v
        r = R.loc[cid] if cid in R.index else None
        ok = v.notna(); hg = ok & (df.group == "HG3")
        scores.append(dict(proposal_id=i, score_column_in_this_csv=SHORT[i], pipeline_column_id=cid, tool=tool, meaning=expl,
                           direction="higher value = higher activity" if spearmanr(v[ok], df.log10_kcat_over_KM[ok])[0] > 0 else "lower (more negative) value = higher activity",
                           rho_all19=spearmanr(v[ok], df.log10_kcat_over_KM[ok])[0], rho_HG3=spearmanr(v[hg], df.log10_kcat_over_KM[hg])[0],
                           rho_tmol_run1=(r.rho_run1 if (r is not None and tool == "tmol") else np.nan), rho_tmol_run2=(r.rho_run2 if (r is not None and tool == "tmol") else np.nan)))
    df = df.sort_values("kcat_over_KM_M-1s-1")
    df.to_csv(out / "kemp_reference_dataset_selected.csv", index=False, float_format="%.5g")
    S = pd.DataFrame(scores); S.to_csv(out / "kemp_selected_scores_list.csv", index=False, float_format="%.3g")
    print(df.shape, "->", out / "kemp_reference_dataset_selected.csv"); print(S[["proposal_id", "score_column_in_this_csv", "pipeline_column_id", "direction", "rho_all19", "rho_HG3"]].round(2).to_string(index=False))

    # ---- plots (replace the previous figures; keep them in figures/_previous)
    fg = out / "figures"; old = fg / "_previous"; old.mkdir(exist_ok=True)
    for f in fg.glob("*.png"):
        shutil.move(str(f), str(old / f.name))
    from .plot_best import grid
    cids = [PICKS[i - 1][0] for i in ids]
    Rall = R.copy()
    stats = {}
    for c in cids:
        r = Rall.loc[c]; ex_ = f"HG3-only ρ = {r.rho_HG3_only:+.2f}"
        if np.isfinite(r.rho_run1) and PICKS[[p[0] for p in PICKS].index(c)][2] == "tmol":
            ex_ += f"\ntmol runs: {r.rho_run1:+.2f} / {r.rho_run2:+.2f}"
        stats[c] = dict(rho=r.rho_all19_tmolmean, q_BH=r.q_BH_all19, rho_loo_min=r.rho_loo_min, rho_loo_max=r.rho_loo_max, rho_within_family=r.rho_within_family, extra=ex_)
    grid(cids, M, lab, stats, fg / "selected_scores_all19.png", "Selected ranking scores vs activity: all 19 proteins",
         ncols=4, sub="Spearman ρ vs log10 kcat/KM (tmol values = mean of two relax runs). q = BH over 2,174 screened candidates; LOPO = ρ range dropping any one protein; within-family = ranks inside HG3 / 1A53 / KABLE.")
    MH = pd.read_csv(work / "report/everything_tmolmean_HG3/candidate_matrix.csv", index_col=0); CH = pd.read_csv(work / "report/everything_tmolmean_HG3/candidates_all.csv").set_index("candidate")
    labh = lab[lab.family == "HG3"]
    sh = {c: dict(rho=CH.loc[c, "rho"], q_BH=CH.loc[c, "q_BH"], rho_loo_min=CH.loc[c, "rho_loo_min"], rho_loo_max=CH.loc[c, "rho_loo_max"], rho_within_family=np.nan) for c in cids if c in CH.index}
    grid([c for c in cids if c in MH.columns], MH, labh, sh, fg / "selected_scores_HG3.png", "Selected ranking scores vs activity: HG3 lineage only (n = 12)",
         ncols=4, sub="Same scores within the HG3 scaffold. q = BH within the HG3 screen (2,146 candidates). n = 12: hypothesis-generating.")
    print("figures ->", fg)


if __name__ == "__main__":
    main()
