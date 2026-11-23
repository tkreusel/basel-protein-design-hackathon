"""Stage 4b (CPU): compare the 1x Boltz-2 TSA refolds with the RFD4 designs and gate them.

usage: boltz_gate.py --boltz-out DIR --seq sequences.csv --designs gate.csv --out boltz1_eval.csv [--survivors survivors.txt]
                     [--min-keep 10] [--thr plddt=80 ca=2 ...]

Per prediction: Boltz confidence (ipTM, ligand ipTM, complex pLDDT, ipDE; protein pLDDT from the CIF B-factors), CA RMSD to the
design, ligand RMSD after CA superposition, RMSD of the fixed motif tips, and the pocket/base/nitro metrics recomputed on the
predicted complex of the designed sequence (so the sequence-dependent terms are judged on the real sequence).

Strict gate (config/gates.yaml): pLDDT, CA RMSD, ligand ipTM, ligand RMSD, motif-tip RMSD and the refold geometry all pass.
Fallback: if fewer than --min-keep sequences of the run pass, the best of the rest are added (most criteria passed, then
ligand RMSD, then ligand ipTM); they must at least fold (pLDDT and CA RMSD pass; sheet runs are exempt).
boltz_tier 1 = strict, 2 = fallback; survivors = tier 1 and 2.
"""
import argparse
import glob
import json
import os

os.environ.setdefault("CLUSTER", "local")
import numpy as np
import pandas as pd

from . import config

CRITERIA = ["ok_plddt", "ok_ca", "ok_liptm", "ok_ligand", "ok_motif", "geom_ok_refold"]


def plddt_scale(protein_plddt):
    """Boltz writes pLDDT as 0-1 in the B-factor column; some tools use 0-100. Decide from the data."""
    return 100.0 if protein_plddt.max() > 1.5 else 1.0


def apply_gate(t, thr=None, min_keep=None, design_gate=None, motifs=None):
    """Add the pass flags, tiers and fallback to a table of per-prediction measurements."""
    cfg = config.load("gates")
    thr = {**cfg["boltz_gate"], **(thr or {})}
    min_keep = thr["min_keep_sequences"] if min_keep is None else min_keep
    g = (design_gate or cfg["design_gate"])
    geo, seq = g["geometry"], g["sequence"]
    motifs = motifs or config.load("motifs")["motifs"]
    t = t.copy()
    polar_allowed = t.motif.map(lambda m: motifs[m]["extra_polar_near_base"])

    scale = plddt_scale(t.protein_plddt)
    t["ok_plddt"] = t.protein_plddt >= thr["plddt_min"] / 100 * scale
    t["ok_ca"] = t.ca_rmsd <= thr["ca_rmsd_max"]
    t["ok_liptm"] = t.ligand_iptm >= thr["ligand_iptm_min"]
    t["ok_ligand"] = t.ligand_rmsd_ca_aligned <= thr["ligand_rmsd_max"]
    t["ok_motif"] = t.motif_tip_rmsd_ca_aligned <= thr["motif_tip_rmsd_max"]
    t["geom_ok_refold"] = (t.b_volume_ratio.between(*geo["volume_ratio"]) & (t.b_base_sasa_rel <= geo["base_sasa_rel_max"])
                           & (t.b_base_sc_polar_4p5 <= polar_allowed) & (t.b_base_cation_6 <= seq["base_cation_6_max"])
                           & (t.b_nitro_hb_far >= seq["nitro_hb_far_min"]) & (t.b_nitro_hb_near_base <= seq["nitro_hb_near_base_max"]))
    t["boltz_gate_strict"] = t[CRITERIA].all(axis=1)
    t["n_crit_ok"] = t[CRITERIA].sum(axis=1)
    t["boltz_tier"] = np.where(t.boltz_gate_strict, 1, 0)

    need = min_keep - int(t.boltz_gate_strict.sum())
    if need > 0:
        folds = (t.ok_plddt & t.ok_ca) | t.arm.str.contains("sheet")
        pool = t[(t.boltz_tier == 0) & folds].sort_values(["n_crit_ok", "ca_rmsd", "ligand_rmsd_ca_aligned"], ascending=[False, True, True])
        t.loc[pool.head(need).index, "boltz_tier"] = 2
    t["boltz_gate"] = t.boltz_tier > 0
    return t


def strictness(t):
    n_all = t[CRITERIA].sum(axis=1)
    return pd.DataFrame([dict(criterion=c, pass_rate=round(float(t[c].mean()), 3), fails_only_this=int(((n_all == len(CRITERIA) - 1) & ~t[c]).sum()))
                         for c in CRITERIA])


def measure(cif, boltz_json, design_cif, fixed, metrics):
    """Measurements of one prediction against its design. `fixed` = 'A127 A50' (motif residues of the design)."""
    import biotite.structure as bs
    from . import structure
    tips = structure.residue_tips()
    pred, design = structure.load_cif(cif, extra_fields=["b_factor"]), structure.load_cif(design_cif)
    pp, pl, dp, dl = pred[~pred.hetero], pred[pred.hetero], design[~design.hetero], design[design.hetero]
    pca, dca = pp[pp.atom_name == "CA"], dp[dp.atom_name == "CA"]
    conf = json.load(open(boltz_json))
    row = dict(length=len(pca), complex_plddt=conf.get("complex_plddt"), ptm=conf.get("ptm"), iptm=conf.get("iptm"), ligand_iptm=conf.get("ligand_iptm"),
               confidence=conf.get("confidence_score"), complex_iplddt=conf.get("complex_iplddt"), complex_ipde=conf.get("complex_ipde"),
               complex_pde=conf.get("complex_pde"), protein_plddt=float(np.mean(pp.b_factor)), ligand_plddt=float(np.mean(pl.b_factor)))
    if len(pca) == len(dca):
        _, tr = bs.superimpose(dca, pca)
        row["ca_rmsd"] = float(bs.rmsd(dca, tr.apply(pca)))
        shared = [n for n in dl.atom_name if n in set(pl.atom_name)]
        row["ligand_rmsd_ca_aligned"] = structure.rmsd_xyz(tr.apply(structure.atom_xyz(pl, shared)), structure.atom_xyz(dl, shared))
        tip_rmsd = []
        for token in str(fixed).split():
            rid = int(token[1:])
            ref_res, pred_res = dp[dp.res_id == rid], pp[pp.res_id == rid]
            names = tips.get(ref_res.res_name[0], ())
            if names and all(n in set(pred_res.atom_name) for n in names):
                tip_rmsd.append(structure.rmsd_xyz(structure.atom_xyz(ref_res, names), tr.apply(structure.atom_xyz(pred_res, names))))
        row["motif_tip_rmsd_ca_aligned"] = float(np.max(tip_rmsd)) if tip_rmsd else np.nan
    row.update({"b_" + k: v for k, v in {**metrics["pocket"].compute(pred), **metrics["base"].compute(pred), **metrics["channel"].compute(pred)}.items()})
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--boltz-out", required=True, help="boltz_out/6NT (searched recursively)")
    ap.add_argument("--seq", required=True)
    ap.add_argument("--designs", required=True, help="gate.csv of stage 2 (needs design and cif columns)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--survivors")
    ap.add_argument("--min-keep", type=int)
    ap.add_argument("--thr", nargs="*", default=[], help="override thresholds, e.g. ligand_rmsd_max=2.5")
    a = ap.parse_args()

    from rfproteina.metrics.custom_metrics.kemp_base import KempBaseMetric
    from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
    from rfproteina.metrics.custom_metrics.kemp_channel import KempChannelMetric
    metrics = dict(base=KempBaseMetric(), pocket=KempPocketMetric(target_ratio=0.5, tol=0.10), channel=KempChannelMetric())
    seqs = pd.read_csv(a.seq).set_index("seq_id")
    designs = pd.read_csv(a.designs).drop_duplicates("design").set_index("design")

    rows = []
    for f in sorted(glob.glob(a.boltz_out + "/**/predictions/*/*_model_0.cif", recursive=True)):
        pred = os.path.basename(os.path.dirname(f))
        sid = pred.rsplit("__", 1)[0]
        if sid not in seqs.index:
            continue
        design = seqs.loc[sid, "design"]
        try:
            row = measure(f, f"{os.path.dirname(f)}/confidence_{pred}_model_0.json", designs.loc[design, "cif"], seqs.loc[sid, "fixed"], metrics)
            rows.append(dict(seq_id=sid, design=design, source=seqs.loc[sid, "source"], variant=seqs.loc[sid, "variant"],
                             motif=designs.loc[design, "motif"], arm=designs.loc[design, "arm"], **row))
        except Exception as exc:
            print("fail", pred, repr(exc)[:150], flush=True)

    thr = {}
    for kv in a.thr:
        k, v = kv.split("=")
        thr[k] = float(v)
    t = apply_gate(pd.DataFrame(rows), thr, a.min_keep)
    t.to_csv(a.out, index=False)
    s = strictness(t)
    s.to_csv(a.out.replace(".csv", "_strictness.csv"), index=False)
    print(s.to_string(index=False))
    if a.survivors:
        open(a.survivors, "w").write("\n".join(t[t.boltz_gate].seq_id) + "\n")
    print(len(t), "predictions;", int(t.boltz_gate_strict.sum()), "strict pass,", int((t.boltz_tier == 2).sum()), "fallback")


if __name__ == "__main__":
    main()
