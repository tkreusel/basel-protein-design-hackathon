"""Literature-calibrated ranking: combine the judging-criteria stand-in with the metrics that explain kcat/KM in Maxi's
univariate regressions (screenshot, ranked by leave-one-out Q2).

Per design and per metric: z = (x - mean_lit) / sd_lit over the 19 literature enzymes (kemp_eval_dataset_full.csv), clipped to +-3,
oriented by the sign of the Pearson r in the screenshot (positive z = looks more like an active enzyme), weighted by LOO Q2.
S_lit = sum(w z) / sum(w) over the metrics available for that design; L = Phi(S_lit) in (0,1).
Final score = judging stand-in (scoring/rank_v2.py: geometry^2 * iPDE * burial * consistency, geometry dominates) * L
(multiplicative, as SUBMISSION.md says the terms multiply).

Metrics used (screenshot rows; l = bz.sub-prod.cat_to_pocket_pae (-, .549) added 10:20 from Maxi's handoff_metrics S12 and our features_extra that can be computed for BOTH Maxi's designs and ours; exact duplicates counted once):
  a ipsae.p2l|sub (+, .467)   b ligand iptm sub-prod (+, .411)   c iptm tsa-sub (-, .404)   d site.pae_cat_to_lig|sub (-, .398)
  e site.pae_cat_to_lig tsa-sub (+, .342)   f ipsae.p2l tsa-sub (-, .325)   g rlxcst.field.CH.env_nobase.rdd4 (-, .325)
  h tmol E_protein_in_complex sub-prod (-, .281)   i ipsae.p2l sub-prod (+, .271)   j iptm|sub (+, .259)   k tmol int.lk_ball tsa_pk-prod (+, .258)
NOT used (need the features_extra table, not available for Maxi's designs): bz.sub-prod.cat_to_pocket_pae (Q2 .549),
pl2.tsa_holo.donor_med (Q2 .513, r +0.775: LARGER PLACER donor distance goes with activity, which contradicts the judging criteria; likely an outlier effect),
and the three ipsae.iPAE12 rows (iPAE12 is not in Maxi's design tables).
Usage: python -I lit_calibrated_rank.py"""
import csv, json, math, re
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
LIT = list(csv.DictReader(open("/mnt/work/projects/Maxi/pipeline/work/report/deliverables/kemp_eval_dataset_full.csv")))
LIT = [r for r in LIT if r.get("labelled") in ("True", "1", "true")] or LIT
MX = {r["seq_id"]: r for r in csv.DictReader(open("/mnt/work/projects/Maxi/kemp_results_20261122/tables/evaluated_sequences_full_suite.csv"))}
OURS = {r[""] if "" in r else list(r.values())[0]: r for r in csv.DictReader(open(ROOT / "maxi_eval/work/report/summary_entries.csv"))}

f = lambda r, k: (float(r[k]) if r.get(k) not in (None, "", "nan", "NaN") else np.nan)


def derive(b):
    """b: dict of base values -> the 11 metrics"""
    return {
        "a ipsae|sub": b["ip_sub"], "b liptm sub-prod": b["it_sub"] - b["it_prod"], "c liptm tsa-sub": b["it_tsa"] - b["it_sub"],
        "d pae_cat_lig|sub": b["pae_sub"], "e pae_cat_lig tsa-sub": b["pae_tsa"] - b["pae_sub"], "f ipsae tsa-sub": b["ip_tsa"] - b["ip_sub"],
        "g field CH env_nobase": b["field"], "h E_prot sub-prod": b["ep_sub"] - b["ep_prod"], "i ipsae sub-prod": b["ip_sub"] - b["ip_prod"],
        "j liptm|sub": b["it_sub"], "k lk_ball tsa-prod": b["lkb_tsa"] - b["lkb_prod"],
        "l bz cat_pocket_PAE sub-prod": b["s12"]}


SPEC = {"a ipsae|sub": (+1, .467), "b liptm sub-prod": (+1, .411), "c liptm tsa-sub": (-1, .404), "d pae_cat_lig|sub": (-1, .398),
        "e pae_cat_lig tsa-sub": (+1, .342), "f ipsae tsa-sub": (-1, .325), "g field CH env_nobase": (-1, .325), "h E_prot sub-prod": (-1, .281),
        "i ipsae sub-prod": (+1, .271), "j liptm|sub": (+1, .259), "k lk_ball tsa-prod": (+1, .258),
        "l bz cat_pocket_PAE sub-prod": (-1, .549)}

# literature reference distribution
litvals = {k: [] for k in SPEC}
for r in LIT:
    b = dict(ip_tsa=f(r, "ip.tsa.ipsae.p2l"), ip_sub=f(r, "ip.sub.ipsae.p2l"), ip_prod=f(r, "ip.prod.ipsae.p2l"),
             it_tsa=f(r, "ip.tsa.boltz.ligand_iptm"), it_sub=f(r, "ip.sub.boltz.ligand_iptm"), it_prod=f(r, "ip.prod.boltz.ligand_iptm"),
             pae_tsa=f(r, "site.pae_cat_to_lig.tsa"), pae_sub=f(r, "site.pae_cat_to_lig.sub"), field=f(r, "rlxcst.field.CH.env_nobase.rdd4"),
             ep_sub=f(r, "sub.E_protein_in_complex"), ep_prod=f(r, "prod.E_protein_in_complex"),
             lkb_tsa=f(r, "tsa_pk.int.lk_ball"), lkb_prod=f(r, "prod.int.lk_ball"), s12=f(r, "x:bz.sub-prod.cat_to_pocket_pae"))
    for k, v in derive(b).items():
        litvals[k].append(v)
REF = {k: (np.nanmean(v), np.nanstd(v, ddof=1)) for k, v in litvals.items()}


def lit_score(d):
    zs, ws, parts = [], [], {}
    for k, (sgn, q2) in SPEC.items():
        v = d[k]
        if not np.isfinite(v): continue
        z = float(np.clip((v - REF[k][0]) / REF[k][1], -3, 3)) * sgn
        zs.append(z * q2); ws.append(q2); parts[k] = round(z, 2)
    s = sum(zs) / sum(ws) if ws else 0.0
    return s, 0.5 * (1 + math.erf(s / math.sqrt(2))), parts, len(ws)


rows = []
S12 = {r["seq_id"]: f(r, "S12_boltz_cat_to_pocket_PAE_SUB_minus_PROD") for r in csv.DictReader(open("/mnt/work/projects/Maxi/kemp_results_20261122/final_strict_all_gates_ipsae0.8/handoff_metrics/top30_metrics.csv"))}
FX = {r[""]: r for r in csv.DictReader(open(ROOT / "maxi_eval/work/features_extra.csv"))}
# ---- Maxi's overnight top 30 ------------------------------------------------------------------------
M30 = {r["name"]: r for r in csv.DictReader(open(HERE / "maxi_top30_scored.csv"))}
T = Path("/mnt/work/projects/Maxi/kemp_results_20261122/final_strict_all_gates_ipsae0.8/top30")
tab = {r["seq_id"]: r for r in csv.DictReader(open(T / "top30_table.csv"))}
for name, s in M30.items():
    sid = tab and next((k for k, r in tab.items() if int(r["rank30"]) == int(s["rank30"])), None)
    m = MX[sid]
    b = dict(ip_tsa=f(m, "ipsae.tsa"), ip_sub=f(m, "ipsae.sub"), ip_prod=f(m, "ipsae.prod"),
             it_tsa=f(m, "6NT.ligand_iptm.mean"), it_sub=f(m, "H5J.ligand_iptm.mean"), it_prod=f(m, "6VP.ligand_iptm.mean"),
             pae_tsa=f(m, "site.pae_cat_to_lig.tsa"), pae_sub=f(m, "site.pae_cat_to_lig.sub"), field=f(m, "rlxcst.field.CH.env_nobase.rdd4"),
             ep_sub=f(m, "sub.E_protein_in_complex"), ep_prod=f(m, "prod.E_protein_in_complex"),
             lkb_tsa=f(m, "tsa_pk.int.lk_ball"), lkb_prod=f(m, "prod.int.lk_ball"), s12=S12.get(sid, np.nan))
    sc, L, parts, n = lit_score(derive(b))
    rows.append(dict(src="maxi_overnight", id=sid, sequence=s["sequence"], base=s["base"], donor=s["donor"], judge=float(s["total_mult"]),
                     geometry=s["geometry"], b_dist=s["b_dist"], d_dist=s["d_dist"], S_lit=round(sc, 3), L=round(L, 3), n_metrics=n, parts=json.dumps(parts)))
# ---- our evaluated finalists -------------------------------------------------------------------------
V2 = {r["cand_id"]: r for r in csv.DictReader(open(HERE / "ranking_v2.csv"))}
for ent, o in OURS.items():
    cid = ent.replace("_boltz", "")
    def iptm(arm):
        p = ROOT / "maxi_eval/boltz_root"
        d = p / (cid if arm == "tsa" else f"{cid}__{arm}")
        for jf in d.glob("confidence_*.json"):
            return float(json.load(open(jf)).get("ligand_iptm", np.nan))
        return np.nan
    b = dict(ip_tsa=f(o, "ipsae.tsa"), ip_sub=f(o, "ipsae.sub"), ip_prod=f(o, "ipsae.prod"), it_tsa=iptm("tsa"), it_sub=iptm("H5J"), it_prod=iptm("6VP"),
             pae_tsa=f(o, "site.pae_cat_to_lig.tsa"), pae_sub=f(o, "site.pae_cat_to_lig.sub"), field=f(o, "rlxcst.field.CH.env_nobase.rdd4"),
             ep_sub=f(o, "sub.E_protein_in_complex"), ep_prod=f(o, "prod.E_protein_in_complex"),
             lkb_tsa=f(o, "tsa_pk.int.lk_ball"), lkb_prod=f(o, "prod.int.lk_ball"), s12=f(FX.get(cid, {}), "bz.sub-prod.cat_to_pocket_pae"))
    v = V2.get(cid)
    if not v: continue
    sc, L, parts, n = lit_score(derive(b))
    rows.append(dict(src="ours/earlier-maxi", id=cid, sequence=v["sequence"], base=v["base_pos"], donor=v["donor_pos"], judge=float(v["total_mult"]),
                     geometry=v["geometry"], b_dist=v["b_dist"], d_dist=v["d_dist"], S_lit=round(sc, 3), L=round(L, 3), n_metrics=n, parts=json.dumps(parts)))
for r in rows:
    r["final"] = round(r["judge"] * r["L"], 5)
rows.sort(key=lambda r: -r["final"])
with open(HERE / "lit_calibrated_ranking.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("literature reference (mean, sd):", {k: (round(v[0], 3), round(v[1], 3)) for k, v in REF.items()})
for i, r in enumerate(rows, 1):
    print(f"{i:>2} {r['src'][:9]:9} {r['id'][-50:]:50} final={r['final']:.4f} judge={r['judge']:.4f} L={r['L']:.2f} (S_lit {r['S_lit']:+.2f}, n={r['n_metrics']}) b={r['b_dist']} d={r['d_dist']}")
