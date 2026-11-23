"""Re-rank the merged top-80 pool (Maxi overnight top 30 + all our scored designs, top 80 by the judging stand-in), eligibility-filtered.
 1. pool            merged80_pool.py -> merged80_pool.csv (80 unique sequences, judging stand-in = geometry^2 * iPDE * burial * consistency)
 2. eligibility     seq identity <0.70 (PDB seqres + UniRef90) and TM <= 0.70 (Foldseek+US-align vs PDB 30%-cluster reps); results_v5 + results_v6, joined by SEQUENCE
 3. literature L    COMMON metric subset available for every design: a ipsae|sub, b liptm sub-prod, c liptm tsa-sub, d pae_cat_lig|sub, e pae tsa-sub, f ipsae tsa-sub,
                    i ipsae sub-prod, j liptm|sub  (z vs 19 literature enzymes, sign by regression, weight = LOO Q2; g,h,k,l need tmol/PyRosetta/PLACER, not run for all)
                    L = Phi(weighted mean z). Designs with no arm data get L = 0.5 (flagged).
 4. final = judging * L ; identical sequences keep best rank only.
Outputs: rerank_top20.csv (concise judged table), rerank_all.csv, designs_rerank_top20.fasta."""
import csv, json, math, re, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
src_txt = (HERE / "build_final_table_lit.py").read_text().split("# ---- ranking source")[0]
exec(compile(src_txt, "bftl", "exec"))           # method(), submission_name(), SEQ_ORIGIN
f = lambda v: float(v) if v not in (None, "", "nan", "NaN") else np.nan

SPEC = {"a": (+1, .467), "b": (+1, .411), "c": (-1, .404), "d": (-1, .398), "e": (+1, .342), "f": (-1, .325), "i": (+1, .271), "j": (+1, .259)}
LIT = [r for r in csv.DictReader(open("/mnt/work/projects/Maxi/pipeline/work/report/deliverables/kemp_eval_dataset_full.csv")) if r.get("labelled") in ("True", "1", "true")]
def d8(b):
    return {"a": b["ip_sub"], "b": b["it_sub"] - b["it_prod"], "c": b["it_tsa"] - b["it_sub"], "d": b["pae_sub"], "e": b["pae_tsa"] - b["pae_sub"],
            "f": b["ip_tsa"] - b["ip_sub"], "i": b["ip_sub"] - b["ip_prod"], "j": b["it_sub"]}
lv = {k: [] for k in SPEC}
for r in LIT:
    b = dict(ip_tsa=f(r.get("ip.tsa.ipsae.p2l")), ip_sub=f(r.get("ip.sub.ipsae.p2l")), ip_prod=f(r.get("ip.prod.ipsae.p2l")), it_tsa=f(r.get("ip.tsa.boltz.ligand_iptm")),
             it_sub=f(r.get("ip.sub.boltz.ligand_iptm")), it_prod=f(r.get("ip.prod.boltz.ligand_iptm")), pae_tsa=f(r.get("site.pae_cat_to_lig.tsa")), pae_sub=f(r.get("site.pae_cat_to_lig.sub")))
    for k, v in d8(b).items(): lv[k].append(v)
REF = {k: (np.nanmean(v), np.nanstd(v, ddof=1)) for k, v in lv.items()}
def L8(b):
    zs = {}
    for k, (sg, q) in SPEC.items():
        v = d8(b)[k]
        if np.isfinite(v): zs[k] = float(np.clip((v - REF[k][0]) / REF[k][1], -3, 3)) * sg
    if not zs: return None, 0.5, {}
    s = sum(zs[k] * SPEC[k][1] for k in zs) / sum(SPEC[k][1] for k in zs)
    return s, 0.5 * (1 + math.erf(s / math.sqrt(2))), zs

# ---- base-value sources -------------------------------------------------------------------------------
MX = {r["seq_id"]: r for r in csv.DictReader(open("/mnt/work/projects/Maxi/kemp_results_20261122/tables/evaluated_sequences_full_suite.csv"))}
tab = {int(r["rank30"]): r["seq_id"] for r in csv.DictReader(open("/mnt/work/projects/Maxi/kemp_results_20261122/final_strict_all_gates_ipsae0.8/top30/top30_table.csv"))}
M30 = {r["sequence"]: r for r in csv.DictReader(open(HERE / "maxi_top30_scored.csv"))}
V2 = {r["sequence"]: r for r in csv.DictReader(open(HERE / "ranking_v2.csv"))}
IP = {}; ipf = ROOT / "maxi_eval/work/ipsae/ipsae_models.csv"
if ipf.exists():
    for r in csv.DictReader(open(ipf)):
        if r.get("status") == "ok": IP.setdefault((r["entry"].replace("_boltz", ""), r["arm"]), []).append(r)
def ipv(cid, arm, key):
    v = [f(r.get(key)) for r in IP.get((cid, arm), [])]; return float(np.nanmean(v)) if v else np.nan
def base_ours(cid):
    return dict(ip_tsa=ipv(cid, "tsa", "ipsae.p2l"), ip_sub=ipv(cid, "sub", "ipsae.p2l"), ip_prod=ipv(cid, "prod", "ipsae.p2l"), it_tsa=ipv(cid, "tsa", "boltz.ligand_iptm"),
                it_sub=ipv(cid, "sub", "boltz.ligand_iptm"), it_prod=ipv(cid, "prod", "boltz.ligand_iptm"),
                pae_tsa=ipv(cid, "tsa", "site.pae_cat_to_lig"), pae_sub=ipv(cid, "sub", "site.pae_cat_to_lig"))
def base_maxi(sid):
    m = MX[sid]
    return dict(ip_tsa=f(m.get("ipsae.tsa")), ip_sub=f(m.get("ipsae.sub")), ip_prod=f(m.get("ipsae.prod")), it_tsa=f(m.get("6NT.ligand_iptm.mean")), it_sub=f(m.get("H5J.ligand_iptm.mean")),
                it_prod=f(m.get("6VP.ligand_iptm.mean")), pae_tsa=f(m.get("site.pae_cat_to_lig.tsa")), pae_sub=f(m.get("site.pae_cat_to_lig.sub")))

# ---- eligibility ---------------------------------------------------------------------------------------
seq_of, elig, ur = {}, {}, {}
for d, cf in (("results_v5", "candidates_v5.csv"), ("results_v6", "candidates_v6.csv")):
    for r in csv.DictReader(open(ROOT / "eligibility" / cf)): seq_of[(d, r["cand_id"])] = r["sequence"]
    p = ROOT / "eligibility" / d / "eligibility.csv"
    if p.exists():
        for r in csv.DictReader(open(p)):
            if (d, r["cand_id"]) in seq_of: elig[seq_of[(d, r["cand_id"])]] = r
    t = ROOT / "eligibility" / d / "uniref90_hits.tsv"
    if t.exists():
        for line in t.read_text().splitlines():
            x = line.split("\t"); s = seq_of.get((d, x[0]))
            if s:
                pid, aln, qlen = float(x[2]), int(x[3]), int(x[12]); pid = pid / 100 if pid > 1 else pid
                ur[s] = max(ur.get(s, 0.0), pid * aln / qlen)

def safe_method(r):
    try: return method(r)
    except Exception: return f"Maxi earlier round (handoff 2026-11-21): RFD4 beam search with base-only steering rewards (no TMol term), id {r['orig']}; Boltz-2 refolds by us."
UR_READY = (ROOT / 'eligibility/results_v6/uniref90_hits.tsv').exists()
rows = []
for p in csv.DictReader(open(HERE / "merged80_pool.csv")):
    seq = p["sequence"]; cid = p["orig"]
    if p["src"] == "maxi_overnight":
        m = M30[seq]; sid = tab[int(m["rank30"])]; b = base_maxi(sid); cid = sid
        disp = dict(base=m["base"], oxy=m["donor"], geometry=m["geometry"], b=m["b_dist"], d=m["d_dist"], s_ipde=m["s_ipde"], s_burial=m["s_burial"], cons=m["cons_pose_only"], ipsae=m["ipsae"], ipde=m["ipde"])
    else:
        v = V2[seq]; b = base_ours(cid)
        disp = dict(base=v["base_pos"], oxy=v["donor_pos"], geometry=v["geometry"], b=v["b_dist"], d=v["d_dist"], s_ipde=v["s_ipde"], s_burial=v["s_burial"], cons=v["s_consistency"], ipsae="", ipde="")
    S, L, zs = L8(b)
    e = elig.get(seq, {})
    top_r = dict(cand_id=cid, src=p["src"] if not p["src"].startswith("ours") else "ours_shortlist", orig=p["orig"] if p["src"] != "maxi_overnight" else p["orig"], total=0)
    rows.append(dict(name=submission_name(dict(top_r, src=p["src"], orig=p["orig"], cand_id=p["orig"])), source=p["src"], length=len(seq), **disp, judging=float(p["judge"]), L8=round(L, 3),
                     n_metrics=len(zs), final=float(p["judge"]) * L, max_TM=e.get("tm_max", "?"), tm_hit=e.get("tm_hit", ""), pdb_id=e.get("pdb_glob_id", "?"),
                     uniref_id=(round(ur[seq], 3) if seq in ur else (0.0 if UR_READY else "pending")), eligibility=e.get("verdict", "?"), method=safe_method(dict(top_r, src=p["src"], orig=p["orig"], cand_id=p["orig"])), sequence=seq, orig_id=p["orig"]))
def ok(r):
    if r["eligibility"] == "FAIL": return False
    if isinstance(r["uniref_id"], float) and r["uniref_id"] >= 0.70: return False
    try: return float(r["max_TM"]) <= 0.70
    except ValueError: return True
rows.sort(key=lambda r: -r["final"])
final, seen = [], set()
for r in rows:
    if ok(r) and r["sequence"] not in seen: seen.add(r["sequence"]); final.append(r)
for i, r in enumerate(final, 1): r["rank"] = i
cols = ["rank", "name", "length", "base", "oxy", "b", "d", "geometry", "s_ipde", "s_burial", "cons", "judging", "L8", "n_metrics", "final", "max_TM", "tm_hit", "pdb_id", "uniref_id", "eligibility", "source", "method", "sequence", "orig_id"]
for fn, data in (("rerank_all.csv", rows), ("rerank_top20.csv", final[:20])):
    with open(HERE / fn, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(data)
with open(HERE / "designs_rerank_top20.fasta", "w") as fh:
    for r in final[:20]: fh.write(f">{r['name']} base={r['base']} oxy={r['oxy']}\n{r['sequence']}\n")
print(f"pool {len(rows)}; eligibility-failed {sum(1 for r in rows if not ok(r))}; unknown eligibility {sum(1 for r in rows if r['eligibility']=='?')}; no-L {sum(1 for r in rows if r['n_metrics']==0)}")
for r in final[:30]:
    print(f"{r['rank']:>2} {r['name'][:44]:44} len={r['length']:>3} b={r['b']} d={r['d']} geo={r['geometry']} judge={r['judging']:.4f} L8={r['L8']} n={r['n_metrics']} final={r['final']:.4f} TM={r['max_TM']} {r['eligibility']}")
