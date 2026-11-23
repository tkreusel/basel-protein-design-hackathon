"""Assemble the final top-30 table (and the 20-design submission FASTA) from:
  union_top45.json        ranking (stand-in score: geometry^2 * iPDE * burial * consistency; see scoring/rank_v2.py)
  eligibility results     eligibility/results_v3/eligibility.csv (+ UniRef90 hits if present)
  generation method       decoded from the design id (documented per source below)
Usage: python -I build_final_table.py"""
import csv, json, re
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
top = json.load(open(HERE / "union_top45.json"))
elig = {r["cand_id"]: r for r in csv.DictReader(open(ROOT / "eligibility/results_v3/eligibility.csv"))}
ur_tsv = ROOT / "eligibility/results_v3/uniref90_hits.tsv"
ur = {}
if ur_tsv.exists():
    for line in ur_tsv.read_text().splitlines():
        f = line.split("\t"); q = f[0]; pid, aln, qlen = float(f[2]), int(f[3]), int(f[12])
        pid = pid / 100 if pid > 1 else pid
        ur[q] = max(ur.get(q, 0.0), pid * aln / qlen)
pd_parents = {}
for p in sorted((ROOT / "runs").glob("2026-11-22_partial-diffusion_job62*/parents/parents.csv")):
    pd_parents = {r["parent"]: r for r in csv.DictReader(open(p))}

SEQ_ORIGIN = lambda tail: ("RFD4 co-designed (native) sequence" if tail == "rfd4" else
                           "LigandMPNN shell-fixed (motif + residues <=5 A of ligand / <=6 A of base kept)" if tail.startswith("shell") else
                           f"LigandMPNN full redesign #{tail[4:]} (T 0.1; K R H D E N Q S T Y C forbidden within 6 A of base carboxylate; motif fixed)")
REFOLD = "Boltz-2, 6NT, no MSA"


def method(r):
    cid, src = (re.sub(r"^\d\d_", "", r["orig"]) if r["src"] == "maxi_overnight" else r["orig"]), r["src"]
    tail = re.search(r"__?(mpnn\d+|rfd4|shell\d+)$", cid)
    tail = tail.group(1) if tail else "rfd4"
    if src == "maxi_overnight":
        m = re.match(r"hg(nt|wt)__(steer_tmol_ch18|steer_tmol|unsteered)__(b\d+_c\d+)_(no|with)_trp_(\d+)(?:_beam3_(\d+))?", cid)
        motif = "HG3.17 tips Gln50+Asp127+6NT" + (" (no Trp)" if m.group(4) == "no" else " + Trp44")
        arm = {"steer_tmol": "RFD4 reward-guided beam search: base polar/Lys-Arg/SASA/apolar, pocket-ratio, nitro-H-bond, motif-rigidity and clash rewards + TMol dG",
               "steer_tmol_ch18": "same as steer_tmol + solvent-channel bottleneck-radius band 1.7-1.9 A",
               "unsteered": "RFD4 unsteered (conditioning only), geometry-only design gate"}[m.group(2)]
        return (f"Maxi overnight run (big cluster), run {m.group(3)}, lineage {m.group(5)}" + (f", beam member {m.group(6)}" if m.group(6) else "") +
                f". {arm}. Motif {motif}. Design gate -> Boltz-2 x1 gate (pLDDT>=80, CA<=2 A, lig ipTM>=0.6, lig RMSD<=2 A, motif<=1 A) -> x9 more TSA, x10 substrate, x10 product -> ipSAE2>0.8. Sequence: " + SEQ_ORIGIN(tail))
    if src == "maxi_steered":
        m = re.match(r"hg(nt|wt)__steer_base_tmol__(no|with)_trp_(\d+)_beam3_(\d+)", cid)
        return (f"Maxi round 4/5 (2026-11-21, handoff): RFD4 beam search (4 branches x beam 3, scored at steps 30/42/52/60), rewards on base polarity/cation/SASA, pocket ratio, nitro H-bond, motif rigidity, clashes + TMol dG (weight 0.5/20), "
                f"motif HG3.17 tips Gln50+Asp127+6NT ({'no Trp' if m.group(2) == 'no' else '+ Trp44'}), lineage {m.group(3)}, beam member {m.group(4)}, seed 7. Sequence: " + SEQ_ORIGIN(tail) + f". Boltz-2 refolds (no MSA) by us, 5 seeds.")
    if src == "overnight":
        n = re.match(r"kemp4_with_trp_(\d+)", cid).group(1)
        return (f"Our overnight control (job 60): RFD4-Proteina UNSTEERED (conditioning only) with Maxi's with-Trp spec (HG3.17 tips Gln50+Asp127+Trp44 + 6NT, fixed absolute coordinates, unindexed, BURIED, hydrophobic-contact and nitro-O acceptor flags), "
                f"sample {n} of 250 (seed 7); top 100 designs by Maxi-core-criteria count -> LigandMPNN. Sequence: " + SEQ_ORIGIN(tail) + f". {REFOLD}, seed 0, then 5 seeds.")
    if src == "ours_shortlist":
        arm = "full-motif conditioning + BURIED" if cid.startswith("full-") else "tip-atom conditioning + BURIED"
        return (f"Our earlier run (2026-11-21): RFD4-Proteina unsteered, HG3.17 theozyme ({arm}); LigandMPNN (v_32_010_25, T 0.1, 4 sequences/backbone, motif fixed). "
                f"{REFOLD}, refold-funnel-validated shortlist, 5 seeds.")
    if src.startswith("child:"):
        p = src.split(":")[1]; par = pd_parents.get(p, {}).get("best_cand", "?")
        lvl = re.search(r"_f(0\.\d)_", cid).group(1)
        return (f"Our partial diffusion (job 62): RFD4 partial diffusion of parent {par} at partial_trajectory_fraction_steps={lvl} (unsteered, parent conditioning kept); "
                f"CPU backbone filter (>= parent's criteria, burial >= parent-2) -> LigandMPNN. Sequence: " + SEQ_ORIGIN(tail) + f". {REFOLD}, 1 seed then 3 seeds.")
    return "?"


def submission_name(r):
    c = re.sub(r"^\d\d_", "", r["orig"]) if r["src"] == "maxi_overnight" else r["cand_id"]
    if r["src"] == "maxi_overnight": c = "mON_" + c
    c = c.replace("maxiON_", "mON_").replace("hgnt__steer_tmol_ch18__", "nt_ch18_").replace("hgnt__steer_tmol__", "nt_tmol_").replace("hgwt__steer_tmol_ch18__", "wt_ch18_").replace("hgwt__steer_tmol__", "wt_tmol_")
    c = c.replace("hgnt__steer_base_tmol__no_trp_0_", "mR5_nt_").replace("hgwt__steer_base_tmol__with_trp_0_", "mR5_wt_")
    c = re.sub(r"[^A-Za-z0-9_-]", "_", c).replace("__", "_")
    return c[:64]


rows = []
for i, r in enumerate(top, 1):
    cid = r["cand_id"].replace(".", "p")
    e = elig.get(cid, {})
    rows.append(dict(
        rank=i, name=submission_name(r), source=r["src"], length=len(r["sequence"]), base=r["base"], oxy=r["donor"],
        score=round(r["total"], 4), geometry=r["geometry"], base_O_N3_A=r["b"], donor_A=r["d"], iPDE_A=r["ipde"] or "", ipSAE2=r["ipsae"] or "",
        consistency=r["cons"], max_TM=e.get("tm_max", "?"), tm_hit=e.get("tm_hit", ""), pdb_identity=e.get("pdb_glob_id", "?"),
        uniref_identity=round(ur.get(cid, 0.0), 3) if ur_tsv.exists() else "pending", eligibility=e.get("verdict", "?"),
        method=method(r), sequence=r["sequence"], orig_id=r["orig"]))
ok = [r for r in rows if r["eligibility"] != "FAIL" and not (isinstance(r["max_TM"], str) is False and float(r["max_TM"]) > 0.70)]
final = ok[:30]
for i, r in enumerate(final, 1):
    r["final_rank"] = i; r["in_submission_top20"] = i <= 20
with open(HERE / "final_top30.csv", "w") as fh:
    cols = ["final_rank", "in_submission_top20", "name", "source", "length", "base", "oxy", "score", "geometry", "base_O_N3_A", "donor_A", "iPDE_A", "ipSAE2",
            "consistency", "max_TM", "tm_hit", "pdb_identity", "uniref_identity", "eligibility", "method", "sequence", "orig_id"]
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(final)
with open(HERE / "designs_final20.fasta", "w") as fh:
    for r in final[:20]:
        fh.write(f">{r['name']} base={r['base']} oxy={r['oxy']}\n{r['sequence']}\n")
print(f"{len(final)} designs; {sum(1 for r in rows if r not in ok)} removed for eligibility")
for r in final:
    print(f"{r['final_rank']:>2}{'*' if r['in_submission_top20'] else ' '} {r['name'][:46]:46} {r['source'][:9]:9} len={r['length']:>3} score={r['score']:<6} b={r['base_O_N3_A']} d={r['donor_A']} TM={r['max_TM']} {r['eligibility']}")
