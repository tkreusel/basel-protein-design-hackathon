"""Eligibility gate for the submission: <70 % sequence identity to any published enzyme, TM-score <= 0.7 to any
published structure, length <= 260. Local tools only; no sequence or structure leaves the cluster.

Inputs  : candidates CSV with columns  cand_id,sequence,structure   (structure = predicted model, Boltz refold CIF/PDB)
Stages  : 1 length / alphabet
          2 MMseqs2 vs PDB seqres (all published PDB protein chains)            -> CPU, head node, minutes
          3 MMseqs2 vs UniRef90-2022-05 (=targetDB_gpu)                          -> GPU job (eligibility/uniref_search.sbatch);
                                                                                  its TSV is read here if present
          4 Foldseek vs the PDB 30%-cluster representatives (db/pdb30reps), top 20 hits by TM-score, then
            US-align on each hit, BOTH normalisations, higher one reported (TM-score is asymmetric)
Identity: reported as global identity (matches / query length) AND identity over the aligned region with
          coverage; a >= 90 % identical stretch of >= 50 residues is flagged even when global identity passes.
Output  : eligibility/results/eligibility.csv and eligibility.md
Usage   : python -I check_eligibility.py <candidates.csv> [--skip-structure] [--out results]
Run in the kemp-orch env (biotite). Thresholds are the hackathon's (GUIDE / plan); WARN = within 10 % of a limit.
"""
import argparse, csv, subprocess, sys, tempfile, os
from pathlib import Path
import numpy as np
from biotite.structure.io.pdbx import CIFFile, get_structure
from biotite.structure.io.pdb import PDBFile
import biotite.structure.io.pdb as pdbio

HERE = Path(__file__).resolve().parent
ENV = os.environ
MMSEQS = "/mnt/ref/software/mmseqs2-18-8cc5c/bin/mmseqs"
FOLDSEEK = "/mnt/work/repos/foldseek_dl/foldseek/bin/foldseek"
USALIGN = "/mnt/work/repos/usalign/USalign"
ID_LIMIT, TM_LIMIT, LEN_LIMIT = 0.70, 0.70, 260
STD = set("ACDEFGHIKLMNPQRSTVWY")
FMT = "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov"

ap = argparse.ArgumentParser()
ap.add_argument("candidates"); ap.add_argument("--out", default=str(HERE / "results"))
ap.add_argument("--skip-structure", action="store_true"); ap.add_argument("--reuse", action="store_true", help="reuse existing mmseqs/foldseek hit tables"); ap.add_argument("--threads", type=int, default=3)
ap.add_argument("--uniref-tsv", default=str(HERE / "results/uniref90_hits.tsv"))
a = ap.parse_args()
out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
rows = list(csv.DictReader(open(a.candidates)))
tmp = Path(tempfile.mkdtemp(prefix="elig-", dir=HERE / "tmp"))


def run(cmd):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode:
        raise SystemExit(f"FAILED {' '.join(map(str, cmd))[:200]}\n{p.stderr[-800:]}")
    return p.stdout


# ---- 1 length / alphabet --------------------------------------------------------
res = {r["cand_id"]: dict(cand_id=r["cand_id"], length=len(r["sequence"]), alphabet_ok=set(r["sequence"]) <= STD) for r in rows}
fa = tmp / "q.fasta"
fa.write_text("".join(f">{r['cand_id']}\n{r['sequence']}\n" for r in rows))


def summarise(tsv: Path, tag: str):
    best = {}
    stretch = {}
    if tsv.exists():
        for line in tsv.read_text().splitlines():
            f = line.split("\t"); q, t = f[0], f[1]
            pid, aln, qlen, qcov = float(f[2]) / (100 if float(f[2]) > 1 else 1), int(f[3]), int(f[12]), float(f[14])
            g = pid * aln / qlen                                    # global identity over the whole query
            b = best.get(q)
            if b is None or g > b["glob"]:
                best[q] = dict(glob=g, aln_id=pid, cov=qcov, hit=t)
            if pid >= 0.90 and aln >= 50:
                stretch[q] = max(stretch.get(q, 0), aln)
    for cid, r in res.items():
        b = best.get(cid)
        r[f"{tag}_glob_id"] = round(b["glob"], 3) if b else 0.0
        r[f"{tag}_aln_id"] = round(b["aln_id"], 3) if b else ""
        r[f"{tag}_cov"] = round(b["cov"], 2) if b else ""
        r[f"{tag}_hit"] = b["hit"] if b else "no hit with e<1e-3"
        r[f"{tag}_stretch90"] = stretch.get(cid, 0)


# ---- 2 PDB seqres --------------------------------------------------------------
db = HERE / "db/pdbseqresDB"
if db.with_suffix(".dbtype").exists():
    pdb_tsv = out / "pdb_seqres_hits.tsv"
    if not (a.reuse and pdb_tsv.exists()): run([MMSEQS, "easy-search", str(fa), str(db), str(pdb_tsv), str(tmp / "mm"), "-s", "7.5", "--max-seqs", "2000",
         "-e", "1e-3", "--min-seq-id", "0.0", "--format-output", FMT, "--threads", str(a.threads), "-v", "1"])
    summarise(pdb_tsv, "pdb")
else:
    print("WARNING: db/pdbseqresDB missing; PDB identity not checked")

# ---- 3 UniRef90 (GPU job output) ----------------------------------------------
summarise(Path(a.uniref_tsv), "uniref")
have_uniref = Path(a.uniref_tsv).exists()

# ---- 4 structure ---------------------------------------------------------------
def to_pdb(path: str, dst: Path):
    p = Path(path)
    s = get_structure(CIFFile.read(p), model=1) if p.suffix == ".cif" else pdbio.PDBFile.read(p).get_structure(model=1)
    s = s[(~s.hetero) & (s.chain_id == s.chain_id[~s.hetero][0])]
    f = PDBFile(); f.set_structure(s); f.write(dst)


if not a.skip_structure:
    fsdb = HERE / "db/pdb30reps"
    if not fsdb.with_suffix(".dbtype").exists() and not Path(str(fsdb)).exists():
        print("WARNING: db/pdb30reps missing; TM-score not checked")
    else:
        qdir = tmp / "qpdb"; qdir.mkdir()
        for r in rows:
            to_pdb(r["structure"], qdir / f"{r['cand_id']}.pdb")
        hits = out / "foldseek_hits.tsv"
        if not (a.reuse and hits.exists()): run([FOLDSEEK, "easy-search", str(qdir), str(fsdb), str(hits), str(tmp / "fs"), "--alignment-type", "1",
             "--format-output", "query,target,alntmscore,qtmscore,ttmscore,fident,alnlen,evalue", "-e", "10",
             "--max-seqs", "3000", "--threads", str(a.threads), "-v", "1"])
        top = {}
        for line in hits.read_text().splitlines():
            f = line.split("\t"); q = f[0][:-4] if f[0].endswith(".pdb") else f[0]; top.setdefault(q, []).append((float(f[2]), f[1]))
        # exact TM with US-align on the top 20 hits (both normalisations)
        from biotite.structure.io.pdbx import CIFFile as CF
        cache = {}
        for r in rows:
            cid = r["cand_id"]
            lst, seen_e = [], set()
            for sc, tn in sorted(top.get(cid, []), reverse=True):      # best hit per ENTRY (NMR ensembles would flood the top 20)
                e0 = tn.split("_")[0].lower()
                if e0 not in seen_e:
                    seen_e.add(e0); lst.append((sc, tn))
                if len(lst) >= 20: break
            best = (0.0, "")
            for _, tname in lst:
                entry, chain = tname.split("_")[0].split(".cif")[0].lower(), tname.rsplit("_", 1)[-1]
                key = (entry, chain)
                if key not in cache:
                    try:
                        import gzip, io
                        with gzip.open(HERE / f"db/pdb30_cif/{entry}.cif.gz", "rt") as fh:
                            cf = CF.read(fh)
                        st = get_structure(cf, model=1, use_author_fields=True)
                        st = st[(~st.hetero) & (st.chain_id == chain)]
                        pth = tmp / f"t_{entry}_{chain}.pdb"; f = PDBFile(); f.set_structure(st); f.write(pth); cache[key] = pth
                    except Exception:
                        cache[key] = None
                if cache[key] is None:
                    continue
                o = run([USALIGN, str(qdir / f"{cid}.pdb"), str(cache[key]), "-TMscore", "0", "-outfmt", "2"])
                for ln in o.splitlines():
                    if ln.startswith("#") or not ln.strip():
                        continue
                    f2 = ln.split("\t"); tm = max(float(f2[2]), float(f2[3]))
                    if tm > best[0]:
                        best = (tm, f"{entry.upper()}_{chain}")
            res[cid]["tm_max"] = round(best[0], 3); res[cid]["tm_hit"] = best[1]; res[cid]["tm_n_hits_checked"] = len(lst)

# ---- verdicts ----------------------------------------------------------------
for r in res.values():
    why, warn = [], []
    if r["length"] > LEN_LIMIT: why.append("length")
    if not r["alphabet_ok"]: why.append("non-standard residue")
    for tag in ("pdb", "uniref"):
        g = r.get(f"{tag}_glob_id")
        if g not in (None, ""):
            if g >= ID_LIMIT: why.append(f"{tag} identity {g:.2f}")
            elif g >= 0.9 * ID_LIMIT: warn.append(f"{tag} identity {g:.2f}")
        if r.get(f"{tag}_stretch90", 0): warn.append(f"{tag} >=90% stretch of {r[f'{tag}_stretch90']} aa")
        ai = r.get(f"{tag}_aln_id")
        if ai not in (None, "") and ai >= ID_LIMIT and float(r.get(f"{tag}_cov") or 0) >= 0.5: why.append(f"{tag} aligned identity {ai:.2f} over cov {r[f'{tag}_cov']}")
    if "tm_max" in r:
        if r["tm_max"] > TM_LIMIT: why.append(f"TM {r['tm_max']:.2f} to {r['tm_hit']}")
        elif r["tm_max"] > 0.9 * TM_LIMIT: warn.append(f"TM {r['tm_max']:.2f} to {r['tm_hit']}")
    else:
        warn.append("structure not checked")
    if not have_uniref: warn.append("UniRef90 not checked")
    r["verdict"] = "FAIL" if why else ("PASS-with-warnings" if warn else "PASS")
    r["reasons"] = "; ".join(why + warn)

cols = ["cand_id", "verdict", "length", "pdb_glob_id", "pdb_aln_id", "pdb_cov", "pdb_hit", "pdb_stretch90",
        "uniref_glob_id", "uniref_aln_id", "uniref_cov", "uniref_hit", "uniref_stretch90", "tm_max", "tm_hit", "reasons"]
with open(out / "eligibility.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore"); w.writeheader(); w.writerows(res.values())
print(f"{'cand':44}{'verdict':20}{'PDBid':>6}{'URid':>6}{'TM':>6}  reasons")
for r in res.values():
    print(f"{r['cand_id'][:43]:44}{r['verdict']:20}{r.get('pdb_glob_id', ''):>6}{r.get('uniref_glob_id', ''):>6}{r.get('tm_max', ''):>6}  {r['reasons']}")
