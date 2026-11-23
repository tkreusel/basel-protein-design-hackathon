"""Copy the structure files of the current final top 20 (departure/final_designs_name_mapping.csv slots -> rerank_top20.csv) into
/mnt/work/projects/departure/final_top20_structures/{boltz_6NT_refold,rfd4_design}/<SLOT>.cif, plus index.csv.
Sources: the earlier structure folder (runs/2026-11-22_top30-structures_by-score_no-eligibility-filter) where the design is present, else searched directly.
Usage: python -I export_final_structures.py"""
import csv, glob, re, shutil
from pathlib import Path
T = Path("/mnt/work/projects/tobi_v1"); D = Path("/mnt/work/projects/departure")
OLD = T / "runs/2026-11-22_top30-structures_by-score_no-eligibility-filter"
M30 = Path("/mnt/work/projects/Maxi/kemp_results_20261122/final_strict_all_gates_ipsae0.8/top30")
R = T / "runs"
top = list(csv.DictReader(open(T / "submission/rerank_top20.csv")))
fa = [l[1:].split()[0] for l in open(D / ("designs.fasta" if (D / "designs.fasta").exists() else "final_designs.fasta")) if l.startswith(">")]
out = D / "final_top20_structures"
for s in ("boltz_6NT_refold", "rfd4_design"): (out / s).mkdir(parents=True, exist_ok=True)
for f in list(out.glob("*/*.cif")) + list(out.glob("index.csv")): f.unlink()
def old(sub, fid):
    h = [p for p in (OLD / sub).glob("*.cif") if p.stem.split("_", 1)[1] == fid]; return h[0] if h else None
idx = []
for slot, r in zip(fa, top):
    fid = re.sub(r"^\d\d_", "", r["orig_id"]) if r["source"] == "maxi_overnight" else r["orig_id"]
    b, d = old("boltz_6NT_refold", fid), old("rfd4_design", fid)
    if r["source"] == "maxi_overnight" and (b is None or d is None):
        nn = r["orig_id"][:2]; h = glob.glob(f"{M30}/{nn}_*")
        if h: b = b or Path(h[0]) / "boltz_TSA_6NT_best.cif"; d = d or Path(h[0]) / "design_RFD4.cif"
    if b is None:
        h = [p for pat in ("pool-v1_5seeds/seed_0", "overnight-top24_5seeds/seed_0") for p in glob.glob(f"{R}/*{pat}/**/predictions/{fid}/{fid}_model_0.cif", recursive=True)]
        b = Path(h[0]) if h else None
    if b: shutil.copy2(b, out / "boltz_6NT_refold" / f"{slot}.cif")
    if d: shutil.copy2(d, out / "rfd4_design" / f"{slot}.cif")
    idx.append(dict(slot=slot, design_id=fid, rank=r["rank"], base_res=r["base"], oxy_res=r["oxy"], base_O_N3_A=r["b"], donor_A=r["d"], boltz_refold=bool(b), rfd4_design=bool(d)))
w = csv.DictWriter(open(out / "index.csv", "w"), fieldnames=list(idx[0])); w.writeheader(); w.writerows(idx)
(out / "README.md").write_text("""# Final top-20 design structures (eligibility-filtered, current ranking)
- `boltz_6NT_refold/<slot>.cif` Boltz-2 refold with 6NT (what the hosts score): protein chain A, ligand chain B. Maxi overnight designs: his best of 10 TSA samples; others: our seed-0 refold.
- `rfd4_design/<slot>.cif` the RFD4 backbone + ligand the sequence was designed on (missing where not kept; see index.csv).
- Slot names = `final_designs.fasta` (SubstrateVentures01..20); `index.csv` maps slot -> design id, base/oxy residue (1-based), distances.
Catalytic residues: base Asp/Glu = base_res, oxyanion donor = oxy_res. Superpose on chain A CA atoms.
Regenerate after ranking changes: python -I tobi_v1/submission/export_final_structures.py
""")
print(sum(i["boltz_refold"] for i in idx), "refolds,", sum(i["rfd4_design"] for i in idx), "design backbones of", len(idx))
for i in idx:
    if not (i["boltz_refold"] and i["rfd4_design"]): print("incomplete:", i["slot"], i["design_id"])
