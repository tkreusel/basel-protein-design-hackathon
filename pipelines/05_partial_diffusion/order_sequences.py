"""After LigandMPNN: keep only MPNN sequences (the RFD4 co-designed sequence of a
partial-diffusion child is not refolded: it costs a Boltz call and our earlier
refold funnel favoured LigandMPNN sequences), and order them so every parent's
best-ranked children are refolded first. Also writes the candidates file the
host-style scorer reads. Usage: python -I order_sequences.py <run dir>"""
import csv, shutil, sys
from pathlib import Path

run = Path(sys.argv[1])
sel = {r["design"]: r for r in csv.DictReader(open(run / "children/selection.csv"))}
seqs = list(csv.DictReader(open(run / "sequences.csv")))
if not (run / "sequences_all.csv").exists():
    shutil.copy2(run / "sequences.csv", run / "sequences_all.csv")
keep = [r for r in seqs if r["source"].startswith("mpnn")]
keep.sort(key=lambda r: (int(sel[r["design"]]["select_rank"]), r["seq_id"]))
with open(run / "sequences.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(seqs[0])); w.writeheader(); w.writerows(keep)
with open(run / "child_candidates.csv", "w") as fh:
    w = csv.writer(fh); w.writerow(["cand_id", "sequence", "base_pos", "donor_pos", "wall_pos", "source", "priority"])
    for i, r in enumerate(keep):
        b, d, wl = (int(x[1:]) for x in r["fixed"].split())
        w.writerow([r["seq_id"], r["sequence"], b, d, wl, f"child:{sel[r['design']]['parent']}", i])
print(f"{len(keep)} MPNN sequences for {len(sel)} children (dropped {len(seqs) - len(keep)} RFD4-sequence rows)")
