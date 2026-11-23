"""Apply the steering handoff's computable core criteria to our stage-1 runs."""
import sys, csv, glob, warnings
from pathlib import Path
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from steering_criteria import measure
from kemp.record import read_manifest
from biotite.structure.io.pdbx import CIFFile, get_structure

RUNS = {"full motif, no burial (job13)": "runs/2026-11-21_rfd4-backbones_job13_32designs_seed101",
        "full motif + BURIED (job21)": "runs/2026-11-21_rfd4-backbones-BURIED_job21_32designs_seed201",
        "tip atoms + BURIED (job23)": "runs/2026-11-21_rfd4-backbones-TIPATOMS-BURIED_job23_32designs_seed202"}

def locate(run, rec):
    p = rec.paths.get("rfd4_cif")
    if p and Path(p).exists():
        return p
    h = glob.glob(f"{run}/attempts/*_{rec.design_id}_*.cif")   # raw dir was removed
    return h[0] if h else None

def crit(m):
    return dict(lig=m["lig_sasa_frac"] <= 0.10, base_buried=m["base_sasa_rel"] <= 0.10,
                polar=m["base_sc_polar_4p5"] <= 1, nocat=m["base_cation_6"] == 0,
                nitro=(m["nitro_hb_far"] >= 1 and m["nitro_hb_near_base"] == 0),
                clash=m["lig_clash"] == 0)

rows = []
for lab, run in RUNS.items():
    for r in read_manifest(f"{run}/manifest.jsonl"):
        p = locate(run, r)
        if not p:
            print("missing structure:", lab, r.design_id); continue
        m = measure(get_structure(CIFFile.read(p), model=1))
        if not m:
            print("no measurement:", lab, r.design_id); continue
        m.update(run=lab, design_id=r.design_id, our_gate=(r.gates.get("stage1_motif") == "pass"))
        rows.append(m)

cols = sorted({k for r in rows for k in r})
with open("results_steering_criteria_on_ours.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

print(f"{'run':32}{'n':>4}{'ourgate':>8}{'lig':>5}{'baseB':>6}{'pol<=1':>7}{'0cat':>5}{'nitro':>6}{'clash':>6}{'ALL':>5}{'pol+cat':>8}")
for lab in RUNS:
    s = [r for r in rows if r["run"] == lab]; cs = [crit(r) for r in s]
    f = lambda k: sum(c[k] for c in cs)
    print(f"{lab:32}{len(s):>4}{sum(r['our_gate'] for r in s):>8}{f('lig'):>5}{f('base_buried'):>6}{f('polar'):>7}{f('nocat'):>5}{f('nitro'):>6}{f('clash'):>6}"
          f"{sum(all(c.values()) for c in cs):>5}{sum(c['polar'] and c['nocat'] for c in cs):>8}")
print("\nmeans (polar atoms<=4.5A, cations<=6A, base rel SASA, nitro H-bonds far):")
for lab in RUNS:
    s = [r for r in rows if r["run"] == lab]
    print(f"  {lab:32} polar {np.mean([r['base_sc_polar_4p5'] for r in s]):.2f}  cation {np.mean([r['base_cation_6'] for r in s]):.2f}"
          f"  baseSASA {np.mean([r['base_sasa_rel'] for r in s]):.3f}  nitroHBfar {np.mean([r['nitro_hb_far'] for r in s]):.2f}")
