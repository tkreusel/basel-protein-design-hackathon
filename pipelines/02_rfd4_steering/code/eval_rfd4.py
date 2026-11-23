"""Post-hoc evaluation of RFD4 trial outputs with KempPocketMetric (same metric/params as eval_rfd3).
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python eval_rfd4.py  -> out/20261121_eval_rfd4.csv"""
import glob, os, numpy as np, pandas as pd
import biotite.structure.io.pdbx as pdbx
from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
H = os.path.expanduser("~/HACKATHON_ETH/rfd4_runs"); m = KempPocketMetric(); rows = []
for d in sorted(glob.glob(f"{H}/out/20261121_*")):
    arm = os.path.basename(d)[9:]
    if arm in ("smoke",) or not os.path.isdir(d + "/inference_samples"): continue
    for f in sorted(glob.glob(d + "/inference_samples/*/*.cif")):
        a = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False)
        r = m.compute(a); r.update(arm=arm, design=os.path.basename(f), length=len(set(a.res_id[~a.hetero]))); rows.append(r)
df = pd.DataFrame(rows); df.to_csv(f"{H}/out/20261121_eval_rfd4.csv", index=False)
pd.set_option("display.width", 200)
print(df.groupby("arm")[["pocket_volume", "volume_ratio", "volume_ratio_error", "ligand_sasa_frac", "loop_frac_shell", "length"]].agg(["mean", "std", "count"]).round(2).to_string())
print(df[["arm","design","length","volume_ratio","loop_frac_shell"]].round(2).to_string())
