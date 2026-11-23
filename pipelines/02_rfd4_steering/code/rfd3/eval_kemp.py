"""Post-hoc evaluation of RFD3 trial outputs with the SAME metric as RFD4 (rfd4 KempPocketMetric, full atoms, 0.7 A grid).
usage: /home/hd/hd_hd/hd_la299/HACKATHON_ETH/rfd4/.venv/bin/python -I eval_kemp.py   -> out/20261121_eval.csv + summary
Input: out/20261121_<key>_<arm>/raw.tar.gz (or raw/ dir) with *.cif.gz designs. Only the design matching <key> is scored."""
import glob, gzip, io, os, sys, tarfile, numpy as np, pandas as pd
import biotite.structure.io.pdbx as pdbx
from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
H = os.path.expanduser("~/HACKATHON_ETH/rfd3"); m = KempPocketMetric(); rows = []
for d in sorted(glob.glob(f"{H}/out/20261121_*_*")):
    base = os.path.basename(d)[9:]
    if base.startswith("test") or "_" not in base: continue
    key, arm = base.split("_", 1)
    files = []
    if os.path.exists(f"{d}/raw.tar.gz"):
        with tarfile.open(f"{d}/raw.tar.gz") as t:
            files = [(mem.name, t.extractfile(mem).read()) for mem in t if mem.name.endswith(".cif.gz")]
    else:
        files = [(f, open(f, "rb").read()) for f in glob.glob(f"{d}/**/*.cif.gz", recursive=True)]
    for name, raw in files:
        if f"_{key}_150_260_model" not in name: continue
        a = pdbx.get_structure(pdbx.CIFFile.read(io.StringIO(gzip.decompress(raw).decode())), model=1, use_author_fields=False)
        r = m.compute(a); r.update(key=key, arm=arm, design=os.path.basename(name), length=len(set(a.res_id[~a.hetero])))
        rows.append(r)
df = pd.DataFrame(rows); df.to_csv(f"{H}/out/20261121_eval.csv", index=False)
pd.set_option("display.width", 200)
print(df.groupby(["key", "arm"])[["pocket_volume", "volume_ratio", "volume_ratio_error", "ligand_sasa_frac", "loop_frac_shell", "length"]].agg(["mean", "std", "count"]).round(2).to_string())
