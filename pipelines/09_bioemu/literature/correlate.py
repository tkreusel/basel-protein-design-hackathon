"""Correlate per-enzyme metrics with literature log10(kcat/KM).

Labels: /mnt/work/manually_created_dataset.csv (tab separated; kcat/km in M^-1 s^-1).
Metric sources: (a) BioEmu ensemble metrics (literature/metrics/*.json), (b) Boltz-2 / ipSAE2 tables from the
ipsae_scoring_results zip (mean over the 10 models per enzyme).
Sets: HG3-family only (12 enzymes, one scaffold) and all 13 incl. 1A53-core (different scaffold; scaffold is a confound).
Spearman rho with a two-sided permutation-free p (scipy), Benjamini-Hochberg q over ALL tested metrics (exploratory:
~35 metrics x 2 sets, n = 12-13, so |rho| ~ 0.6 is needed even for nominal p < 0.05).
Usage: python -I correlate.py   (kemp-orch env; writes literature/correlations.csv)"""
import csv, json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
Z = HERE.parent / "new_zip/ipsae/ipsae_scoring_results/tables"
rows = list(csv.reader(open("/mnt/work/manually_created_dataset.csv"), delimiter="\t"))[1:]
LAB = {r[1]: float(r[2]) for r in rows}
HG3 = [n for n in LAB if n.startswith("HG3")]
ALL = HG3 + ["1A53-core"]

M = {}                                           # metric -> {enzyme: value}
for p in sorted((HERE / "metrics").glob("*.json")):
    j = json.load(open(p))
    for k, v in j.items():
        if isinstance(v, (int, float)) and k not in ("n", "n_res", "base", "donor"):
            M.setdefault("bioemu." + k, {})[j["name"]] = v
ips = list(csv.DictReader(open(Z / "metrics_model.csv")))
for k in ["confidence_score", "ptm", "iptm", "ligand_iptm", "plddt", "iplddt", "pde", "ipde", "ipsae2_max", "ipsae2_mean",
          "ipsae2_lig_max", "affinity_pred_value", "affinity_probability_binary", "pic50"]:
    d = {}
    for n in ALL:
        v = [float(r[k]) for r in ips if r["sample"] == n and r.get(k) not in ("", None)]
        if v: d[n] = float(np.mean(v))
    if d: M["boltz." + k] = d
        # model-to-model sd = how reproducible the 10 Boltz models are (consistency-like)
    d2 = {}
    for n in ALL:
        v = [float(r[k]) for r in ips if r["sample"] == n and r.get(k) not in ("", None)]
        if len(v) > 1: d2[n] = float(np.std(v))
    if d2: M["boltz.sd_" + k] = d2

out = []
for setname, names in (("HG3-family n=12", HG3), ("all n=13", ALL)):
    for m, d in M.items():
        if m.startswith("bioemu.") and setname.startswith("all"):
            continue      # BioEmu site metrics use HG3 numbering (Asp129): meaningless for the 1A53 scaffold
        ns = [n for n in names if n in d and np.isfinite(d[n])]
        if len(ns) < 8:
            continue
        x = np.array([d[n] for n in ns]); y = np.log10([LAB[n] for n in ns])
        if np.ptp(x) == 0:
            continue
        rho, p = spearmanr(x, y)
        out.append(dict(set=setname, metric=m, n=len(ns), rho=float(rho), p=float(p)))
# selection-aware p: permute the labels, take the max |rho| over the whole metric family (the best-of-many correction)
from scipy.stats import rankdata
rng = np.random.default_rng(0)
def family_perm(setname, names, prefix, nperm=20000):
    ms = [m for m in M if m.startswith(prefix) and all(n in M[m] and np.isfinite(M[m][n]) for n in names)
          and np.ptp([M[m][n] for n in names]) > 0 and not (prefix == "bioemu." and setname.startswith("all"))]
    if len(ms) < 2: return {}
    Xr = np.array([rankdata([M[m][n] for n in names]) for m in ms]); Xr = (Xr - Xr.mean(1, keepdims=True)); Xr /= np.linalg.norm(Xr, axis=1, keepdims=True)
    y = rankdata([np.log10(LAB[n]) for n in names]); 
    P = np.array([rng.permutation(y) for _ in range(nperm)]); P = P - P.mean(1, keepdims=True); P /= np.linalg.norm(P, axis=1, keepdims=True)
    R = np.abs(Xr @ P.T); mx = R.max(0)
    yo = (y - y.mean()) / np.linalg.norm(y - y.mean()); obs = np.abs(Xr @ yo)
    return {m: float((mx >= o - 1e-12).mean()) for m, o in zip(ms, obs)} | {"__family_size": len(ms)}
for setname, names in (("HG3-family n=12", HG3), ("all n=13", ALL)):
    for prefix in ("bioemu.", "boltz."):
        fp = family_perm(setname, names, prefix)
        for o in out:
            if o["set"] == setname and o["metric"] in fp: o["p_selection_aware"] = fp[o["metric"]]; o["family_size"] = fp["__family_size"]
# BH over everything
ps = np.array([o["p"] for o in out]); order = np.argsort(ps); q = np.empty(len(ps))
prev = 1.0
for rank, i in reversed(list(enumerate(order, 1))):
    prev = min(prev, ps[i] * len(ps) / rank); q[i] = prev
for o, qq in zip(out, q): o["q_BH"] = float(qq)
out.sort(key=lambda o: o["p"])
with open(HERE / "correlations.csv", "w") as fh:
    w = csv.DictWriter(fh, fieldnames=["set", "metric", "n", "rho", "p", "q_BH", "p_selection_aware", "family_size"], extrasaction="ignore"); w.writeheader(); w.writerows(out)
print(f"{len(out)} tests ({len(M)} metrics)")
for o in out[:25]:
    print(f"{o['set']:16} {o['metric']:42} n={o['n']:2} rho={o['rho']:+.2f} p={o['p']:.3f} q={o['q_BH']:.2f} p_sel={o.get('p_selection_aware', float('nan')):.3f}")
