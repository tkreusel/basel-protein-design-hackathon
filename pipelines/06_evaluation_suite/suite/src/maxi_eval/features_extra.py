"""Extra features mined from the RAW outputs (per-model PLACER geometry, per-residue tmol energies, Boltz PAE/pLDDT matrices, sequences).

  python -m maxi_eval.features_extra --out work/features_extra.csv

Families (column prefixes):
  pl2.<ens>.*   PLACER ensemble DISTRIBUTION features, ens in {tsa_holo, tsa_apo, sub_holo, prod_holo, apoarm_apo}:
                near-attack-conformation (NAC) fractions at several cut-offs (scaffold-independent, unlike the HG3 pass windows), quantiles, best-prmsd-10% geometry.
  pr.<tag>.*    tmol per-residue ligand-interaction energies for tag in {tsa_pk, sub, prod}: base / donor / wall contributions, residue-class sums
                (charged, polar, aromatic, hydrophobic), concentration of the binding energy (top-1, top-3 residues), residue counts, H-bond partners.
  bz.<arm>.*    Boltz PAE/pLDDT pocket features for arm in {tsa, sub, prod, apo}: pocket pLDDT, pocket-internal PAE (rigidity), pocket-to-rest PAE,
                catalytic-residue PAE, ligand<->pocket PAE.
  cmp.*         composition of the 6 A pocket, sequence length, mutations vs HG3 (HG3 family only: a confound/lineage check).
Each family also gets the arm DIFFERENCES (tsa-sub, tsa-prod, sub-prod, and apo-tsa for bz/pl2) as columns named <prefix>.<feature>.<a>-<b>.
"""
from __future__ import annotations
import argparse, glob, json, os, re, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")

CLS = {"pos": set("KR"), "neg": set("DE"), "polar": set("STNQHC"), "arom": set("FWYH"), "hydro": set("AVILMPG")}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
       "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def protein(entry: str) -> str:
    return entry.replace("_boltz", "")


def seq_of(pdb: str) -> dict[int, str]:
    out = {}
    for l in open(pdb):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            out[int(l[22:26])] = AA3.get(l[17:20].strip(), "X")
    return out


# ---------------------------------------------------------------- PLACER ensembles
def ensemble_features(df: pd.DataFrame) -> dict:
    f = {}
    d = df.get("geo.d_O_X")
    if d is not None and d.notna().sum() > 3:
        f.update(dOX_med=d.median(), dOX_p10=d.quantile(0.1), dOX_min=d.min(), dOX_iqr=d.quantile(0.75) - d.quantile(0.25),
                 dOX_f30=(d <= 3.0).mean(), dOX_f35=(d <= 3.5).mean(), dOX_f40=(d <= 4.0).mean())
        ang = df.get("geo.ang_O_H_X")
        if ang is not None:
            f["nac35"] = ((d <= 3.5) & (ang >= 110)).mean(); f["nac30"] = ((d <= 3.0) & (ang >= 130)).mean()
            f["angOHX_med"] = ang.median()
        if "prmsd" in df:
            top = df.nsmallest(max(2, len(df) // 10), "prmsd"); f["top10_dOX"] = top["geo.d_O_X"].mean()
    for c, k in (("geo.ang_C_O_H", "angCOH_med"), ("geo.donor_best_d", "donor_med"), ("geo.wall_d_centroid", "wall_d_med"), ("geo.contact_min", "contact_min_med"),
                 ("geo.n_contacts", "n_contacts_mean"), ("geo.lig_rel_sasa", "lig_sasa_mean"), ("geo.enclosure", "enclosure_mean"), ("geo.base_sasa_free_frac", "base_sasa_mean"),
                 ("geo.n_shell", "n_shell_mean")):
        if c in df and df[c].notna().any():
            f[k] = df[c].median() if k.endswith("med") else df[c].mean()
    if "geo.donor_best_d" in df and df["geo.donor_best_d"].notna().sum() > 3:
        f["donor_f35"] = (df["geo.donor_best_d"] <= 3.5).mean()
    if "prmsd" in df:
        f.update(prmsd_min=df.prmsd.min(), prmsd_q25=df.prmsd.quantile(0.25), prmsd_sd=df.prmsd.std())
    if "plddt" in df:
        f.update(plddt_max=df.plddt.max(), plddt_med=df.plddt.median())
    for c, k in (("field_CH", "fieldCH_mean"), ("field_LG", "fieldLG_mean")):
        if c in df and df[c].notna().sum() > 2:
            f[k] = df[c].mean(); f[k.replace("mean", "sd")] = df[c].std()
    return f


def placer_block(man) -> pd.DataFrame:
    rows = {}
    for entry, e in man["entries"].items():
        if not entry.endswith("_boltz"):
            continue
        acc: dict[str, list] = {}
        for r in e["replicates"]:
            base = f"work/placer/{entry}/{r['rep_id'].replace(':', '_')}"
            for mode, ens in (("holo", f"{r['arm']}_holo"), ("apo", "tsa_apo" if r["arm"] == "tsa" else "apoarm_apo")):
                if (mode == "holo" and r["arm"] == "apo") or (mode == "apo" and r["arm"] in ("sub", "prod")):
                    continue
                fn = f"{base}_{mode}_models.csv"
                if Path(fn).exists():
                    for k, v in ensemble_features(pd.read_csv(fn)).items():
                        acc.setdefault(f"pl2.{ens}.{k}", []).append(v)
        rows[protein(entry)] = {k: float(np.nanmean(v)) for k, v in acc.items()}
    return pd.DataFrame(rows).T


# ---------------------------------------------------------------- tmol per-residue energies
def residue_features(det: dict, seq: dict[int, str], site: dict) -> dict:
    ip = {int(k): v for k, v in det.get("int_per_res", {}).items()}
    if not ip:
        return {}
    f = {}
    tot = sum(ip.values()); f["int_sum"] = tot
    vals = np.array(sorted(ip.values())); f["int_top1"] = vals[0]; f["int_top3"] = vals[:3].sum(); f["int_top1_frac"] = vals[0] / tot if tot else np.nan
    f["n_res_lt05"] = float((vals < -0.5).sum()); f["n_res_lt1"] = float((vals < -1.0).sum()); f["n_res_gt1"] = float((vals > 1.0).sum()); f["repulsive_sum"] = float(vals[vals > 0].sum())
    b = site.get("base"); w = site.get("wall"); dn = [d for d in (site.get("donors") or []) if not d.get("backbone")]
    if b:
        f["int_base"] = ip.get(b["resid"], np.nan)
    if w:
        f["int_wall"] = ip.get(w["resid"], np.nan)
    if dn:
        f["int_donor"] = ip.get(dn[0]["resid"], np.nan)
    for cls, aas in CLS.items():
        f[f"int_{cls}"] = sum(v for r, v in ip.items() if seq.get(r) in aas)
    hb = det.get("int_hbond_per_res", {})
    f["n_hbond_res"] = float(len(hb)); f["hbond_sum"] = float(sum(hb.values()))
    pe = det.get("per_res_energy", {})
    pocket = site.get("pocket_resids") or []
    pv = np.array([pe[str(r)] for r in pocket if str(r) in pe])
    if len(pv):
        f["pocket_E_sum"] = pv.sum(); f["pocket_E_max"] = pv.max(); f["pocket_n_gt2"] = float((pv > 2).sum())
    return f


def tmol_block(man) -> pd.DataFrame:
    rows = {}
    for entry, e in man["entries"].items():
        if not entry.endswith("_boltz"):
            continue
        tsa = [r for r in e["replicates"] if r["arm"] == "tsa"]
        acc: dict[str, list] = {}
        for r in tsa:
            fn = f"{os.environ.get('TMOL_DIR', 'work/tmol')}/{r['rep_id'].replace(':', '_')}.json"
            if not Path(fn).exists():
                continue
            d = json.load(open(fn)); det = d.get("detail", {})
            if "site" not in det:
                continue
            seq = seq_of(r["pdb"]); site = det["site"]
            for tag in ("tsa_pk", "sub", "prod"):
                if tag in det:
                    for k, v in residue_features(det[tag], seq, site).items():
                        acc.setdefault(f"pr.{tag}.{k}", []).append(v)
        rows[protein(entry)] = {k: float(np.nanmean(v)) for k, v in acc.items()}
    return pd.DataFrame(rows).T


# ---------------------------------------------------------------- Boltz PAE / pLDDT
def blk(P, a, b):
    if len(a) == 0 or len(b) == 0:
        return np.nan
    m = P[np.ix_(a, b)]
    return float(m.mean()) if m.size else np.nan


def bz_features(pae, plddt, n_res, pocket, cat, has_lig) -> dict:
    P = np.asarray(pae); f = {}
    pk = np.array([r - 1 for r in pocket if 0 < r <= n_res]); ca = np.array([r - 1 for r in cat if 0 < r <= n_res])
    rest = np.setdiff1d(np.arange(n_res), pk)
    if plddt is not None:
        pl = np.asarray(plddt)
        f.update(protein_plddt=float(pl[:n_res].mean()), pocket_plddt_mean=float(pl[pk].mean()), pocket_plddt_min=float(pl[pk].min()), cat_plddt_mean=float(pl[ca].mean()) if len(ca) else np.nan)
    sub = P[np.ix_(pk, pk)]; f["pocket_pae_int"] = float((sub.sum() - np.trace(sub)) / max(sub.size - len(pk), 1))
    f["pocket_to_rest_pae"] = blk(P, pk, rest); f["rest_to_pocket_pae"] = blk(P, rest, pk)
    f["cat_pae_int"] = blk(P, ca, ca); f["cat_to_pocket_pae"] = blk(P, ca, pk)
    if has_lig and P.shape[0] > n_res:
        lg = np.arange(n_res, P.shape[0])
        f["lig_pocket_pae"] = blk(P, lg, pk); f["pocket_lig_pae"] = blk(P, pk, lg); f["lig_cat_pae"] = blk(P, lg, ca); f["cat_lig_pae"] = blk(P, ca, lg)
        if plddt is not None:
            f["lig_plddt"] = float(np.asarray(plddt)[lg].mean())
    return f


def boltz_block(man) -> pd.DataFrame:
    rows = {}
    for entry, e in man["entries"].items():
        if not entry.endswith("_boltz"):
            continue
        tsa = next((r for r in e["replicates"] if r["arm"] == "tsa" and r.get("site")), None)
        if tsa is None:
            continue
        s = tsa["site"]; pocket = s["pocket_resids"]
        cat = [x["resid"] for x in (s.get("base"), s.get("wall")) if x] + [d["resid"] for d in s["donors"] if not d["backbone"]][:1]
        acc: dict[str, list] = {}
        for r in e["replicates"]:
            if not r.get("boltz_pae"):
                continue
            z = np.load(r["boltz_pae"]); pae = z[z.files[0]]
            pl = None
            if r.get("boltz_plddt") and Path(r["boltz_plddt"]).exists():
                q = np.load(r["boltz_plddt"]); pl = q[q.files[0]]
            for k, v in bz_features(pae, pl, r["n_res"], pocket, cat, r["arm"] != "apo").items():
                acc.setdefault(f"bz.{r['arm']}.{k}", []).append(v)
        rows[protein(entry)] = {k: float(np.nanmean(v)) for k, v in acc.items()}
    return pd.DataFrame(rows).T


# ---------------------------------------------------------------- composition / lineage
def composition_block(man) -> pd.DataFrame:
    rows, seqs = {}, {}
    for entry, e in man["entries"].items():
        if not entry.endswith("_boltz"):
            continue
        tsa = next((r for r in e["replicates"] if r["arm"] == "tsa" and r.get("site")), None)
        if tsa is None:
            continue
        sq = seq_of(tsa["pdb"]); seqs[protein(entry)] = "".join(sq[i] for i in sorted(sq))
        pk = [sq[r] for r in tsa["site"]["pocket_resids"] if r in sq]; n = max(len(pk), 1)
        row = {"cmp.length": float(len(sq)), "cmp.n_pocket": float(len(pk)), "cmp.pocket_charge": float(sum(aa in "KR" for aa in pk) - sum(aa in "DE" for aa in pk))}
        for cls, aas in CLS.items():
            row[f"cmp.pocket_frac_{cls}"] = sum(aa in aas for aa in pk) / n
        row["cmp.pocket_n_arom"] = float(sum(aa in "FWY" for aa in pk)); row["cmp.pocket_n_gly"] = float(sum(aa == "G" for aa in pk))
        rows[protein(entry)] = row
    hg3 = seqs.get("HG3")
    for p, sq in seqs.items():
        if p.startswith("HG3") and hg3 and len(sq) == len(hg3):
            rows[p]["cmp.n_mut_vs_HG3"] = float(sum(a != b for a, b in zip(sq, hg3)))
    return pd.DataFrame(rows).T


# ---------------------------------------------------------------- differences
def add_diffs(F: pd.DataFrame) -> pd.DataFrame:
    new = {}
    def feats(prefix_fmt, arms):
        sets = [{c[len(prefix_fmt.format(a)) :] for c in F.columns if c.startswith(prefix_fmt.format(a))} for a in arms]
        return set.intersection(*sets) if sets else set()
    # pr.<tag>.  tags: tsa_pk, sub, prod
    for pa, pb in (("tsa_pk", "sub"), ("tsa_pk", "prod"), ("sub", "prod")):
        for s in feats("pr.{}.", (pa, pb)):
            new[f"pr.{pa}-{pb}.{s}"] = F[f"pr.{pa}.{s}"] - F[f"pr.{pb}.{s}"]
    for pa, pb in (("tsa", "sub"), ("tsa", "prod"), ("sub", "prod"), ("apo", "tsa"), ("apo", "sub"), ("apo", "prod")):
        for s in feats("bz.{}.", (pa, pb)):
            new[f"bz.{pa}-{pb}.{s}"] = F[f"bz.{pa}.{s}"] - F[f"bz.{pb}.{s}"]
    for pa, pb in (("tsa_holo", "sub_holo"), ("tsa_holo", "prod_holo"), ("sub_holo", "prod_holo"), ("tsa_apo", "tsa_holo"), ("apoarm_apo", "tsa_holo"), ("apoarm_apo", "tsa_apo")):
        for s in feats("pl2.{}.", (pa, pb)):
            new[f"pl2.{pa}-{pb}.{s}"] = F[f"pl2.{pa}.{s}"] - F[f"pl2.{pb}.{s}"]
    return pd.concat([F, pd.DataFrame(new, index=F.index)], axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs"); ap.add_argument("--out", default="work/features_extra.csv")
    a = ap.parse_args()
    man = json.loads((Path(a.inputs) / "manifest.json").read_text())
    blocks = [placer_block(man), tmol_block(man), boltz_block(man), composition_block(man)]
    F = pd.concat(blocks, axis=1)
    F = add_diffs(F)
    F = F.loc[:, F.notna().sum() >= min(12, max(len(F) - 1, 1))]   # (was a fixed 12: dropped everything for small batches)
    F.to_csv(a.out)
    pref = pd.Series([c.split(".")[0] for c in F.columns]).value_counts().to_dict()
    print(f"{F.shape[0]} proteins x {F.shape[1]} extra features  {pref} -> {a.out}")

if __name__ == "__main__":
    main()
