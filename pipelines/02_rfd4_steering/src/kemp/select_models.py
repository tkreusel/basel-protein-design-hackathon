"""Stage 5 (CPU): curate the Boltz forest that the evaluation suite reads, and write the per-sequence report.

usage: select_models.py --boltz-dirs D1 D2 ... --out DIR --seq sequences.csv --designs gate.csv [--ids survivors.txt] [--rank-key pocket_rmsd]

All samples of a (sequence, ligand) pair are pooled across runs (1x + 9x TSA; 10x substrate; 10x product). For every sample:
  scrmsd       CA RMSD to the RFD4 design (global CA superposition)
  pocket_rmsd  backbone (N, CA, C, O) RMSD of the pocket residues (any heavy atom within 8 A of the design's 6NT) after superposing on them
plus the Boltz confidence scalars. Samples are ranked per pair (--rank-key; pocket_rmsd lowest first by default, ligand_iptm highest first);
rank 0 is `model_0`, which PLACER, tmol and PyRosetta use (ipSAE2 uses all samples). Symlinks DIR/<seq>__<LIG>/..._model_<rank>.*.

Outputs: ranking.csv (one row per sample), boltz_model_summary.csv (per pair: count/mean/sd/min/max), stage5_report.csv
(one row per sequence: per ligand mean and sd of the scalars, plus TSA-minus-substrate and TSA-minus-product differences of the means).
The report is descriptive; nothing is ranked by it.
"""
import argparse
import glob
import json
import os
import re

os.environ.setdefault("CLUSTER", "local")
import numpy as np
import pandas as pd

SCALARS = ("scrmsd", "pocket_rmsd", "iptm", "ligand_iptm", "complex_ipde", "complex_plddt", "complex_iplddt")
SUMMARY_COLUMNS = ("scrmsd", "pocket_rmsd", "confidence_score", "iptm", "ligand_iptm", "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde", "ptm")
BACKBONE = ("N", "CA", "C", "O")
POCKET_A = 8.0


class DesignReference:
    """CA atoms and pocket backbone of each RFD4 design, loaded on first use."""

    def __init__(self, designs):
        self.paths, self._cache = designs["cif"].to_dict(), {}

    def get(self, design):
        if design not in self._cache:
            from scipy.spatial import cKDTree
            from . import structure
            atoms = structure.load_cif(self.paths[design])
            prot, lig = atoms[~atoms.hetero], atoms[atoms.hetero]
            near = np.unique(prot.res_id[cKDTree(lig.coord).query(prot.coord)[0] <= POCKET_A])
            self._cache[design] = (prot[prot.atom_name == "CA"], prot[np.isin(prot.atom_name, BACKBONE) & np.isin(prot.res_id, near)])
        return self._cache[design]


def sample_rmsds(cif, ca_ref, bb_ref):
    import biotite.structure as bs
    from . import structure
    prot = structure.load_cif(cif)
    prot = prot[~prot.hetero]
    ca = prot[prot.atom_name == "CA"]
    bb = prot[np.isin(prot.atom_name, BACKBONE) & np.isin(prot.res_id, np.unique(bb_ref.res_id))]
    out = {"scrmsd": np.nan, "pocket_rmsd": np.nan}
    if len(ca) == len(ca_ref):
        sup, _ = bs.superimpose(ca_ref, ca)
        out["scrmsd"] = float(bs.rmsd(ca_ref, sup))
    if len(bb) == len(bb_ref) and len(bb_ref):
        sup, _ = bs.superimpose(bb_ref, bb)
        out["pocket_rmsd"] = float(bs.rmsd(bb_ref, sup))
    return out


def rank_key(rank_by):
    higher_better = rank_by in ("ligand_iptm", "confidence_score")
    def key(sample):
        v = sample.get(rank_by)
        if higher_better:
            return -(v or 0)
        return v if v is not None and v == v else 1e9
    return key


def stage5_report(ranking, seqs):
    """One row per sequence from the per-sample table (columns seq_id, ligand, scalars)."""
    rows = []
    for sid, g in ranking.groupby("seq_id"):
        r = dict(seq_id=sid, design=seqs.loc[sid, "design"], variant=seqs.loc[sid, "variant"])
        for lig, h in g.groupby("ligand"):
            r[f"{lig}.n"] = len(h)
            for c in SCALARS:
                if c in h:
                    r[f"{lig}.{c}.mean"], r[f"{lig}.{c}.sd"] = h[c].mean(), h[c].std()
            r[f"{lig}.best_pocket_rmsd"], r[f"{lig}.best_scrmsd"] = h.pocket_rmsd.min(), h.scrmsd.min()
        for other, tag in (("H5J", "sub"), ("6VP", "prod")):
            for c in ("iptm", "ligand_iptm", "complex_ipde", "scrmsd", "pocket_rmsd"):
                if f"6NT.{c}.mean" in r and f"{other}.{c}.mean" in r:
                    r[f"delta_{c}.TSA_minus_{tag}"] = r[f"6NT.{c}.mean"] - r[f"{other}.{c}.mean"]
        rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--boltz-dirs", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seq", required=True)
    ap.add_argument("--designs", required=True)
    ap.add_argument("--ids")
    ap.add_argument("--rank-key", default="pocket_rmsd")
    a = ap.parse_args()

    keep = set(open(a.ids).read().split()) if a.ids else None
    seqs = pd.read_csv(a.seq).set_index("seq_id")
    refs = DesignReference(pd.read_csv(a.designs).drop_duplicates("design").set_index("design"))

    pool = {}
    for root in a.boltz_dirs:
        for cj in glob.glob(root + "/**/predictions/*/confidence_*_model_*.json", recursive=True):
            d = os.path.dirname(cj)
            m = re.match(r"confidence_(?P<name>.+)_model_(?P<k>\d+)\.json$", os.path.basename(cj))
            name, k = m["name"], int(m["k"])
            sid = name.rsplit("__", 1)[0]
            if (keep and sid not in keep) or sid not in seqs.index or not os.path.exists(f"{d}/{name}_model_{k}.cif"):
                continue
            conf = json.load(open(cj))
            sample = dict(src_dir=d, k=k, seq_id=sid, ligand=name.rsplit("__", 1)[1], **{kk: v for kk, v in conf.items() if isinstance(v, (int, float))})
            try:
                sample.update(sample_rmsds(f"{d}/{name}_model_{k}.cif", *refs.get(seqs.loc[sid, "design"])))
            except Exception as exc:
                print("rmsd fail", name, repr(exc)[:100])
            pool.setdefault(name, []).append(sample)

    rows = []
    for name, samples in pool.items():
        samples.sort(key=rank_key(a.rank_key))
        os.makedirs(f"{a.out}/{name}", exist_ok=True)
        for rank, x in enumerate(samples):
            for prefix, ext in (("", "cif"), ("confidence_", "json"), ("pae_", "npz"), ("plddt_", "npz")):
                src = f"{x['src_dir']}/{prefix}{name}_model_{x['k']}.{ext}"
                dst = f"{a.out}/{name}/{prefix}{name}_model_{rank}.{ext}"
                if os.path.exists(src) and not os.path.lexists(dst):
                    os.symlink(src, dst)
            rows.append(dict(name=name, rank=rank, src_k=x["k"], **{k: v for k, v in x.items() if k != "k"}))
    ranking = pd.DataFrame(rows)
    ranking.to_csv(a.out + "/ranking.csv", index=False)
    cols = [c for c in SUMMARY_COLUMNS if c in ranking]
    summary = ranking.groupby("name")[cols].agg(["count", "mean", "std", "min", "max"])
    summary.columns = [f"{c}|{s}" for c, s in summary.columns]
    summary.to_csv(a.out + "/boltz_model_summary.csv")
    report = stage5_report(ranking, seqs)
    report.to_csv(a.out + "/stage5_report.csv", index=False)
    print(len(pool), "protein-ligand sets,", len(ranking), "samples,", len(report), "sequences -> stage5_report.csv")


if __name__ == "__main__":
    main()
