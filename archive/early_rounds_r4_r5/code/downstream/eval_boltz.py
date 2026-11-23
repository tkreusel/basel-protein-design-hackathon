"""Evaluate Boltz-2 refolds against the RFD4 designs -> runs/20261121_r5_downstream/boltz_eval.csv (GPU job for TMol; run in the rfd4 venv).
Per prediction: Boltz confidence (complex pLDDT, pTM, ipTM, ligand ipTM; protein/ligand pLDDT from CIF B-factors), CA RMSD to the design (Kabsch on all CA),
ligand heavy-atom RMSD after CA superposition, motif-tip RMSD after CA superposition (fixed motif residues), plus the same chemistry metrics as the designs
(pocket ratio, ligand SASA, base SASA/polar/cation, nitro H-bonds, channel bottleneck, TMol dG) computed on the PREDICTED complex.
usage: python -I eval_boltz.py"""
import glob, json, os, numpy as np, pandas as pd, biotite.structure as s, biotite.structure.io.pdbx as pdbx
from rfproteina.metrics.custom_metrics.kemp_base import KempBaseMetric
from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
from rfproteina.metrics.custom_metrics.kemp_channel import KempChannelMetric
from rfproteina.metrics.custom_metrics.kemp_tmol import KempTMolMetric
H = os.path.expanduser("~/HACKATHON_ETH"); D = H + "/runs/20261121_r5_downstream"
mb, mp, mc = KempBaseMetric(), KempPocketMetric(target_ratio=0.5, tol=0.10), KempChannelMetric()
tm = KempTMolMetric(region_a="is_ligand", region_b="is_protein", ligand_params_files=[H + "/inputs/6NT.tmol"])
seqs = pd.read_csv(D + "/sequences.csv").set_index("seq_id"); sel = pd.read_csv(D + "/selection.csv").set_index("design")
TIPS = {"ASP": ("CB", "CG", "OD1", "OD2"), "GLN": ("CG", "CD", "OE1", "NE2"), "GLU": ("CG", "CD", "OE1", "OE2"), "SER": ("CA", "CB", "OG"),
        "TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")}
def load(f):
    a = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False, extra_fields=["b_factor"]); return a[(a.element != "H") & ~np.isnan(a.coord).any(1)]
rows = []
for f in sorted(glob.glob(D + "/boltz_out/chunk_*/boltz_results_*/predictions/*/*_model_0.cif")):
    sid = os.path.basename(f)[:-len("_model_0.cif")]; design = seqs.loc[sid, "design"]
    try:
        pr = load(f); de = load(sel.loc[design, "cif"]); cj = json.load(open(f.replace(".cif", ".json").replace(sid + "_model_0", "confidence_" + sid + "_model_0")))
        pp, pl = pr[~pr.hetero], pr[pr.hetero]; dp, dl = de[~de.hetero], de[de.hetero]
        pca, dca = pp[pp.atom_name == "CA"], dp[dp.atom_name == "CA"]
        r = dict(seq_id=sid, design=design, source=seqs.loc[sid, "source"], motif=sel.loc[design, "motif"], arm=sel.loc[design, "arm"], length=len(pca),
                 complex_plddt=cj.get("complex_plddt"), ptm=cj.get("ptm"), iptm=cj.get("iptm"), ligand_iptm=cj.get("ligand_iptm"), confidence=cj.get("confidence_score"),
                 protein_plddt=float(np.mean(pp.b_factor)), ligand_plddt=float(np.mean(pl.b_factor)))
        if len(pca) == len(dca):
            _, tr = s.superimpose(dca, pca); r["ca_rmsd"] = float(s.rmsd(dca, tr.apply(pca)) if hasattr(tr, "apply") else np.nan)
            ligmap = [n for n in dl.atom_name if n in set(pl.atom_name)]
            X = np.array([pl.coord[pl.atom_name == n][0] for n in ligmap]); Y = np.array([dl.coord[dl.atom_name == n][0] for n in ligmap])
            r["ligand_rmsd_ca_aligned"] = float(np.sqrt(((tr.apply(X) - Y) ** 2).sum(1).mean()))
            acc = []
            for tok in str(seqs.loc[sid, "fixed"]).split():
                rid = int(tok[1:]); q = dp[dp.res_id == rid]; w = pp[pp.res_id == rid]; names = TIPS.get(q.res_name[0], ())
                if names and all(n in set(w.atom_name) for n in names):
                    A = np.array([q.coord[q.atom_name == n][0] for n in names]); B = tr.apply(np.array([w.coord[w.atom_name == n][0] for n in names]))
                    acc.append(np.sqrt(((A - B) ** 2).sum(1).mean()))
            r["motif_tip_rmsd_ca_aligned"] = float(np.max(acc)) if acc else np.nan
        r.update(**{k: v for k, v in {**mp.compute(pr), **mb.compute(pr), **mc.compute(pr)}.items()})
        try: r["tmol_ddg"] = tm.compute(pr)["interface_ddg"]
        except Exception: r["tmol_ddg"] = np.nan
        rows.append(r)
    except Exception as e:
        print("fail", sid, repr(e)[:150], flush=True)
pd.DataFrame(rows).to_csv(D + "/boltz_eval.csv", index=False); print(len(rows), "predictions evaluated")
