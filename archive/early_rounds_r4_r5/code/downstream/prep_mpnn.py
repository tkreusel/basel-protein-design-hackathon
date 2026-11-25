"""LigandMPNN (ligand_mpnn v_32_010_25, T=0.1, 4 sequences) for each selected design, motif residues fixed, polar/charged residues forbidden near the base.
Per design: PDB (chain A protein, chain B 6NT) -> fixed_residues = located motif residues; omit KRHDENQSTYC at non-motif residues within 6 A of the base carboxylate O.
usage: ~/HACKATHON_ETH/rfd4/.venv/bin/python -I prep_mpnn.py [--limit N]   (calls the ligandmpnn_env python)
Outputs: runs/20261121_r5_downstream/mpnn/<design>/{in.pdb, omit.json, seqs/*.fa}; sequences.csv (design, seq_id, source, sequence)."""
import glob, json, os, subprocess, sys, numpy as np, pandas as pd, biotite.structure as s, biotite.structure.io.pdb as pdb, biotite.structure.io.pdbx as pdbx
from scipy.spatial import cKDTree
from rfproteina.metrics.custom_metrics.kemp_chem import KempChemMetric
H = os.path.expanduser("~/HACKATHON_ETH"); D = H + "/runs/20261121_r5_downstream"; LM = os.path.expanduser("~/Workspace_Maxi/MLfBW/LigandMPNN")
PY = os.path.expanduser("~/.conda/envs/ligandmpnn_env/bin/python")
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
OMIT = "KRHDENQSTYC"
chg = KempChemMetric(ref_pdb=H + "/inputs/kemp_theozyme_GDM_KE_1483.pdb")
chh = KempChemMetric(ref_pdb=H + "/inputs/HG3_17_kemp_motif.pdb", base_res="ASP", base_match=("CB", "CG", "OD1", "OD2"), base_o=("OD1", "OD2"), hb_res="GLN", hb_match=("CG", "CD", "OE1", "NE2"), hb_atom="NE2")
ref_h = pdb.PDBFile.read(H + "/inputs/HG3_17_kemp_motif.pdb").get_structure(model=1); ref_h = ref_h[ref_h.element != "H"]
sel = pd.read_csv(D + "/selection.csv"); lim = int(sys.argv[sys.argv.index("--limit") + 1]) if "--limit" in sys.argv else None
rows = []
for _, r in (sel.head(lim) if lim else sel).iterrows():
    out = f"{D}/mpnn/{r.design}"; os.makedirs(out, exist_ok=True)
    a = pdbx.get_structure(pdbx.CIFFile.read(r.cif), model=1, use_author_fields=False); a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]
    lig, p = a[a.hetero], a[~a.hetero]
    # chain A protein with sequential numbering, ligand chain B
    p.chain_id[:] = "A"; lig.chain_id[:] = "B"; lig.res_id[:] = 1
    pdb_f = pdb.PDBFile(); pdb_f.set_structure(s.concatenate([p, lig])); pdb_f.write(f"{out}/in.pdb")
    hg = "hg" in r.motif; ch = chh if hg else chg
    rl = ch.ref[ch.ref.hetero]; names = [n for n in rl.atom_name if n in set(lig.atom_name)]
    _, tr = s.superimpose(np.array([lig.coord[lig.atom_name == n][0] for n in names]), np.array([rl.coord[rl.atom_name == n][0] for n in names]))
    fixed = []
    spec = ([("ASP", ("CB", "CG", "OD1", "OD2")), ("GLN", ("CG", "CD", "OE1", "NE2"))] + ([("TRP", ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"))] if "with-trp" in r.motif else [])) if hg \
        else [("GLU", ("CG", "CD", "OE1", "OE2")), ("SER", ("CA", "CB", "OG"))]
    for res, match in spec:
        rr = ch.ref[(ch.ref.res_name == res) & ~ch.ref.hetero] if not (hg and res == "TRP") else ref_h[(ref_h.res_name == "TRP") & ~ref_h.hetero]
        X = tr.apply(np.array([rr.coord[rr.atom_name == n][0] for n in match])); best = (9, None)
        for rid in np.unique(p.res_id[p.res_name == res]):
            q = p[(p.res_id == rid) & (p.res_name == res)]
            try: Y = np.array([q.coord[q.atom_name == n][0] for n in match])
            except IndexError: continue
            dd = float(np.sqrt(((X - Y) ** 2).sum(1).mean()))
            if dd < best[0]: best = (dd, rid)
        fixed.append((res, best[1], best[0]))
    fixed_ids = [x[1] for x in fixed if x[1] is not None]
    base_res = "ASP" if hg else "GLU"; brid = [x[1] for x in fixed if x[0] == base_res][0]
    bo = p[(p.res_id == brid) & np.isin(p.atom_name, ("OD1", "OD2", "OE1", "OE2"))]
    pt = cKDTree(p.coord); near = {int(p.res_id[i]) for i in np.unique(np.concatenate(pt.query_ball_point(bo.coord, 6.0)))} - set(fixed_ids)
    json.dump({f"A{i}": OMIT for i in sorted(near)}, open(f"{out}/omit.json", "w"))
    seq = "".join(AA3.get(n, "X") for n in p[p.atom_name == "CA"].res_name)
    rows.append(dict(design=r.design, seq_id=f"{r.design}__rfd4", source="rfd4", sequence=seq, fixed=" ".join(f"A{i}" for i in fixed_ids), n_omit=len(near)))
    cmd = [PY, f"{LM}/run.py", "--model_type", "ligand_mpnn", "--checkpoint_ligand_mpnn", f"{LM}/model_params/ligandmpnn_v_32_010_25.pt", "--pdb_path", f"{out}/in.pdb",
           "--out_folder", out, "--number_of_batches", "4", "--batch_size", "1", "--temperature", "0.1", "--seed", "111", "--fixed_residues", " ".join(f"A{i}" for i in fixed_ids),
           "--omit_AA_per_residue", f"{out}/omit.json", "--ligand_mpnn_use_atom_context", "1", "--save_stats", "0"]
    if not glob.glob(f"{out}/seqs/*.fa"): subprocess.run(cmd, cwd=LM, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    fa = open(glob.glob(f"{out}/seqs/*.fa")[0]).read().split(">")[1:]
    for k, rec in enumerate(fa[1:], 1):          # record 0 = input sequence
        sq = rec.split("\n", 1)[1].replace("\n", "").split(":")[0] if rec.count("\n") else ""
        rows.append(dict(design=r.design, seq_id=f"{r.design}__mpnn{k}", source=f"mpnn{k}", sequence=sq, fixed=" ".join(f"A{i}" for i in fixed_ids), n_omit=len(near)))
    print(r.design, "fixed", fixed, "omit", len(near), flush=True)
pd.DataFrame(rows).to_csv(D + "/sequences.csv" if not lim else D + "/sequences_test.csv", index=False)
