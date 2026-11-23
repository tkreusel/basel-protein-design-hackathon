"""Analyse BioEmu-2 equilibrium ensembles of protein + 6NT (Kemp substrate analogue).

Each run in bioemu/outputs/<name>/ holds samples.xtc (N independent samples, NOT a
time series: there is no time axis, no autocorrelation, "motion" means how the
ensemble spreads) and topology.cif (protein chain A + ligand chain B, heavy atoms).

Everything is expressed in the frame of the designed (RFD4) structure: each sample is
superposed on the design's CA atoms using an iteratively chosen rigid core (residues
whose spread is below 1.5x the median), so RMSF measures internal flexibility and not
overall tumbling. Ligand motion is then the ligand's motion RELATIVE TO the protein core.

Per run it writes bioemu/analysis/<run>/
  frames.csv      one row per sample: CA RMSD, Rg, ligand pose RMSD / COM shift / tilt,
                  host-style approach geometry (the same ten measurements and scores as
                  scoring/host_score.py), base burial, ligand SASA, contacts, PCA scores
  residues.csv    per residue: RMSF (core-aligned), RMSD to design, ligand-contact
                  fraction, helix / strand fraction, motif flag
  summary.json    headline numbers
  ensemble_20frames.pdb, medoid.pdb   for PyMOL / ChimeraX

Usage: python -I analyze_ensembles.py <run name> [<run name> ...]   (kemp-orch env)
"""

import csv
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from biotite.structure import annotate_sse, sasa
from biotite.structure.io.pdb import PDBFile
from biotite.structure.io.pdbx import CIFFile, get_structure
import biotite.structure.io.xtc as xtc
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scoring"))
from host_geometry import BASE_METRICS, DONOR_METRICS, burial_count, geometry  # noqa: E402
import host_score  # noqa: E402  (geom_scores; main is guarded)

HAND = ROOT / "imports/steering_handoff_20261121/20261121_rfd4_steering_handoff/runs"
SEQ_CSV = HAND / "20261121_r5_downstream/sequences.csv"
DESIGNS = HAND / "20261121_r5_downstream/designs_in"
CONTACT_CUTOFF = 4.0
BOUND_RMSD = 3.0     # ligand heavy-atom RMSD to the designed pose (core-aligned) counted as "bound in the pocket"
LIG_ATOMS = ["N3", "C3A", "N2", "C4", "N1", "C7A", "C7", "C6", "C5", "NO1", "O21", "O11"]
RING = ["C3A", "C4", "C5", "C6", "C7", "C7A", "N1", "N2", "N3"]


def kabsch_apply(P, Q, w=None):
    """Rotation/translation taking P onto Q (both (n,3)); returns R, pc, qc with  (x-pc)@R+qc."""
    pc, qc = P.mean(0), Q.mean(0)
    v, _, wt = np.linalg.svd((P - pc).T @ (Q - qc))
    if np.linalg.det(v @ wt) < 0:
        v[:, -1] *= -1
    return v @ wt, pc, qc


def align_all(coords, ref_ca, ca_idx, core):
    out = np.empty_like(coords)
    for i, X in enumerate(coords):
        R, pc, qc = kabsch_apply(X[ca_idx][core], ref_ca[core])
        out[i] = (X - pc) @ R + qc
    return out


def normal(pts):
    c = pts.mean(0)
    n = np.linalg.svd(pts - c)[2][-1]
    return n / np.linalg.norm(n)


def analyse(run: str) -> dict:
    d = HERE / "outputs" / run
    out = HERE / "analysis" / run; out.mkdir(parents=True, exist_ok=True)
    top = get_structure(CIFFile.read(d / "topology.cif"), model=1)
    traj = xtc.XTCFile.read(str(d / "samples.xtc")).get_structure(template=top)
    X = traj.coord.astype(np.float64)                      # (N, atoms, 3) Angstrom
    N = len(X)
    prot = (top.chain_id == "A"); lig = (top.chain_id == "B")
    ca_idx = np.where(prot & (top.atom_name == "CA"))[0]
    n_res = len(ca_idx)
    res_ids = top.res_id[ca_idx]
    # sanity: units are Angstrom (N-CA bond ~1.46)
    n_i = np.where(prot & (top.atom_name == "N"))[0][1]
    ca_i = ca_idx[1]
    bond = float(np.linalg.norm(X[0, n_i] - X[0, ca_i]))
    assert 1.3 < bond < 1.6, f"unexpected N-CA length {bond:.2f}: units?"

    # ---- design reference, base / donor positions ------------------------------
    sr = next(r for r in csv.DictReader(open(SEQ_CSV)) if r["seq_id"] == run)
    fixed = [int(x[1:]) for x in sr["fixed"].split()]
    base_res, donor_res = fixed[0], fixed[1]
    seq3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
            "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
            "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    seq = "".join(seq3[n] for n in top.res_name[ca_idx])
    assert seq == sr["sequence"], "topology sequence differs from sequences.csv"
    assert seq[base_res - 1] == "D" and seq[donor_res - 1] == "Q"
    design_name = run.rsplit("__", 1)[0]
    dstruct = get_structure(CIFFile.read(DESIGNS / f"{design_name}.cif"), model=1)
    d_ca = dstruct[(~dstruct.hetero) & (dstruct.atom_name == "CA")]
    d_lig = dstruct[dstruct.res_name == "6NT"]
    assert len(d_ca) == n_res, f"design has {len(d_ca)} residues, ensemble {n_res}"
    ref_ca = d_ca.coord.astype(np.float64)
    lig_names = list(top.atom_name[lig])
    d_lig_xyz = np.array([d_lig.coord[list(d_lig.atom_name).index(n)] for n in lig_names], dtype=np.float64)

    # ---- iterative rigid core alignment on the design CA ------------------------
    core = np.ones(n_res, bool)
    for _ in range(4):
        Xa = align_all(X, ref_ca, ca_idx, core)
        dev = np.sqrt(((Xa[:, ca_idx] - ref_ca) ** 2).sum(-1).mean(0))     # per-residue RMS dev from design
        core = dev < 1.5 * np.median(dev)
    Xa = align_all(X, ref_ca, ca_idx, core)
    CA = Xa[:, ca_idx]
    mean_ca = CA.mean(0)
    rmsf = np.sqrt(((CA - mean_ca) ** 2).sum(-1).mean(0))              # per residue, core-aligned
    rmsd_design_res = np.sqrt(((CA - ref_ca) ** 2).sum(-1).mean(0))
    ca_rmsd_design_all = np.sqrt(((CA - ref_ca) ** 2).sum(-1).mean(1))
    ca_rmsd_design_core = np.sqrt(((CA[:, core] - ref_ca[core]) ** 2).sum(-1).mean(1))
    ca_rmsd_mean = np.sqrt(((CA - mean_ca) ** 2).sum(-1).mean(1))
    rg = np.sqrt(((CA - CA.mean(1, keepdims=True)) ** 2).sum(-1).mean(1))

    # ---- ligand motion relative to the protein core ------------------------------
    L = Xa[:, lig]                                                     # (N, 12, 3)
    lig_rmsd = np.sqrt(((L - d_lig_xyz) ** 2).sum(-1).mean(1))         # vs designed pose
    com = L.mean(1); com_d = d_lig_xyz.mean(0)
    com_shift = np.linalg.norm(com - com_d, axis=1)
    lig_rmsf = np.sqrt(((L - L.mean(0)) ** 2).sum(-1).mean(0))         # per ligand atom
    ring_ix = [lig_names.index(n) for n in RING]
    n_d = normal(d_lig_xyz[ring_ix])
    tilt = np.array([np.degrees(np.arccos(np.clip(abs(normal(f[ring_ix]) @ n_d), 0, 1))) for f in L])
    # ligand internal flexibility: heavy-atom RMSD after superposing ligand on ligand
    ligint = np.empty(N)
    for i, f in enumerate(L):
        R, pc, qc = kabsch_apply(f, d_lig_xyz)
        ligint[i] = np.sqrt((((f - pc) @ R + qc - d_lig_xyz) ** 2).sum(-1).mean())

    # displacement of the catalytic atoms from the designed positions (core-aligned): is a bad base
    # distance due to the ligand moving, or to the Asp / Gln side chains moving?
    def atom_xyz(struct, rid, names, protein=True):
        m = (struct.res_id == rid) & np.isin(struct.atom_name, names) & ((~struct.hetero) if protein else True)
        return struct.coord[m].astype(np.float64).mean(0)
    dsp = {"shift_N3": ("lig", "N3"), "shift_N1": ("lig", "N1"), "shift_baseO": ("base", ("OD1", "OD2")),
           "shift_donorNE2": ("donor", ("NE2",))}
    shifts = {}
    for k, (kind, nm) in dsp.items():
        if kind == "lig":
            ref = d_lig_xyz[lig_names.index(nm)]; cur = L[:, lig_names.index(nm)]
        else:
            rid = base_res if kind == "base" else donor_res
            ref = atom_xyz(dstruct, rid, nm)
            ix = np.where(prot & (top.res_id == rid) & np.isin(top.atom_name, nm))[0]
            cur = Xa[:, ix].mean(1)
        shifts[k] = np.linalg.norm(cur - ref, axis=1)

    # ---- per-sample catalytic geometry / burial / SASA / contacts ----------------
    frame = top.copy()
    bcols = ["b_dist", "b_ang_O", "b_ang_NH", "b_oop", "b_dih", "d_dist", "d_ang_N", "d_ang_ax", "d_oop", "d_dih",
             "geom_base", "geom_donor", "geometry", "burial_CG", "lig_sasa_rel", "n_contact_res", "nitro_hb",
             "base_ligN3_min"]
    rows = []
    iso = top[lig].copy()
    prot_atoms = np.where(prot)[0]; lig_atoms = np.where(lig)[0]
    res_of_atom = top.res_id[prot_atoms]
    contact_counts = np.zeros(n_res); contact_bound = np.zeros(n_res)
    bound = lig_rmsd <= BOUND_RMSD                         # ligand in (or very near) the designed pose
    min_dist = np.empty(N)
    sse_counts = np.zeros((n_res, 3))
    donor_mask = np.isin(top.atom_name[prot_atoms], ["N", "NE2", "ND2", "NZ", "NE", "NH1", "NH2", "OG", "OG1", "OH",
                                                     "ND1", "NE1", "SG"])
    nitro_ix = [lig_names.index(n) for n in ("O21", "O11")]
    frame.chain_id = np.where(lig, "B", "A").astype("U4")
    for i in range(N):
        frame.coord = Xa[i].astype(np.float32)
        try:
            m = geometry(frame, chain="A", base_resid=base_res, donor_resid=donor_res, lig_chain="B", lig_resid=1)
            g = host_score.geom_scores(m)
            m.update({k: g[k] for k in ("geom_base", "geom_donor", "geometry")})
        except ValueError:
            m = {k: np.nan for k in bcols}
        b = burial_count(frame, chain="A", base_resid=base_res)
        m["burial_CG"] = b["burial_CG"]
        iso.coord = frame.coord[lig]
        s_c = sasa(frame, vdw_radii="ProtOr")[lig_atoms].sum()
        s_i = sasa(iso, vdw_radii="ProtOr").sum()
        m["lig_sasa_rel"] = float(s_c / s_i)
        P = Xa[i][prot_atoms]
        tree = cKDTree(P)
        dists, idx = tree.query(L[i], k=1)
        near = set()
        for li in range(len(lig_atoms)):
            near.update(res_of_atom[j] for j in tree.query_ball_point(L[i][li], CONTACT_CUTOFF))
        min_dist[i] = float(dists.min())
        for r in near:
            contact_counts[int(r) - 1] += 1
            if bound[i]:
                contact_bound[int(r) - 1] += 1
        m["n_contact_res"] = len(near)
        hb = 0
        for oi in nitro_ix:
            dd = np.linalg.norm(P[donor_mask] - L[i][oi], axis=1)
            hb += int((dd < 3.3).any())
        m["nitro_hb"] = hb
        m["base_ligN3_min"] = float(np.linalg.norm(
            Xa[i][prot_atoms][(res_of_atom == base_res) & np.isin(top.atom_name[prot_atoms], ("OD1", "OD2"))]
            - L[i][lig_names.index("N3")], axis=1).min())
        sse = annotate_sse(frame[(~frame.hetero) & (frame.chain_id == "A")])
        for k, c in enumerate(("a", "b", "c")):
            sse_counts[:, k] += (sse == c)
        rows.append(m)

    # ---- PCA of the core-aligned CA cloud --------------------------------------
    flat = (CA - mean_ca).reshape(N, -1)
    U, S, Vt = np.linalg.svd(flat, full_matrices=False)
    var = S ** 2 / (S ** 2).sum()
    pcs = U[:, :3] * S[:3]

    # ---- write ----------------------------------------------------------------
    cols = ["frame", "ca_rmsd_design_all", "ca_rmsd_design_core", "ca_rmsd_to_mean", "rg_ca", "lig_pose_rmsd",
            "lig_com_shift", "lig_tilt_deg", "lig_internal_rmsd", "lig_min_dist", "bound", "shift_N3", "shift_N1", "shift_baseO", "shift_donorNE2"] + bcols + ["pc1", "pc2", "pc3"]
    with open(out / "frames.csv", "w") as fh:
        w = csv.writer(fh); w.writerow(cols)
        for i, m in enumerate(rows):
            w.writerow([i, ca_rmsd_design_all[i], ca_rmsd_design_core[i], ca_rmsd_mean[i], rg[i], lig_rmsd[i],
                        com_shift[i], tilt[i], ligint[i], min_dist[i], int(bound[i])] + [shifts[k][i] for k in dsp] + [m.get(k, np.nan) for k in bcols] + list(pcs[i]))
    with open(out / "residues.csv", "w") as fh:
        w = csv.writer(fh); w.writerow(["res", "aa", "rmsf", "rmsd_to_design", "core", "ligand_contact_frac",
                                        "contact_frac_bound", "helix_frac", "strand_frac", "role"])
        for k in range(n_res):
            role = "base" if k + 1 == base_res else "donor" if k + 1 == donor_res else ""
            w.writerow([k + 1, seq[k], f"{rmsf[k]:.3f}", f"{rmsd_design_res[k]:.3f}", int(core[k]),
                        f"{contact_counts[k] / N:.3f}", f"{contact_bound[k] / max(bound.sum(), 1):.3f}", f"{sse_counts[k, 0] / N:.3f}", f"{sse_counts[k, 1] / N:.3f}", role])
    # representative structures
    medoid = int(np.argmin(ca_rmsd_mean))
    pick = np.linspace(0, N - 1, 20).astype(int)
    f = PDBFile(); f.set_structure(traj[pick]); f.write(out / "ensemble_20frames.pdb")
    f = PDBFile(); f.set_structure(traj[medoid:medoid + 1]); f.write(out / "medoid.pdb")
    # (these are in the xtc's own frame: re-write them in the aligned frame)
    for name, ids in (("ensemble_20frames.pdb", pick), ("medoid.pdb", [medoid])):
        t = traj[list(ids)].copy(); t.coord = Xa[list(ids)].astype(np.float32)
        f = PDBFile(); f.set_structure(t); f.write(out / name)

    fr = {c: np.array([r.get(c, np.nan) for r in rows], float) for c in bcols}
    q = lambda v, p: float(np.nanpercentile(v, p))  # noqa: E731
    summ = dict(
        run=run, n_samples=N, n_residues=n_res, base=f"Asp{base_res}", donor=f"Gln{donor_res}",
        n_core=int(core.sum()),
        rmsf_median=float(np.median(rmsf)), rmsf_max=float(rmsf.max()), rmsf_max_res=int(rmsf.argmax() + 1),
        rmsf_motif=dict(base=float(rmsf[base_res - 1]), donor=float(rmsf[donor_res - 1])),
        ca_rmsd_design_core=dict(median=q(ca_rmsd_design_core, 50), p90=q(ca_rmsd_design_core, 90)),
        ca_rmsd_design_all=dict(median=q(ca_rmsd_design_all, 50), p90=q(ca_rmsd_design_all, 90)),
        rg=dict(mean=float(rg.mean()), sd=float(rg.std())),
        ligand=dict(pose_rmsd_median=q(lig_rmsd, 50), pose_rmsd_p90=q(lig_rmsd, 90),
                    frac_pose_within_1A=float((lig_rmsd < 1.0).mean()), frac_pose_within_2A=float((lig_rmsd < 2.0).mean()),
                    com_shift_median=q(com_shift, 50), tilt_median=q(tilt, 50), internal_rmsd_median=q(ligint, 50),
                    rmsf_by_atom={n: float(v) for n, v in zip(lig_names, lig_rmsf)}),
        geometry=dict(b_dist_median=q(fr["b_dist"], 50), b_dist_p10_p90=[q(fr["b_dist"], 10), q(fr["b_dist"], 90)],
                      d_dist_median=q(fr["d_dist"], 50), d_dist_p10_p90=[q(fr["d_dist"], 10), q(fr["d_dist"], 90)],
                      geometry_score_median=q(fr["geometry"], 50), geometry_score_p10_p90=[q(fr["geometry"], 10), q(fr["geometry"], 90)],
                      frac_b_dist_le_3p5=float(np.nanmean(fr["b_dist"] <= 3.5)),
                      frac_d_dist_le_3p5=float(np.nanmean(fr["d_dist"] <= 3.5)),
                      frac_both_le_3p5=float(np.nanmean((fr["b_dist"] <= 3.5) & (fr["d_dist"] <= 3.5))),
                      frac_geometry_ge_0p6=float(np.nanmean(fr["geometry"] >= 0.6))),
        burial_CG=dict(median=q(fr["burial_CG"], 50), p10_p90=[q(fr["burial_CG"], 10), q(fr["burial_CG"], 90)]),
        lig_sasa_rel=dict(median=q(fr["lig_sasa_rel"], 50), p90=q(fr["lig_sasa_rel"], 90),
                          frac_le_0p15=float(np.nanmean(fr["lig_sasa_rel"] <= 0.15))),
        nitro_hb_mean=float(np.nanmean(fr["nitro_hb"])),
        occupancy=dict(bound_rmsd_cut=BOUND_RMSD, frac_bound=float(bound.mean()),
                       frac_pose_within_5A=float((lig_rmsd < 5).mean()),
                       frac_ligand_no_contact=float((min_dist > CONTACT_CUTOFF).mean()),
                       lig_min_dist_median=q(min_dist, 50),
                       com_shift_pctiles_0_10_25_50_75_90_100=[q(com_shift, p) for p in (0, 10, 25, 50, 75, 90, 100)]),
        bound_only=(dict(
            n=int(bound.sum()),
            b_dist_median=q(fr["b_dist"][bound], 50), b_dist_p10_p90=[q(fr["b_dist"][bound], 10), q(fr["b_dist"][bound], 90)],
            d_dist_median=q(fr["d_dist"][bound], 50), d_dist_p10_p90=[q(fr["d_dist"][bound], 10), q(fr["d_dist"][bound], 90)],
            geometry_median=q(fr["geometry"][bound], 50), geometry_p10_p90=[q(fr["geometry"][bound], 10), q(fr["geometry"][bound], 90)],
            geometry_max=float(np.nanmax(fr["geometry"][bound])),
            frac_b_dist_le_3p5=float(np.nanmean(fr["b_dist"][bound] <= 3.5)),
            frac_d_dist_le_3p5=float(np.nanmean(fr["d_dist"][bound] <= 3.5)),
            frac_both_le_3p5=float(np.nanmean((fr["b_dist"][bound] <= 3.5) & (fr["d_dist"][bound] <= 3.5))),
            lig_sasa_rel_median=q(fr["lig_sasa_rel"][bound], 50), burial_CG_median=q(fr["burial_CG"][bound], 50),
            nitro_hb_mean=float(np.nanmean(fr["nitro_hb"][bound])), pose_rmsd_median=q(lig_rmsd[bound], 50),
            ligand_rmsf_by_atom={n: float(v) for n, v in zip(lig_names, np.sqrt(((L[bound] - L[bound].mean(0)) ** 2).sum(-1).mean(0)))},
            shift_from_design_median={k: q(v[bound], 50) for k, v in shifts.items()},
            tilt_median=q(tilt[bound], 50), internal_rmsd_median=q(ligint[bound], 50),
            contacts=sorted(((int(r), seq[r - 1], float(contact_bound[r - 1] / bound.sum())) for r in range(1, n_res + 1)),
                            key=lambda t: -t[2])[:12]) if bound.sum() >= 5 else None),
        helix_frac_overall=float(sse_counts[:, 0].sum() / (N * n_res)),
        strand_frac_overall=float(sse_counts[:, 1].sum() / (N * n_res)),
        pca_var_explained=[float(v) for v in var[:5]],
        pc1_corr=dict({c: float(np.corrcoef(pcs[:, 0], fr[c])[0, 1]) for c in ("b_dist", "d_dist", "geometry", "burial_CG")}
                      | {"lig_pose_rmsd": float(np.corrcoef(pcs[:, 0], lig_rmsd)[0, 1])}),
        corr=dict(b_dist_vs_d_dist=float(np.corrcoef(fr["b_dist"], fr["d_dist"])[0, 1]),
                  lig_rmsd_vs_geometry=float(np.corrcoef(lig_rmsd, fr["geometry"])[0, 1])),
        contacts=sorted(((int(r), seq[r - 1], float(contact_counts[r - 1] / N)) for r in range(1, n_res + 1)),
                        key=lambda t: -t[2])[:12],
        medoid_frame=medoid,
    )
    (out / "summary.json").write_text(json.dumps(summ, indent=2))
    return summ


if __name__ == "__main__":
    for run in sys.argv[1:]:
        s = analyse(run)
        print(f"{run}: {s['n_samples']} samples, ligand pose RMSD median {s['ligand']['pose_rmsd_median']:.2f} A, "
              f"b_dist {s['geometry']['b_dist_median']:.2f}, d_dist {s['geometry']['d_dist_median']:.2f}")
