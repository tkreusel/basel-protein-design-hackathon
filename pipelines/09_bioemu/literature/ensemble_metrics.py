"""Per-enzyme BioEmu-2 ensemble metrics for the literature Kemp eliminases (and our designs), for correlation
with kcat/KM. One process per enzyme; writes literature/metrics/<name>.json + frames_<name>.csv.

Catalytic positions in the 312-residue BioEmu sequences of the HG3 family (crystal 7K4Z numbering + 2):
base Asp129 (= Asp127), oxyanion donor residue 52 (= residue 50: Gln in most variants, Lys in HG3, His in HG3.3),
wall Trp46 (= Trp44). Our designs use their own indices (base 77, donor 50). 1A53-core is a different scaffold:
only site-independent metrics are meaningful there.

Samples are independent equilibrium draws. All frames are superposed on the CA atoms of an iteratively chosen
rigid core (residues with spread < 1.5 x median) so ligand motion is relative to the protein.
Usage: python -I ensemble_metrics.py <name> [--base N --donor N]
"""
import argparse, json, sys, warnings
from pathlib import Path
import numpy as np
import yaml
from biotite.structure import sasa
from biotite.structure.io.pdbx import CIFFile, get_structure
import biotite.structure.io.xtc as xtc
from scipy.spatial import cKDTree

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
from host_geometry import RING, angle, dihedral, _u, _plane_normal, burial_count, circ_mean, circ_sd  # noqa: E402
OUT = ROOT / "bioemu/new_zip/bioemu_outputs/outputs"
REF = yaml.safe_load(open(ROOT / "configs/host_reference.yaml"))["metrics"]
SC_SIDE = {"N", "O", "OXT", "C", "CA", "CB"}


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    v, _, wt = np.linalg.svd((P - pc).T @ (Q - qc))
    if np.linalg.det(v @ wt) < 0:
        v[:, -1] *= -1
    return v @ wt, pc, qc


def run(name, base, donor, outdir):
    d = OUT / name
    top = get_structure(CIFFile.read(d / "topology.cif"), model=1)
    X = xtc.XTCFile.read(str(d / "samples.xtc")).get_structure(template=top).coord.astype(np.float64)
    N = len(X)
    prot, lig = top.chain_id == "A", top.chain_id == "B"
    pi, li = np.where(prot)[0], np.where(lig)[0]
    ca_idx = np.where(prot & (top.atom_name == "CA"))[0]; n_res = len(ca_idx)
    lnames = list(top.atom_name[lig]); L = lambda a: lnames.index(a)  # noqa: E731
    # ---- alignment on an iterative core, against the ensemble mean ----------
    ref = X[0][ca_idx].copy(); core = np.ones(n_res, bool)
    for _ in range(4):
        Xa = np.empty_like(X)
        for i, f in enumerate(X):
            R, pc, qc = kabsch(f[ca_idx][core], ref[core]); Xa[i] = (f - pc) @ R + qc
        CA = Xa[:, ca_idx]; ref = CA.mean(0)
        dev = np.sqrt(((CA - ref) ** 2).sum(-1).mean(0)); core = dev < 1.5 * np.median(dev)
    rmsf = np.sqrt(((CA - ref) ** 2).sum(-1).mean(0))
    Lx = Xa[:, li]
    rid = top.res_id; an = top.atom_name
    sel = lambda r, names: np.where(prot & (rid == r) & np.isin(an, names))[0]  # noqa: E731
    bo = sel(base, ["OD1", "OD2"]); bcg = sel(base, ["CG"])
    dsc = np.where(prot & (rid == donor) & ~np.isin(an, list(SC_SIDE)) & np.isin(top.element, ["N", "O"]))[0]
    n3, n1 = L("N3"), L("N1")
    iso = top[lig].copy()
    frame = top.copy()
    s_iso = sasa(iso, vdw_radii="ProtOr").sum() if True else 1
    cols = {k: np.full(N, np.nan) for k in ["b_dist", "d_dist", "b_ang_O", "b_ang_NH", "b_oop", "b_dih", "sasa_rel",
                                              "n_contact", "burial_CG", "chi1", "chi2", "geom_base"]}
    ring_ix = [L(a) for a in RING]
    for i in range(N):
        P = Xa[i]; Lp = P[li]
        ob = P[bo][np.argmin(np.linalg.norm(P[bo] - Lp[n3], axis=1))]
        cg = P[bcg][0]
        c3a, n2, c4 = Lp[L("C3A")], Lp[L("N2")], Lp[L("C4")]
        v = ob - Lp[n3]; nh = _u(-(_u(c3a - Lp[n3]) + _u(n2 - Lp[n3])))
        nrm = _plane_normal(Lp[ring_ix])
        cols["b_dist"][i] = np.linalg.norm(v)
        cols["b_ang_O"][i] = angle(cg, ob, Lp[n3])
        cols["b_ang_NH"][i] = np.degrees(np.arccos(np.clip(np.dot(_u(v), nh), -1, 1)))
        cols["b_oop"][i] = np.degrees(np.arcsin(abs(np.dot(_u(v), nrm))))
        cols["b_dih"][i] = abs(dihedral(ob, Lp[n3], c3a, c4))
        if len(dsc):
            cols["d_dist"][i] = np.linalg.norm(P[dsc] - Lp[n1], axis=1).min()
        # burial of the base environment
        cols["burial_CG"][i] = (np.linalg.norm(P[ca_idx] - cg, axis=1) <= 12).sum() - 1
        frame.coord = P.astype(np.float32); iso.coord = frame.coord[li]
        cols["sasa_rel"][i] = sasa(frame, vdw_radii="ProtOr")[li].sum() / s_iso
        tree = cKDTree(P[pi])
        cols["n_contact"][i] = len({j for a in Lp for j in tree.query_ball_point(a, 4.0)})
        # base side-chain torsions (preorganisation of the nucleophile)
        n_, ca_, cb_, cg_ = (P[sel(base, [a])[0]] for a in ("N", "CA", "CB", "CG"))
        od1 = P[sel(base, ["OD1"])[0]]
        cols["chi1"][i] = dihedral(n_, ca_, cb_, cg_); cols["chi2"][i] = dihedral(ca_, cb_, cg_, od1)
    # host-style base geometry score (calibrated on HG3 crystals)
    gs = []
    for k in ("b_dist", "b_ang_O", "b_ang_NH", "b_oop", "b_dih"):
        r = REF[k]; z = (cols[k] - r["mean"]) / r["sigma"]; gs.append(np.exp(-0.5 * z * z))
    cols["geom_base"] = np.mean(gs, 0)
    sa, bd = cols["sasa_rel"], cols["b_dist"]
    pocket = (sa <= 0.25) & (bd <= 8.0)               # ligand buried AND near the base
    nac = pocket & (bd <= 3.5)                        # near-attack: base O within 3.5 A of N3
    # ligand pose spread inside the pocket (after protein alignment)
    spread = np.nan
    if pocket.sum() >= 10:
        mean_pose = Lx[pocket].mean(0)
        spread = float(np.median(np.sqrt(((Lx[pocket] - mean_pose) ** 2).sum(-1).mean(1))))
    # pocket CA rmsf: residues touching the ligand in >= 50 % of pocket frames
    pk_rmsf = np.nan
    if pocket.sum() >= 10:
        cnt = np.zeros(n_res)
        res_of = rid[pi]
        for i in np.where(pocket)[0][:200]:
            t = cKDTree(Xa[i][pi]); near = {int(res_of[j]) for a in Lx[i] for j in t.query_ball_point(a, 4.0)}
            for r in near: cnt[r - 1] += 1
        m = cnt >= 0.5 * min(200, pocket.sum()); pk_rmsf = float(rmsf[m].mean()) if m.any() else np.nan
    q = lambda v, p: float(np.nanpercentile(v, p)) if np.isfinite(v).any() else np.nan  # noqa: E731
    sub = lambda v, m: v[m] if m.sum() else np.array([np.nan])  # noqa: E731
    res = dict(
        name=name, n=N, n_res=n_res, base=base, donor=donor, donor_aa=str(top.res_name[sel(donor, ['CA'])[0]]),
        p_pocket=float(pocket.mean()), p_nac=float(nac.mean()), p_nac_given_pocket=float(nac.sum() / max(pocket.sum(), 1)),
        p_b35=float((bd <= 3.5).mean()), p_b40=float((bd <= 4.0).mean()), b_dist_p10=q(bd, 10), b_dist_med=q(bd, 50),
        b_dist_pocket_med=q(sub(bd, pocket), 50),
        d_dist_pocket_med=q(sub(cols["d_dist"], pocket), 50),
        p_d35_given_pocket=float(np.nanmean(sub(cols["d_dist"], pocket) <= 3.5)) if len(dsc) else np.nan,
        geom_base_mean=float(np.nanmean(cols["geom_base"])), geom_base_pocket=float(np.nanmean(sub(cols["geom_base"], pocket))),
        geom_base_p90=q(cols["geom_base"], 90),
        sasa_rel_med=q(sa, 50), contact_med=q(cols["n_contact"], 50), burial_CG_med=q(cols["burial_CG"], 50),
        pose_spread_pocket=spread, pocket_ca_rmsf=pk_rmsf, global_ca_rmsf_med=float(np.median(rmsf)),
        base_chi1_sd=circ_sd(cols["chi1"]), base_chi2_sd=circ_sd(cols["chi2"]),
        base_chi1_sd_pocket=circ_sd(sub(cols["chi1"], pocket)) if pocket.sum() >= 5 else np.nan,
        base_chi2_sd_pocket=circ_sd(sub(cols["chi2"], pocket)) if pocket.sum() >= 5 else np.nan,
        lig_com_sd=float(np.sqrt(((Lx.mean(1) - Lx.mean(1).mean(0)) ** 2).sum(1).mean())),
    )
    outdir.mkdir(exist_ok=True)
    (outdir / f"{name}.json").write_text(json.dumps(res, indent=1))
    np.savez_compressed(outdir / f"frames_{name}.npz", **cols, pocket=pocket, nac=nac)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("name")
    ap.add_argument("--base", type=int, default=129); ap.add_argument("--donor", type=int, default=52)
    a = ap.parse_args()
    r = run(a.name, a.base, a.donor, HERE / "metrics")
    print(a.name, f"p_pocket {r['p_pocket']:.2f} p_nac {r['p_nac']:.3f} b_dist_med {r['b_dist_med']:.1f} geom_base {r['geom_base_mean']:.2f}")
