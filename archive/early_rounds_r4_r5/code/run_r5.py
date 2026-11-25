"""Round-5 driver (one process per job). usage: run_r5.py <campaign> <seed> <arm> [<arm> ...]  with arm in {unsteered, steer_base, steer_base_tmol, steer_base_channel}
Outputs: runs/<campaign>/structures/<arm>/s<seed>/ ; designs/ flat copies named <motif>__<arm>__s<seed>_<rep>_beam<k>.cif ; metrics/<arm>_s<seed>.csv
Per arm: generate -> flatten to structures/<arm>/designs/<motif>__<arm>__<lineage>.cif -> metrics/<arm>.csv (all metrics incl. TMol).
Env: N_REP (replicates per steered arm; unsteered gets 2x)."""
import os, sys, glob, shutil, time, traceback
os.environ.setdefault("CLUSTER", "local")
import numpy as np, pandas as pd
import rfd4_proteina as rfd4
H = os.path.expanduser("~/HACKATHON_ETH")
MK = {"20261121_r4_hg317-tips_with-trp": ("hgwt", "kemp4_hg317-tips_with-trp.yaml", "hg"), "20261121_r4_hg317-tips_no-trp": ("hgnt", "kemp4_hg317-tips_no-trp.yaml", "hg"),
      "20261121_r4_gdm1483-glu-ser": ("gdm", "kemp4_gdm1483-glu-ser.yaml", "gdm")}
def evaluate(camp, arm, seed):
    import biotite.structure as s, biotite.structure.io.pdbx as pdbx
    from atomworks.io.parser import parse
    from rfproteina.metrics.custom_metrics.kemp_base import KempBaseMetric
    from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
    from rfproteina.metrics.custom_metrics.kemp_chem import KempChemMetric
    from rfproteina.metrics.custom_metrics.kemp_tmol import KempTMolMetric
    from rfproteina.metrics.custom_metrics.kemp_channel import KempChannelMetric
    mk, spec, kind = MK[camp]
    ch = KempChemMetric(ref_pdb=f"{H}/inputs/HG3_17_kemp_motif.pdb", base_res="ASP", base_match=("CB", "CG", "OD1", "OD2"), base_o=("OD1", "OD2"), hb_res="GLN",
                        hb_match=("CG", "CD", "OE1", "NE2"), hb_atom="NE2") if kind == "hg" else KempChemMetric(ref_pdb=f"{H}/inputs/kemp_theozyme_GDM_KE_1483.pdb")
    mb, mp, mch = KempBaseMetric(), KempPocketMetric(target_ratio=0.5, tol=0.10), KempChannelMetric()
    tm = KempTMolMetric(region_a="is_ligand", region_b="is_protein", ligand_params_files=[f"{H}/inputs/6NT.tmol"])
    d = f"{H}/runs/{camp}/structures/{arm}/s{seed}"; os.makedirs(f"{d}/designs", exist_ok=True); rows = []
    for f in sorted(glob.glob(f"{d}/inference_samples/*/*_generated.cif")):
        stem = os.path.basename(os.path.dirname(f)); name = f"{mk}__{arm}__s{seed}_{stem.replace('kemp4_', '').replace('_bsc', '_beam')}"
        shutil.copy(f, f"{d}/designs/{name}.cif")
        a = pdbx.get_structure(pdbx.CIFFile.read(f), model=1, use_author_fields=False); a = a[(a.element != "H") & ~np.isnan(a.coord).any(1)]
        c = ch.compute(a); r = {**mp.compute(a), **mb.compute(a), **mch.compute(a), "motif_base_rmsd": c["base_rmsd"], "motif_partner_rmsd": c["hb_rmsd"],
                                "serOG_or_Gln_to_ligand_dist": c["hb_partner_dist"]}
        sse = np.array(s.annotate_sse(a[~a.hetero])); r.update(campaign=camp, arm=arm, seed=seed, design=name, length=len(sse), helix=float((sse == "a").mean()), strand=float((sse == "b").mean()))
        try: r["tmol_ddg"] = tm.compute(parse(f)["assemblies"]["1"][0])["interface_ddg"]
        except Exception as e: r["tmol_ddg"] = float("nan"); print("tmol fail", name, repr(e)[:120], flush=True)
        rows.append(r)
    os.makedirs(f"{H}/runs/{camp}/metrics", exist_ok=True); pd.DataFrame(rows).to_csv(f"{H}/runs/{camp}/metrics/{arm}_s{seed}.csv", index=False)
    print(f"[t] metrics {camp}/{arm}/s{seed}: {len(rows)} designs", flush=True)
def main():
    camp, seed, arms = sys.argv[1], int(sys.argv[2]), sys.argv[3:]; mk, spec, kind = MK[camp]; N = int(os.environ.get("N_REP", "3"))
    model = rfd4.Proteina.from_pretrained("nvidia/RFD4-Proteina", revision="a06e2509c7033b010fc43614e9fd840303588502")
    for arm in arms:
        t = time.time()
        try:
            od = f"{H}/runs/{camp}/structures/{arm}/s{seed}"
            if arm == "unsteered":
                ov = [f"inputs={H}/inputs/{spec}", f"out_dir={od}", f"run_name_={mk}_unsteered_s{seed}", "~metrics", f"num_replicates={2*N}"]
            else:
                exp = f"r5_{mk}_{arm}" if arm == "steer_base_channel" else f"r4_{mk}_{arm}"
                ov = [f"inference_experiment={exp}", f"out_dir={od}", f"run_name_={mk}_{arm}_s{seed}", f"num_replicates={N}"]
            rfd4.generate(model, overrides=ov + ["batch_size=4", "compile_lengths=null", f"seed={seed}"])
            print(f"[t] generated {camp}/{arm} {time.time()-t:.0f}s", flush=True); evaluate(camp, arm, seed)
        except Exception:
            traceback.print_exc()
        print(f"[t] arm {camp}/{arm} total {time.time()-t:.0f}s", flush=True)
if __name__ == "__main__":
    main()
