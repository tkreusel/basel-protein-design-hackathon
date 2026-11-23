"""Run several RFD4 arms in ONE process (cold import is the expensive part). Usage: run_arms.py <arm> [<arm> ...]
arms: unsteered | pocket | pocket_loops | tmol | tmolcheck. N = env N_REP (replicates), batch 4."""
import os, sys, time, traceback, glob, torch
os.environ.setdefault("CLUSTER", "local")
import rfd4_proteina as rfd4
def main():
    R = os.path.expanduser("~/HACKATHON_ETH/rfd4_runs"); N = os.environ.get("N_REP", "4")
    t0 = time.time(); model = rfd4.Proteina.from_pretrained("nvidia/RFD4-Proteina", revision="a06e2509c7033b010fc43614e9fd840303588502") \
        if hasattr(rfd4, "Proteina") else None
    print(f"[t] model load {time.time()-t0:.0f}s", flush=True)
    for arm in sys.argv[1:]:
        t = time.time()
        try:
            if arm == "tmolcheck":
                import biotite.structure.io.pdbx as pdbx
                from rfproteina.metrics.physics.tmol_block_pair_metric import TMolBlockPairMetric
                m = TMolBlockPairMetric(region_a="is_ligand", region_b="is_protein", ligand_params_files=[R+"/inputs/6NT.tmol"])
                for f in sorted(glob.glob(f"{R}/out/20261121_smoke/inference_samples/*/*.cif")):
                    from atomworks.io.parser import parse; a = parse(f)["assemblies"]["1"][0]
                    t1 = time.time(); print(f.split("/")[-1], m.compute(a), f"{time.time()-t1:.1f}s", flush=True)
                continue
            if arm == "evaltmol":
                import pandas as pd
                from atomworks.io.parser import parse
                from rfproteina.metrics.custom_metrics.kemp_tmol import KempTMolMetric
                m = KempTMolMetric(region_a="is_ligand", region_b="is_protein", ligand_params_files=[R+"/inputs/6NT.tmol"]); rows = []
                for f in sorted(glob.glob(f"{R}/out/20261121_{os.environ.get('EVAL_GLOB', '*')}/inference_samples/*/*.cif")):
                    a = parse(f)["assemblies"]["1"][0]
                    try: v = m.compute(a)["interface_ddg"]
                    except Exception as e: v = float("nan"); print("fail", f, repr(e)[:150])
                    rows.append(dict(arm=f.split("/out/20261121_")[1].split("/")[0], design=os.path.basename(f), tmol_ddg=v))
                pd.DataFrame(rows).to_csv(f"{R}/out/20261121_eval_rfd4_tmol{os.environ.get('EVAL_TAG', '')}.csv", index=False); print(pd.DataFrame(rows).groupby("arm").tmol_ddg.describe().round(2).to_string())
                continue
            if arm.startswith("k3u_"):
                ov = [f"inputs={R}/inputs/kemp3_rfd4_{arm[4:]}.yaml", f"out_dir={R}/out/20261121_{arm}", f"run_name_={arm}", "~metrics"]
            elif arm.startswith("k3_"):
                ov = [f"inference_experiment=kemp3_{arm[3:]}"]
            elif arm.startswith("k2u_"):
                ov = [f"inputs={R}/inputs/kemp2_rfd4_{arm[4:]}.yaml", f"out_dir={R}/out/20261121_{arm}", f"run_name_={arm}", "~metrics"]
            elif arm.startswith("k2_"):
                ov = [f"inference_experiment=kemp2_{arm[3:]}"]
            elif arm == "unsteered":
                ov = [f"inputs={R}/inputs/kemp_rfd4.yaml", f"out_dir={R}/out/20261121_unsteered", "run_name_=unsteered", "~metrics"]
            else:
                ov = [f"inference_experiment=kemp_{arm}"]
            ov += [f"num_replicates={int(N)*2 if arm.startswith(('k2u_', 'k3u_')) else int(N)}", "batch_size=4", "compile_lengths=null", f"seed=7"]
            rfd4.generate(model, overrides=ov)
        except Exception:
            traceback.print_exc()
        print(f"[t] arm {arm} {time.time()-t:.0f}s", flush=True)

if __name__ == "__main__":
    main()
