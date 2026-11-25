"""Run several FK-steered RFD3 arms in ONE process (cold GPFS import is the expensive part).
usage: run_arms_rfd3.py <key> <arm> [<arm>...]   arm: unsteered | pocket | pocket_loops ; env N_PART (particles), N_TRAJ"""
import os, sys, time, traceback
sys.path.insert(0, "/gpfs/bwfor/work/ws/hd_la299-MPIP/rfdiff3_structure_guidance/src")
from guidance.config import load_inference_config
from guidance.run_inference import run_inference
H = os.path.expanduser("~/HACKATHON_ETH/rfd3"); key = sys.argv[1]
NP, NT = os.environ.get("N_PART", "12"), os.environ.get("N_TRAJ", "1")
for arm in sys.argv[2:]:
    t = time.time()
    try:
        cfgf = f"{H}/configs/kemp_{'pocket' if arm == 'unsteered' else arm}.yaml"
        ov = [f"example_key={key}_150_260", f"design_inputs={H}/inputs/kemp_rfd3_{key}_150_260.json", f"out_dir={H}/out/20261121_{key}_{arm}", f"rfd3.diffusion_batch_size={NP}", f"num_trajectories={NT}"]
        if arm == "unsteered": ov.append("resampling.potential_mode=blind")
        run_inference(load_inference_config(cfgf, ov))
    except Exception:
        traceback.print_exc()
    print(f"[t] arm {key}/{arm} {time.time()-t:.0f}s", flush=True)
