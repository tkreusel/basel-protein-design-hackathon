"""PLACER ensembles for every TSA replicate: holo (ligand predicted) and apo (no ligand).

RUNS IN THE placer_env (python 3.10, torch 2.3, dgl) - no biotite/rdkit there, so this file only uses numpy/pandas and
PLACER itself. All analysis is in analyze_placer.py (RFD4/tmol env).

  python -m maxi_eval.run_placer --inputs work/inputs --out work/placer -n 50 [--entries A,B] [--modes holo,apo]
Resumable: a replicate/mode with an existing CSV is skipped.
"""
from __future__ import annotations
import argparse, json, os, sys, time, traceback, warnings
from pathlib import Path
warnings.filterwarnings("ignore")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", default="work/inputs")
    ap.add_argument("--out", default="work/placer")
    ap.add_argument("-n", "--nsamples", type=int, default=50)
    ap.add_argument("--entries", default=None)
    ap.add_argument("--max-replicates", type=int, default=None)
    ap.add_argument("--modes", default="holo,apo")
    ap.add_argument("--arms", default=None, help="arms to run, separated by + or , (tsa,sub,prod,apo); default all. NOTE: with sbatch --export use + (a comma splits the export list)")
    ap.add_argument("--placer-root", default=os.environ.get("PLACER_ROOT", "/mnt/work/repos/placer"))
    ap.add_argument("--threads", type=int, default=None)
    a = ap.parse_args()
    sys.path.insert(0, a.placer_root)
    import torch
    if a.threads:
        torch.set_num_threads(a.threads)
    import PLACER
    placer = PLACER.PLACER()
    print("PLACER loaded; cuda:", torch.cuda.is_available(), flush=True)
    man = json.loads((Path(a.inputs) / "manifest.json").read_text())
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    for entry, e in man["entries"].items():
        if a.entries and entry not in a.entries.split(","):
            continue
        # every arm (6NT / H5J / 6VP) gets a holo ensemble; the apo ensemble does not depend on the ligand, so it is run once
        # per entry from the TSA replicate (apo needs the catalytic base, taken from that replicate's detected site)
        reps = list(e["replicates"])
        if a.arms:
            reps = [r for r in reps if r["arm"] in a.arms.replace("+", ",").split(",")]
        if a.max_replicates:
            seen: dict = {}
            keep = []
            for r in reps:
                key = (r["arm"], r["source"])
                seen[key] = seen.get(key, 0) + 1
                if seen[key] <= a.max_replicates:
                    keep.append(r)
            reps = keep
        for r in reps:
            code = {"tsa": "6NT", "sub": "H5J", "prod": "6VP"}.get(r["arm"])     # None for the ligand-free apo arm
            for mode in a.modes.split(","):
                if mode == "apo" and r["arm"] not in ("tsa", "apo"):
                    continue
                if mode == "holo" and r["arm"] == "apo":
                    continue
                label = r["rep_id"].replace(":", "_") + f"_{mode}"
                prefix = out / entry / label
                prefix.parent.mkdir(exist_ok=True)
                if Path(f"{prefix}.csv").exists():
                    print("skip (done)", label); continue
                t = time.time()
                try:
                    pin = PLACER.PLACERinput()
                    pin.pdb(r["pdb"])
                    pin.name(label)
                    if code:
                        pin.ligand_reference({code: "CCD"})
                        pin.ignore_ligand_hydrogens(True)
                    if mode == "holo":
                        pin.predict_ligand([("B", code, 1)])
                    else:  # apo: drop ligand, crop around the catalytic base
                        base = (r.get("site") or {}).get("base")
                        if not base:     # apo-arm replicate has no ligand/site: take the catalytic base from the entry's TSA replicate
                            base = next(((x.get("site") or {}).get("base") for x in e["replicates"] if x["arm"] == "tsa" and (x.get("site") or {}).get("base")), None)
                        if not base:
                            raise RuntimeError("apo mode needs a detected catalytic base (no carboxylate near the ligand)")
                        pin.exclude_sm(True)
                        pin.target_res(("A", int(base["resid"])))
                    outputs = placer.run(pin, a.nsamples)
                    PLACER.protocol.dump_output(outputs, str(prefix), rerank=None)
                    print(f"done {label} {time.time()-t:.0f}s", flush=True)
                except Exception as exc:  # fail loud, continue
                    print(f"FAILED {label}: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
                    Path(f"{prefix}.FAILED").write_text(traceback.format_exc())

if __name__ == "__main__":
    main()
