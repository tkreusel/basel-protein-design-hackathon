"""Readiness + cost check for a Boltz output folder (no GPU, no jobs, reads files only).

  python -m maxi_eval.check_boltz_ready --root /mnt/work/projects/Maxi/boltz_outputs [--require-pae]

Reports per protein which of the three arms (6NT tsa / H5J sub / 6VP prod) have a finished model (+ PAE/confidence for ipSAE2), lists
incomplete proteins, and estimates run time for the full suite from measured speeds. Exit code 0 = every protein has all three arms.
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path
from collections import defaultdict

RX = re.compile(r"^(?P<entry>.+?)(?:__|_)(?P<arm>tsa|sub|prod|6NT|H5J|6VP|analogon|analogue|analog|substrate|product|apo)$")
ARM = {"tsa": "tsa", "sub": "sub", "prod": "prod", "6NT": "tsa", "H5J": "sub", "6VP": "prod", "analogon": "tsa", "analogue": "tsa", "analog": "tsa",
       "substrate": "sub", "product": "prod", "apo": "apo"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--require-pae", action="store_true")
    ap.add_argument("--nsamples", type=int, default=20)
    a = ap.parse_args()
    from .prepare_inputs import infer_arm_from_cif
    found: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    unnamed = []
    for cif in sorted(Path(a.root).rglob("*_model_*.cif")):
        m = re.match(r"^(?P<name>.+)_model_(?P<k>\d+)\.cif$", cif.name)
        name = m["name"]; mm = RX.match(name)
        entry = (mm["entry"] if mm else name).lower()
        arm = infer_arm_from_cif(cif)                      # the ligand actually in the file decides, exactly as prepare_inputs does
        if arm is None and mm and ARM[mm["arm"]] == "apo":
            arm = "apo"
        if arm is None:
            unnamed.append(cif.name); continue
        if mm and ARM[mm["arm"]] != arm:
            print(f"  [warn] {cif.name}: name says {ARM[mm['arm']]} but the ligand is {arm}")
        ok_pae = (cif.parent / f"pae_{name}_model_{m['k']}.npz").exists() and (cif.parent / f"confidence_{name}_model_{m['k']}.json").exists()
        found[entry][arm].append(ok_pae)
    n = len(found)
    complete = [e for e, arms in found.items() if all(x in arms for x in ("tsa", "sub", "prod"))]
    print(f"{n} proteins found; {len(complete)} have all three arms (6NT, H5J, 6VP); {len(unnamed)} models with no 6NT/H5J/6VP ligand (ignored)")
    inc = {e: sorted(set(("tsa", "sub", "prod")) - set(arms)) for e, arms in found.items() if e not in complete}
    print(f"  apo predictions: {sum(1 for arms in found.values() if 'apo' in arms)} of {n} proteins (optional; enables apo-vs-holo preorganisation metrics)")
    for e, miss in list(inc.items())[:30]:
        print(f"  incomplete: {e}: missing {', '.join(miss)}")
    if a.require_pae:
        nopae = [e for e, arms in found.items() for arm, l in arms.items() if not all(l)]
        print(f"  models lacking PAE/confidence files: {len(nopae)} proteins affected")
    models = sum(len(l) for arms in found.values() for l in arms.values()) or len(unnamed)
    N = max(n, len({re.sub(r'_model_\d+\.cif$', '', c) for c in unnamed}), 1)
    placer_h = N * (3 + 1) * a.nsamples * 2.6 / 3600       # holo x3 arms + apo x1, 2.6 s/model measured on an L4
    tmol_h = N * (7 + 3) / 60                              # ~7 min TSA protocol + ~3 min for the two native arms, 8 cores
    ros_h = N * (6 + 3) / 60                               # ~6 min + native arms, ONE core each
    print(f"\nestimated time for {N} proteins ({models} models):")
    print(f"  PLACER (1 GPU)            ~{placer_h:.1f} h   [{N*4} ensembles x {a.nsamples} models x 2.6 s]")
    print(f"  tmol (CPU, 2 procs)       ~{tmol_h/2:.1f} h   [serial {tmol_h:.1f} h]")
    print(f"  PyRosetta (CPU, 3 procs)  ~{ros_h/3:.1f} h   [serial {ros_h:.1f} h]")
    print(f"  PLACER analysis (CPU)     ~{N*4*1.3/60:.1f} h; ipSAE2 + aggregate: minutes")
    sys.exit(0 if len(complete) == n and n > 0 else 1)

if __name__ == "__main__":
    main()
