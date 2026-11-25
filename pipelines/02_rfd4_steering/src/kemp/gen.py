"""Stage 1 (GPU): RFD4 generation and per-design metrics for one batch.

usage: gen.py --motif hgnt|hgwt --arm <arm> --batch N --out DIR [--n-rep R] [--seed S] [--chunk 25] [--batch-size 4]

A batch is 300 designs: steered arms 100 replicates x 3 beam members, unsteered 300 replicates x 1 (config/arms.yaml).
Generation runs in chunks of `--chunk` replicates, each a separate rfd4.generate call, so a killed job resumes at the
first unfinished chunk. Output:
    DIR/designs/<motif>__<arm>__bNN_cXXXX_<rfd4 name>.cif    flat copy of every generated design
    DIR/metrics.csv                                           one row per design (geometry, pocket, channel, TMol dG, SSE)
"""
import argparse
import glob
import os
import shutil
import time
import traceback

os.environ.setdefault("CLUSTER", "local")
import numpy as np
import pandas as pd

from . import config

CHUNK_SEED_STRIDE = 100_000


def seed_for(batch, arm, motif, chunk_index=0, arms_cfg=None):
    arms_cfg = arms_cfg or config.load("arms")
    base = 1000 * batch + 10 * arms_cfg["arms"][arm]["seed_offset"] + arms_cfg["motif_seed_offset"][motif]
    return base + CHUNK_SEED_STRIDE * chunk_index


def replicates_for(arm, arms_cfg=None):
    return (arms_cfg or config.load("arms"))["arms"][arm]["replicates"]


def hydra_overrides(motif, arm, batch, out_dir, n_rep, seed, batch_size):
    """Overrides for rfd4.generate. Unsteered arms drop the metrics group (conditioning only)."""
    base = config.arm_family(arm)
    spec = config.INPUTS / config.motif_spec(motif, sheet=arm.endswith("_sheet"))
    common = [f"inputs={spec}", f"out_dir={out_dir}", f"run_name_={motif}_{arm}_b{batch}", f"num_replicates={n_rep}"]
    if config.is_steered(arm):
        common.insert(0, f"inference_experiment=kc_{motif}_{base}")
    else:
        common.insert(3, "~metrics")
    return common + [f"batch_size={batch_size}", "compile_lengths=null", f"seed={seed}"]


def build_metrics():
    """The metric objects used for the per-design table (same settings as config/arms.yaml's reward metrics)."""
    from rfproteina.metrics.custom_metrics.kemp_base import KempBaseMetric
    from rfproteina.metrics.custom_metrics.kemp_pocket import KempPocketMetric
    from rfproteina.metrics.custom_metrics.kemp_chem import KempChemMetric
    from rfproteina.metrics.custom_metrics.kemp_tmol import KempTMolMetric
    from rfproteina.metrics.custom_metrics.kemp_channel import KempChannelMetric
    motifs = config.load("motifs")
    tips = {k: v["tips"] for k, v in motifs["residues"].items()}
    chem = KempChemMetric(ref_pdb=str(config.PIPELINE / motifs["reference_pdb"]), base_res="ASP", base_match=tuple(tips["ASP"]),
                          base_o=tuple(motifs["base_oxygens"]), hb_res="GLN", hb_match=tuple(tips["GLN"]), hb_atom=motifs["donor_atom"])
    band = config.load("gates")["design_gate"]["channel_band"]
    radii = config.load("arms")["channel_band_radii"]
    return dict(chem=chem, base=KempBaseMetric(), pocket=KempPocketMetric(target_ratio=0.5, tol=0.10),
                channel=KempChannelMetric(lo=band[0], hi=band[1], radii=tuple(radii)),
                tmol=KempTMolMetric(region_a="is_ligand", region_b="is_protein", ligand_params_files=[str(config.INPUTS / "6NT.tmol")]))


def design_metrics(cif, motif, arm, name, m):
    import biotite.structure as bs
    import biotite.structure.io.pdbx as pdbx
    from atomworks.io.parser import parse
    atoms = pdbx.get_structure(pdbx.CIFFile.read(cif), model=1, use_author_fields=False)
    atoms = atoms[(atoms.element != "H") & ~np.isnan(atoms.coord).any(axis=1)]
    chem = m["chem"].compute(atoms)
    row = {**m["pocket"].compute(atoms), **m["base"].compute(atoms), **m["channel"].compute(atoms),
           "motif_base_rmsd": chem["base_rmsd"], "motif_partner_rmsd": chem["hb_rmsd"], "serOG_or_Gln_to_ligand_dist": chem["hb_partner_dist"]}
    sse = np.array(bs.annotate_sse(atoms[~atoms.hetero]))
    row.update(motif=motif, arm=arm, design=name, length=len(sse), helix=float((sse == "a").mean()), strand=float((sse == "b").mean()))
    try:
        row["tmol_ddg"] = m["tmol"].compute(parse(cif)["assemblies"]["1"][0])["interface_ddg"]
    except Exception as exc:                                   # a few chain-gap structures fail in TMol; the design is kept without a value
        row["tmol_ddg"] = float("nan")
        print("tmol fail", name, repr(exc)[:120], flush=True)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--motif", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--batch", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-rep", type=int)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--chunk", type=int, default=25, help="replicates per generate() call (progress and resume granularity)")
    ap.add_argument("--batch-size", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true", help="print the overrides of the first chunk and exit")
    a = ap.parse_args()

    n_rep = a.n_rep or replicates_for(a.arm)
    out = os.path.abspath(a.out)
    if a.dry_run:
        seed = a.seed if a.seed is not None else seed_for(a.batch, a.arm, a.motif)
        print("\n".join(hydra_overrides(a.motif, a.arm, a.batch, f"{out}/chunks/c0000", min(a.chunk, n_rep), seed, a.batch_size)))
        return

    os.environ.setdefault("PROJECT_ROOT", str(config.PIPELINE / "rfd4"))
    os.environ.setdefault("DATA_PATH", str(config.PIPELINE / "rfd4" / "data_path"))
    os.environ.setdefault("HF_HUB_OFFLINE", "1")
    import rfd4_proteina as rfd4
    model_cfg = config.load("arms")["model"]
    metrics = build_metrics()
    os.makedirs(out + "/designs", exist_ok=True)
    os.makedirs(out + "/chunks", exist_ok=True)
    model = rfd4.Proteina.from_pretrained(model_cfg["weights"], revision=model_cfg["revision"])
    t0 = time.time()
    for ci, c0 in enumerate(range(0, n_rep, a.chunk)):
        n = min(a.chunk, n_rep - c0)
        chunk_dir, done = f"{out}/chunks/c{c0:04d}", f"{out}/chunks/c{c0:04d}/DONE"
        if os.path.exists(done):
            continue
        shutil.rmtree(chunk_dir, ignore_errors=True)
        seed = (a.seed if a.seed is not None else seed_for(a.batch, a.arm, a.motif)) + CHUNK_SEED_STRIDE * ci
        t = time.time()
        try:
            rfd4.generate(model, overrides=hydra_overrides(a.motif, a.arm, a.batch, chunk_dir, n, seed, a.batch_size))
            gen_s = time.time() - t
            rows = []
            for f in sorted(glob.glob(f"{chunk_dir}/inference_samples/*/*_generated.cif")):
                stem = os.path.basename(os.path.dirname(f)).replace("kemp4_", "").replace("_bsc", "_beam")
                name = f"{a.motif}__{a.arm}__b{a.batch:02d}_c{c0:04d}_{stem}"
                shutil.copy(f, f"{out}/designs/{name}.cif")
                try:
                    rows.append(design_metrics(f, a.motif, a.arm, name, metrics) | {"batch": a.batch, "chunk": c0, "seed": seed})
                except Exception:
                    print("metrics fail", name)
                    traceback.print_exc()
            pd.DataFrame(rows).to_csv(f"{chunk_dir}/metrics.csv", index=False)
            open(done, "w").write(f"{len(rows)} designs gen {gen_s:.0f}s total {time.time() - t:.0f}s\n")
            print(f"[t] chunk {c0}+{n}: {len(rows)} designs, generate {gen_s:.0f}s, metrics {time.time() - t - gen_s:.0f}s", flush=True)
        except Exception:
            traceback.print_exc()
    parts = [pd.read_csv(f) for f in sorted(glob.glob(f"{out}/chunks/c*/metrics.csv"))]
    if parts:
        pd.concat(parts, ignore_index=True).to_csv(out + "/metrics.csv", index=False)
    print(f"[t] batch done {a.motif}/{a.arm}/b{a.batch}: {sum(len(x) for x in parts)} designs in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
