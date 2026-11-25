"""Rebuild tests/fixtures from the campaign results (needs the run archives and the results tables; only for maintainers).

usage: python tests/make_fixtures.py ARCHIVE_DIR RESULTS_DIR

ARCHIVE_DIR  folder with <motif>__<arm>__bNN.tar.gz (kemp_results/run_archives)
RESULTS_DIR  kemp_results (tables/)
The fixtures are small excerpts of what the campaign produced, so the tests can check that the code still reproduces it.
"""
import shutil
import sys
import tarfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

FIX = Path(__file__).parent / "fixtures"
RUNS = ["hgnt__steer_tmol__b01", "hgwt__unsteered__b01", "hgnt__unsteered_sheet__b51"]


def trimmed_cif(src, dst):
    import biotite.structure.io.pdbx as pdbx
    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(src)), model=1, use_author_fields=False)
    atoms = atoms[(atoms.element != "H") & ~np.isnan(atoms.coord).any(axis=1)]
    f = pdbx.CIFFile()
    pdbx.set_structure(f, atoms, data_block="x")
    f.write(str(dst))


def main(archives, results):
    archives, results = Path(archives), Path(results)
    for name in RUNS + ["hgnt__unsteered__b01"]:
        tmp = Path(tempfile.mkdtemp())
        tarfile.open(archives / f"{name}.tar.gz").extractall(tmp, filter="data")
        run = tmp / name
        out = FIX / name
        out.mkdir(parents=True, exist_ok=True)
        if name in RUNS:
            for f in ("design_metrics.csv", "boltz1_eval.csv", "survivors.txt"):
                shutil.copy(run / f, out / f)
            gate = pd.read_csv(run / "gate.csv")
            gate[["design", "pass_geom", "pass_seq", "pass_core", "gate_tier", "lineage", "lineage_rank", "gate_pass", "n_ok", "channel_ok", "channel18_ok"]].to_csv(out / "expected_gate.csv", index=False)
        if name == RUNS[0]:
            for f in ("ranking.csv", "stage5_report.csv", "sequences.csv"):
                shutil.copy(run / f, out / f)
        if name in (RUNS[0], RUNS[1]):
            seqs = pd.read_csv(run / "sequences.csv")
            sel = pd.read_csv(run / "sel.csv")
            pick = sorted(sel.design, key=lambda d: (run / "designs" / f"{d}.cif").stat().st_size)[:2]
            seqs[seqs.design.isin(pick)].to_csv(out / "sequences_all.csv", index=False)
            (out / "designs").mkdir(exist_ok=True)
            for d in pick:
                trimmed_cif(run / "designs" / f"{d}.cif", out / "designs" / f"{d}.cif")
        if name == "hgnt__unsteered__b01":
            strict = pd.read_csv(results / "final_strict_all_gates_ipsae0.8" / "strict_all_gates_with_geometry.csv")
            strict = strict[(strict.motif == "hgnt") & (strict.arm == "unsteered") & (strict.batch == "b01")]
            (out / "designs").mkdir(exist_ok=True)
            (out / "best_boltz").mkdir(exist_ok=True)
            seqs = pd.read_csv(run / "sequences.csv")
            seqs[seqs.design.isin(strict.design)].to_csv(out / "sequences_all.csv", index=False)
            for x in strict.itertuples():
                trimmed_cif(run / "designs" / f"{x.design}.cif", out / "designs" / f"{x.design}.cif")
                trimmed_cif(run / "best_boltz" / f"{x.seq_id}__6NT_model_0.cif", out / "best_boltz" / f"{x.seq_id}__6NT_model_0.cif")
            cols = ["seq_id", "design", "motif", "design_tips_rmsd_theozyme", "boltz_best_tips_rmsd_theozyme", "asp_O_to_N3_A", "gln_NE2_to_N1_A",
                    "n_polar_atoms_3p5_of_N1", "base_CB_within_12A_boltz", "base_CB_within_12A_design"]
            strict[cols].to_csv(out / "expected_geometry.csv", index=False)
        shutil.rmtree(tmp)

    tables = results / "tables"
    full = pd.read_csv(tables / "evaluated_sequences_full_suite.csv", low_memory=False)            # all 528 evaluated sequences
    final = pd.read_csv(tables / "evaluated_with_final_metrics.csv", low_memory=False)            # non-sheet sequences with S07 / S12
    need = ["seq_id", "motif", "arm", "variant", "eval_set", "boltz_tier", "ipsae.tsa", "6NT.iptm.mean", "6NT.iptm.sd", "6NT.n", "H5J.iptm.mean", "H5J.iptm.sd", "H5J.n",
            "6VP.iptm.mean", "6VP.iptm.sd", "6VP.n", "placer.apo.frac_pass_core", "placer.holo.frac_pass_core", "ros.tsa.dG_sep"]
    s07 = ["S07_tmol_polar_residue_int_TSA_minus_SUB", "S12_boltz_cat_to_pocket_PAE_SUB_minus_PROD"]
    raw = ["tsa_pk.int.lk_ball_iso", "prod.int.lk_ball_iso", "sub.int.fa_elec", "prod.int.fa_elec", "tsa_pk.int.total", "tsa_pk.E_total", "prod.E_total",
           "tsa_pk.buns.n_unsat_ligand", "prod.buns.n_unsat_ligand"]
    panel = ["boltz_ipde.tsa", "6NT.pocket_rmsd.mean", "ros.delta_native.tsa_minus_sub.dG_sep", "delta.tsa_minus_sub.total"]
    ev = full[[c for c in dict.fromkeys(need + raw + panel) if c in full]].merge(final[["seq_id"] + s07], on="seq_id", how="left")
    ev.to_csv(FIX / "evaluated_inputs.csv", index=False)
    full[["seq_id", "z_tsa_vs_sub", "z_tsa_vs_prod", "flag_strict_boltz", "flag_placer_signal", "flag_tsa_over_sub", "flag_tsa_over_prod", "flag_shortlist"]].to_csv(FIX / "expected_flags.csv", index=False)
    fs = pd.read_csv(tables / "final_scoring.csv", low_memory=False)
    fs[["seq_id", "composite_z", "composite_rank", "final_rank", "n_metrics"]].to_csv(FIX / "expected_final_scoring.csv", index=False)
    tc = pd.read_csv(tables / "top_candidates.csv")
    tc[["seq_id", "panel_rank", "panel_n"]].to_csv(FIX / "expected_panel_rank.csv", index=False)
    cand = pd.read_csv(tables / "candidates.csv", usecols=["seq_id", "motif", "arm", "batch", "eval_set"], low_memory=False)
    ext = cand[(cand.eval_set == "extended") & cand.apply(lambda r: f"{r.motif}__{r.arm}__{r.batch}" in RUNS + ["hgnt__steer_tmol_ch18__b04"], axis=1)]
    ext.to_csv(FIX / "expected_extended.csv", index=False)


def compress():
    """gzip the larger CSVs (pandas reads .csv.gz directly)."""
    import gzip
    for f in FIX.rglob("*.csv"):
        if f.stat().st_size > 40_000:
            with open(f, "rb") as src, open(str(f) + ".gz", "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as dst:
                shutil.copyfileobj(src, dst)
            f.unlink()


if __name__ == "__main__":
    main(*sys.argv[1:3])
    compress()
