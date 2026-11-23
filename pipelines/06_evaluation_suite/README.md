# 06 Evaluation suite (ipSAE2, PLACER, tmol, PyRosetta)

`suite/` holds the code, configuration, jobs and tests of the evaluation suite; `PIPELINE.md` describes the stages, the tmol protocol per replicate, what is native and what a proxy,
and findings along the way; `RESULTS.md` and `PLACER_RESULTS.md` report the validation runs on the HG3 crystals.

Changes made for the large campaign (`../02_rfd4_steering`):

- `env.sh` is relocatable: `PIPE_ROOT` is the folder it lives in, and the interpreters (`PY_TMOL`, `PY_PLACER`, `PY_ROSETTA`, `IPSAE_BIN`, `PLACER_ROOT`) must be exported by the caller (the campaign does this in `config/site.env`). The Slurm scripts no longer contain absolute paths.
- `config/entries.design_run.yaml`: template for a run with no crystal entries, where the proteins come only from a Boltz output folder.
- `aggregate.py`, `run_ipsae.py`: Boltz ipDE, ligand ipTM and ipLDDT are carried into the summary, and min and max over samples are kept.
- `analyze_placer.py`: one model that cannot be analysed no longer costs the run its PLACER summary.
- `geometry.py`: missing benzo or wall-ring atoms (PLACER pocket crops) are skipped instead of raising.
- `run_pyrosetta.py`: a failed Sc calculation no longer discards dG_sep, packstat and the hydrogen-bond terms.
- `run_placer.py`: `--arms` accepts `+` as separator (a comma splits `sbatch --export` lists).

```
Boltz-2 refolds (6NT / H5J / 6VP per protein)
   └─ jobs/00_prepare_and_ipsae.sh   standardise inputs (protein chain A, ligand chain B), ipSAE2 for every arm        CPU, minutes
        ├─ jobs/10_placer*.sbatch    PLACER holo (each ligand) + apo ensembles, 20–50 models                           GPU
        ├─ jobs/20_tmol / 21_tmol_cpu  tmol protocol per replicate (raw, rlx, rlxcst, tsa_pk/sub/prod, …)             head-node CPU, sharded
        ├─ jobs/40_pyrosetta / 41     relax, interface dG_sep, Sc (unreliable), packstat, BUNS, ΔΔG                      CPU
        ├─ jobs/30_placer_analysis    ensemble statistics + tmol electric field                                         GPU/CPU
        └─ aggregate → work/report/{summary_entries.csv, summary_replicates.csv, report.md}
jobs/run_full_suite.sh --root <boltz folder> [--submit]    dry-run by default; prints readiness, queue, plan
```

`our_additions/` holds what was added to run the suite on the hackathon finalists:

| File | Role |
|---|---|
| `stage_finalists.py` | copies the 6NT refolds of the finalists into `boltz_root/<id>/` and writes Boltz inputs for the substrate (H5J) and product (6VP) arms |
| `stage_extra.py` | the same for the extra pool designs never evaluated |
| `arms_then_placer.sbatch`, `arms2.sbatch` | Boltz arm refolds → prepare + ipSAE2 → PLACER; second batch of arms only |
| `geometry_wall_ring_fallback.patch` | earlier fix for the same crash (superseded by the `geometry.py` change listed above) |
| `plot_best.patch`, `env.sh.as-used`, `finalists.csv`, `finalists_report.md` | plotting tweak, environment as used, the 12 finalists, the suite report for them |

Caveats: composite/rank are uncalibrated heuristics; tmol relax is stochastic (~5 REU); Sc unreliable; PLACER is a weak, noisy feature
(HG3.17 crystal holo ≈ 0.05; within-backbone variance ≈ 50 %); 6NT is a neutral TSA analogue and 6VP in the CCD is the neutral phenol.
