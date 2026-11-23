# 08 Provenance and data locations

This repository holds code, configuration, documentation and small result tables. Raw outputs are large and stay on the clusters.

## Where each repo folder came from

| Repo path | Original location (cluster) | Notes |
|---|---|---|
| `core/src`, `core/configs`, `core/tests`, `core/ligands`, `core/env/paths.sh` | `/mnt/work/projects/tobi_v1/{src,configs,tests,ligands,env}` | `core/tests/test_geometry.py` modified: honours `KEMP_ROOT` so tests run from a workspace |
| `pipelines/01_backbones_sequence_refold/{jobs,inputs}` | `tobi_v1/{jobs,inputs}` | stage 1–5 job scripts and RFD4 specs |
| `pipelines/02_rfd4_steering/{config,src,jobs,tests,figures,inputs,rfd4_patch,results_summary}` | `~/MPIP/kemp_campaign` on the DAIS cluster (code that ran the campaign, refactored; see its README) | gates, seeds, configs and rankings reproduce the campaign tables (`tests/`) |
| `archive/early_rounds_r4_r5` | `tobi_v1/imports/steering_handoff_20261121/` (round 4/5 steering code and metrics as received on 2026-11-21) | superseded by `pipelines/02_rfd4_steering` |
| `pipelines/04_overnight_unsteered_control` | `tobi_v1/overnight` | |
| `pipelines/05_partial_diffusion` | `tobi_v1/partial` | |
| `pipelines/06_evaluation_suite/suite` | `/mnt/work/projects/Maxi/pipeline` (code, config, jobs, tests, docs; `work/` and logs excluded) | with the campaign fixes listed in `pipelines/06_evaluation_suite/README.md` |
| `pipelines/06_evaluation_suite/our_additions` | `tobi_v1/maxi_eval/{stage_finalists.py,stage_extra.py,jobs/arms*.sbatch,finalists.csv,env.sh}` + two patches | `*.patch` = our `diff -u` against the suite's `geometry.py` and `plot_best.py` |
| `pipelines/07_host_style_scoring` | `tobi_v1/scoring` | + candidate lists |
| `pipelines/08_eligibility` | `tobi_v1/eligibility` (scripts, candidate lists v5/v6, final `eligibility.csv`, PDB entry list) | databases not included |
| `pipelines/09_bioemu` | `tobi_v1/bioemu/*.py`, `tobi_v1/bioemu/literature`, `/mnt/work/projects/Wasabi/run_bioemu.slurm`, `/mnt/work/projects/bioemu_traj/*.py` | trajectories and per-metric JSONs not included |
| `pipelines/11_selection_and_submission` | `tobi_v1/submission`, `/mnt/work/projects/departure` | scripts, ranking tables, FASTA, mapping |
| `results/final_top20_structures` | `/mnt/work/projects/departure/final_top20_structures` | |
| `docs/reference/*` | `tobi_v1/docs/*.md`, `tobi_v1/runs/README.md`, hosts' `SUBMISSION.md`, our funnel overview | verbatim originals, only cluster account strings redacted |
| `docs/env_versions/` | `pip freeze` of each environment on 2026-11-23 | |

## What is deliberately not in the repo (and how to get it back)

| Data | Size / where | How to regenerate |
|---|---|---|
| Model weights (RFD4, Boltz-2, LigandMPNN, PLACER, BioEmu) | organisers' `/mnt/ref/data_weights`, `~/.cache`, HF | provided; see `core/env/paths.sh` |
| Raw RFD4 outputs, Boltz outputs (`runs/`) | tens of GB under `tobi_v1/runs` (34 GB project total) | rerun with the seeds in [03](03_REPRODUCING.md) |
| Campaign data (10,500 backbones, evaluation tables, sequences) | scratch of the DAIS cluster; one archive per run (`jobs/archive_runs.sh`); summary tables in `pipelines/02_rfd4_steering/results_summary/` | rerun the pipeline with the seeds in `config/arms.yaml`; `python -m kemp.export_results` rebuilds the tables from the runs |
| Evaluation-suite results (`work/`) | `Maxi/pipeline/work`, `tobi_v1/maxi_eval/work` | rerun `run_full_suite.sh` on Boltz outputs |
| Literature crystal structures | public PDB | `scripts/fetch_benchmark.sh` |
| PDB seqres DB, UniRef90 GPU DB, Foldseek PDB30 DB | `tobi_v1/eligibility/db` (1.5 GB), organisers' UniRef90 | recipes in [02](02_ENVIRONMENTS_AND_CLUSTER.md) |
| BioEmu trajectories | `bioemu_traj/`, `Wasabi/outputs` | `run_bioemu.slurm` |
| Literature activity table | `/mnt/work/manually_created_dataset.csv`, `kemp_eval_dataset_full.csv` | supplied by the team; needed by the literature factor |
| AWS credentials / keys | never stored here | — |

## Literature-structures index
`docs/reference/LITERATURE_STRUCTURES.md` maps 97 Kemp structures in the PDB to the literature designs, with 6NT presence verified from coordinates,
and `core/src/find_kemp_structures.sh` reproduces the RCSB search behind it.

## Reference documents (read these before re-deriving anything)

| Doc | Why |
|---|---|
| `reference/TOOL_NOTES.md` | tool interfaces and the RFD4 conditioning recipe |
| `reference/USAGE.md`, `reference/CLUSTER.md` | environments, commands, measured cluster facts |
| `reference/MILESTONE1.md` | geometry stack validated on the HG series, with the failing criterion reported honestly |
| `reference/HIS_ZN_MOTIFS.md` | His/Zn motif exploration (not used for the final designs) |
| `reference/SUBMISSION_hosts_spec.md` | the hosts' rules and criteria |
| `reference/RUN_FOLDER_CONVENTION.md` | how run folders are laid out and what "success" means per stage |
| `reference/pipeline_funnel_overview.md` | the pitch funnel (superseded by `06`, kept verbatim) |

The original staged plan (`KEMP_PIPELINE_PLAN.md`) and a comparison of the round 4/5 steering code with our pipeline were removed from the working tree; both are in the git history (commit d7ab76f).
