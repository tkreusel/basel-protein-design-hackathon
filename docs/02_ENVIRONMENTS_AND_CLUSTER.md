# 02 — Environments, cluster and path assumptions

## Cluster in one table (details: [reference/CLUSTER.md](reference/CLUSTER.md))

| Item | Fact |
|---|---|
| Scheduler | Slurm, one partition `gpu` backed by dynamic AWS nodes (1–2 min cold start) |
| Nodes | `gpu-dy-l4small-[1,2]` (4 vCPU / 16 GB), `gpu-dy-l4-1` (8 vCPU / 32 GB), `gpu-dy-l4big-1` (16 vCPU / 64 GB — **outside our allocation, never use**) |
| GPU | 1 × NVIDIA L4 (23 GB) per node, driver 550.127.08 = **CUDA 12.4** (CUDA 13 torch builds cannot see the GPU) |
| Team budget | **8 vCPU of GPU workers at once**: one `l4-1` *or* both `l4small`; an AWS quota, not a Slurm one (`VcpuLimitExceeded` shows as a `down#` node). Pin pairs with `--nodelist=gpu-dy-l4small-1/2`. |
| Head node | no GPU; fine for CPU stages (tmol, PyRosetta, scoring, Foldseek, US-align). Shared with the team. |
| Filesystems | `/mnt/work` (writable workbench, put every cache here), `/mnt/ref` (read-only reference: tools, weights, UniRef90), `/` small |
| Internet | compute and head nodes have **outbound** access (PyPI, files.rcsb.org). We only ever *downloaded*; no project data was sent out |
| Typical request | `--partition=gpu --gres=gpu:l4:1 --cpus-per-task=4 --mem=14G` (use `--mem=30G` and node `l4-1` for MMseqs2 vs UniRef90) |

## Environments

| Name | Location (as on the cluster) | Used for | How it was made |
|---|---|---|---|
| `kemp-orch` | `<project>/env/kemp-orch` | geometry, scoring, ranking, eligibility glue, tests (CPU) | `/mnt/ref/software/envs/env-rfd3/bin/python -m venv --system-site-packages env/kemp-orch` then `pip install pytest` (so it inherits biotite 1.4.0, numpy, pandas, scipy from `env-rfd3`, Python 3.12) |
| RFD4 + tmol | `/mnt/work/repos/RFD4-Proteina/.venv-cu12` | RFD4 generation (GPU), tmol (GPU or CPU), analysis | provided clone; **use `.venv-cu12`, the default `.venv` is CUDA-13 and blind to this GPU** |
| `env-rfd3` | `$ENV_ROOT/env-rfd3` | RFD3 fallback (not used for final designs) | provided |
| `env-ligandmpnn` | `$ENV_ROOT/env-ligandmpnn` | LigandMPNN | provided (checkpoint `v_32_010_25`) |
| `env-boltz2` | `$ENV_ROOT/env-boltz2` | Boltz-2 (`boltz predict`) | provided; weights in `$BOLTZ_ASSETS` |
| `placer_env` | `~/.conda/envs/placer_env` | PLACER (`/mnt/work/repos/placer`) | created by the evaluation-suite owner; see the PLACER repo README |
| `ipsae2` | `~/.conda/envs/ipsae2` (`ipsae` binary) | ipSAE2 | created by the suite owner; `/mnt/work/repos/ipsae2` |
| `pyrosetta` | `/mnt/work/envs/pyrosetta` | PyRosetta 2026.47, CPU only | official installer **with the env's `bin` first on PATH** (see [07](07_LESSONS_AND_PITFALLS.md)); academic licence |
| `bioemu` (conda) | conda env `bioemu`, repo `/mnt/work/repos/bioemu2` | BioEmu-2 ensembles | provided repo |
| MMseqs2 | `/mnt/ref/software/mmseqs2-18-8cc5c/bin/mmseqs` | PDB-seqres and UniRef90 identity | provided; UniRef90-2022-05 **GPU** DB (`$MMSEQS_DB`) |
| Foldseek / US-align | `/mnt/work/repos/foldseek_dl/foldseek/bin/foldseek`, `/mnt/work/repos/usalign/USalign` | TM-score eligibility | downloaded binaries |

The tool paths the organisers provide come from `/mnt/ref/hackathon/scripts/reference-paths.sh`; `core/env/paths.sh` sources it
(and pins every cache under the project, never `$HOME`) so a reference update cannot silently desync our scripts.

### The one trap
`/mnt/work/repos/RFD4-Proteina/.venv` (default name) has torch `cu132` and cannot use the GPU. A marker file inside says so.
Do not run the repo's `setup.sh`/`make setup` again.

### RFD4 invocation rules that bit us
Run from the RFD4 repo root (the dir with `.project-root`); Hydra syntax is `key=value`, never `--key value`; `inputs` must be
`.cif`/`.cif.gz`/`.bcif` (a `.pdb` is rejected); `recycle_strategy=medium/midpeak` is the documented enzyme schedule;
use `skip_existing=true` so a killed run resumes; set `TMPDIR` to node-local disk in GPU jobs. Full recipe, including the theozyme
conditioning that took three attempts: [reference/TOOL_NOTES.md](reference/TOOL_NOTES.md), [reference/USAGE.md](reference/USAGE.md)
(shipped as `docs/reference/USAGE.md`).

## Path map: what to change to re-run elsewhere

Scripts were copied as run. Hard-coded locations you will need to repoint:

| Appears in scripts as | Meaning | Change to |
|---|---|---|
| `/mnt/work/projects/tobi_v1` (also via `$KEMP_ROOT`) | our project root (runs, inputs, logs, caches) | your checkout of this repo + a writable scratch dir for `runs/`, `cache/` |
| `/mnt/work/projects/Maxi/...` | evaluation-suite and campaign data (his folder) | `pipelines/06_evaluation_suite/suite` and wherever you keep his result tables |
| `/mnt/work/repos/{RFD4-Proteina,bioemu2,placer,ipsae2,usalign,foldseek_dl}` | tool installs | your installs |
| `/mnt/ref/...` | read-only reference library (weights, UniRef90, organisers' scripts) | the equivalent locations |
| `$ENV_ROOT`, `$MMSEQS`, `$MMSEQS_DB`, `$BOLTZ_ASSETS`, `$LIGANDMPNN_*` | defined by `reference-paths.sh` | define them yourself |
| `/home/ec2-user/.conda/envs/...` | conda envs of the suite owner | your conda prefix |

`core/env/paths.sh` is the single place for our own variables; most `pipelines/*/*.sbatch` start with
`source /mnt/work/projects/tobi_v1/env/paths.sh` — replace that line with the path to your copy.

## Reference data we had to build

| Data | How | Used by |
|---|---|---|
| PDB seqres MMseqs2 DB | `pdb_seqres.txt.gz` downloaded from wwPDB, filtered to `mol:protein`, `mmseqs createdb` | eligibility stage 2 |
| UniRef90-2022-05 GPU DB | provided by the organisers (`$MMSEQS_DB`) | eligibility stage 3, GPU job `uniref_search.sbatch` |
| Foldseek DB of PDB 30 %-cluster representatives (40,064 entries) | the official Foldseek PDB download was blocked (HTTP 429), so we downloaded mmCIFs of one representative per 30 % sequence cluster from `files.rcsb.org` (entry list: `pipelines/08_eligibility/db_recipe/pdb30_entries.txt`, helper `fetch_one.sh`) and ran `foldseek createdb` | eligibility stage 4 |
| 6NT/H5J/6VP ligand parameters | built from CCD codes by the suite (`build_ligands`, `build_rosetta_params`) | tmol, PyRosetta |
| Literature enzyme table | `/mnt/work/manually_created_dataset.csv` (name → kcat/KM etc.) and the suite's `kemp_eval_dataset_full.csv` | literature calibration |

## Pinned versions

`docs/env_versions/` holds `pip freeze` for each environment we could introspect and `tool_repo_commits.txt` (git commits of the RFD4-Proteina, bioemu2, placer and ipsae2 clones; MMseqs2 `18-8cc5c`, Foldseek build hash). Model weights are the organisers' copies under `/mnt/ref/data_weights` (Boltz-2 `boltz2_conf.ckpt`, LigandMPNN `v_32_010_25`).
