# Cluster discovery (plan §2.0)

Measured 2026-11-21 on the head node and on GPU node `gpu-dy-l4-1` via job 2
(`jobs/probe_gpu.sbatch`, log `logs/probe-2.log`). Nothing here is assumed; every
row was observed. Re-run that job to re-verify.

## Scheduler

SLURM. `sbatch` / `squeue` / `sinfo` at `/opt/slurm/bin`.

One partition, `gpu` (default, no time limit), backed by **dynamic** nodes
(`gpu-dy-` prefix — they power on when a job is queued, so expect a cold-start
delay of a minute or two on the first submission):

| Node | Instance | Logical CPUs | RAM | Use? |
| --- | --- | --- | --- | --- |
| `gpu-dy-l4small-[1-2]` | g6.xlarge | 4 | 16 GiB | yes — the pair gives 2 GPUs |
| `gpu-dy-l4-1` | g6.2xlarge | 8 | 32 GiB | yes — **required for MMseqs2** (~28 GB RAM) |
| `gpu-dy-l4big-1` | g6.4xlarge | 16 | 64 GiB | **NO** — above our allocation; jobs queue forever |

**Team budget: 8 logical CPUs of GPU worker at a time** — either one
`l4-1` *or* both `l4small` nodes, not both at once. This is shared across
teammates, so check `squeue` before submitting.

### The budget is an AWS vCPU quota, and Slurm does not know about it

This is not a soft guideline. It is enforced by EC2, and the failure is
confusing:

```
gpu-dy-l4small-1    down#    (Code:VcpuLimitExceeded)
```

What happens: you submit two GPU jobs, Slurm puts one on `l4-1`
(g6.2xlarge = **8 vCPUs**, the whole quota) and tries to launch an `l4small`
(g6.xlarge = 4 vCPUs) for the other. AWS refuses, Slurm marks the node
`down#`, and the second job bounces back to `PD (BeginTime)` on a retry
backoff. The symptom looks like a stuck job sitting in `CF` with its clock
advancing.

**To actually get two GPUs, BOTH jobs must be pinned to the l4small pair**
(4 + 4 = 8). Touching `l4-1` at all forfeits parallelism:

```bash
sbatch --nodelist=gpu-dy-l4small-1 ...   # shard 0
sbatch --nodelist=gpu-dy-l4small-2 ...   # shard 1
```

Left to itself the scheduler picks `l4-1` first, which is exactly the
combination that cannot work.

**Also: an idle node still holds its vCPUs.** After a job finishes, `l4-1`
stays `idle` (powered on) through a cooldown before ParallelCluster
terminates it, and during that window a new `l4small` still cannot launch. So
a follow-up job may wait several minutes even with nothing running. Working GPU request syntax,
copied from the verified example:

```
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l4:1
#SBATCH --cpus-per-task=4      # 8 only on l4-1
#SBATCH --mem=12G              # 30G only on l4-1
```

Do **not** use `--constraint=g6.2xlarge` unless you specifically need the large
host; it pins you to the single contended node.

## GPU and driver — the finding that shaped the pipeline

Every node: **1 × NVIDIA L4, 23034 MiB VRAM**, driver **550.127.08 = CUDA 12.4**.

The driver is the binding constraint. CUDA *minor* version compatibility means a
torch built against any CUDA **12.x** runs on this driver, but a CUDA **13.x**
build does not (it needs driver ≥ 580).

| Environment | torch | CUDA build | GPU visible? |
| --- | --- | --- | --- |
| `env-rfd3` (provided) | 2.7.1 | 12.6 | yes |
| `env-ligandmpnn` (provided) | 2.7.1 | 12.6 | yes |
| `env-boltz2` (provided) | 2.7.1 | 12.6 | yes |
| `repos/RFD4-Proteina/.venv` as cloned | 2.14.1 | **13.2** | **no** |

RFD4-Proteina as installed is CPU-only on this cluster. See
[TOOL_NOTES.md](TOOL_NOTES.md) for the fix and its status.

## Filesystem

| Mount | Size | Free | Role |
| --- | --- | --- | --- |
| `/mnt/work` | 492 G | 459 G | writable workbench — **all caches and weights go here** |
| `/mnt/ref` | 295 G | 59 G | read-only reference library |
| `/` (incl. `$HOME`) | 100 G | 69 G | OS; small, shared — keep caches out |

No quota command found. `$HOME/.cache/huggingface` already holds 1.7 GB of RFD4
weights from the clone; `env/paths.sh` redirects `HF_HOME` and every other cache
under `/mnt/work/projects/tobi_v1/cache` so nothing further lands on `/`.

Head node and workers share `/mnt/ref` and `/mnt/work`.

## Internet egress — both head node and compute nodes

**Compute nodes have outbound internet.** Probed from `gpu-dy-l4-1`:
`pypi.org` → HTTP 200, `files.rcsb.org` → HTTP 200.

This retires the risk the plan calls "the single most common way a hackathon
pipeline dies at 2am". We do not need to pre-stage weights from a login node,
and a compute job may `pip install` and fetch PDBs directly. Still prefer
pre-staging for anything large, so a transient network failure does not waste a
queue slot.

## Modules, containers, conda

- No `module` system in use for these tools; no `apptainer`, `singularity` or
  `docker`. **Consequence:** RFD4's CUDA-12 and CPU-only *containers* are not an
  option, and any `*_refolding` metric tag in RFD4 (which needs an AF3
  container) is unavailable. Native installs only.
- `conda` and `mamba` present at `/mnt/ref/software/miniforge3`. Activate via
  `source "$CONDA_ROOT/etc/profile.d/conda.sh"` (done by `env/paths.sh`).
- `uv` is used inside the RFD4 repo for its own `.venv`.

## Provided software and weights (inherited, not duplicated)

`env/paths.sh` sources `/mnt/ref/hackathon/scripts/reference-paths.sh`, which
exports:

| Variable | Path |
| --- | --- |
| `ENV_ROOT` | `/mnt/ref/software/envs` (`env-rfd3`, `env-ligandmpnn`, `env-boltz2`) |
| `LIGANDMPNN_ROOT` | `/mnt/ref/software/LigandMPNN` |
| `MMSEQS` | `/mnt/ref/software/mmseqs2-18-8cc5c/bin/mmseqs` |
| `MMSEQS_DB` | `/mnt/ref/data_weights/mmseqs2_db/targetDB_gpu` |
| `RFD3_CHECKPOINT` | `/mnt/ref/data_weights/rfd3/rfd3_latest.ckpt` |
| `LIGANDMPNN_CHECKPOINT` | `/mnt/ref/data_weights/ligandmpnn/ligandmpnn_v_32_010_25.pt` |
| `BOLTZ_ASSETS` | `/mnt/ref/data_weights/boltz2` |

**Not installed anywhere:** tmol, PLACER, BioEmu-2, Foldseek, TM-align/US-align,
RFdiffusion2. All are reachable on PyPI/GitHub from both node classes.

## Open items

- **`SUBMISSION.md` is missing.** `GUIDE.md` states it sits "alongside this
  guide" and defines the FASTA format, the catalytic-residue header fields, the
  eligibility rules and a local validator. `/mnt/ref/hackathon/docs/` contains
  only `GUIDE.md` and `GUIDE.html`. **Ask an organiser** — it specifies the
  deliverable, and plan §6.7's compliance gate should be written against it
  rather than against our reading of the brief.
- The organisers' setup verification (job 1) was still RUNNING at 17+ min when we
  started, occupying `gpu-dy-l4-1`. Confirm it reached
  `VERIFIED: ...` via `bash /mnt/work/kemp_example/status.sh` before trusting the
  provided example as a reference.
