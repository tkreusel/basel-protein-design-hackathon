# How to use the environments

Every command below was run and verified on 2026-11-21. Project root is
`/mnt/work/projects/tobi_v1`.

## The one-liner

```bash
cd /mnt/work/projects/tobi_v1
source env/paths.sh
```

`paths.sh` is idempotent and safe to source repeatedly. It sets the tool paths
(inherited from the organisers' `reference-paths.sh`), puts every cache on
`/mnt/work`, and — importantly — points `HF_HOME` at the RFD4 weights.

## Which environment for what

| Env | Path | Use for |
|---|---|---|
| **kemp-orch** | `env/kemp-orch/bin/python` | geometry, scoring, benchmark, tests. CPU only, runs on the head node. |
| **RFD4 + tmol** | `/mnt/work/repos/RFD4-Proteina/.venv-cu12/bin/python` | backbone generation (stage 1), physics scoring (stage 3). **GPU — submit via Slurm.** |
| env-rfd3 | `$ENV_ROOT/env-rfd3/bin/python` | RFD3 fallback for stage 1 |
| env-ligandmpnn | `$ENV_ROOT/env-ligandmpnn/bin/python` | sequence design (stage 2) |
| env-boltz2 | `$ENV_ROOT/env-boltz2/bin/boltz` | fold + interface (stage 5) |

### The one trap

`/mnt/work/repos/RFD4-Proteina/.venv` — **the default-named venv — is broken.**
It has torch `cu132` and cannot see the GPU on this driver. Always use
**`.venv-cu12`**. Do not run `setup.sh` or `make setup` again; they target the
broken one. A `DO_NOT_USE_READ_THIS.txt` marker sits inside it.

## 1. Orchestration / scoring (CPU, head node)

```bash
cd /mnt/work/projects/tobi_v1
source env/paths.sh

# the Milestone 1 benchmark table
env/kemp-orch/bin/python -I src/milestone1.py .

# score any motif or design against the calibrated thresholds
env/kemp-orch/bin/python -I src/check_motif.py . inputs/HG3_17_kemp_motif.pdb

# the test suite
PYTHONPATH=src env/kemp-orch/bin/python -m pytest tests/ -q
```

`PYTHONPATH=src` is needed for `pytest` only; `paths.sh` already exports it for
everything else. Use `python -I` for scripts that read downloaded data, so the
interpreter ignores the current directory.

## 2. RFD4 (GPU — always through Slurm)

**The head node has no GPU.** `nvidia-smi` fails there. Never run generation
interactively on it.

Two rules that will bite you otherwise:

1. **Run from the RFD4 repo root** (the directory with `.project-root`).
2. **Hydra syntax: `key=value`, never `--key value`.** `--inputs x` is rejected.

A working invocation, copied from the verified smoke job
(`jobs/rfd4_smoke.sbatch`, log `logs/rfd4smoke-7.log`):

```bash
RFD4=/mnt/work/repos/RFD4-Proteina
cd "$RFD4"
VIRTUAL_ENV="$RFD4/.venv-cu12" PATH="$RFD4/.venv-cu12/bin:$PATH" \
  "$RFD4/.venv-cu12/bin/rfd4-generate" \
      inputs=<your.cif or spec.yaml> \
      num_replicates=8 batch_size=1 \
      recycle_strategy=medium/midpeak \
      out_dir="$KEMP_RUNS/myrun" \
      seed=1 skip_existing=true
```

- `recycle_strategy=medium/midpeak` is the documented schedule for enzyme /
  small-molecule (MCSA) tasks. Do not use the `thinking_*` aliases — they are
  front-loaded, not task-aware.
- Confirm your longest input works at `batch_size=1` before raising it.
- `skip_existing=true` makes a run resumable after a kill.
- `inputs` must be `.cif`, `.cif.gz` or `.bcif` — **a `.pdb` is rejected.** Use
  `src/motif_to_cif.py` to convert.

Submit it:

```bash
sbatch --chdir=/mnt/work/projects/tobi_v1 \
       --output=/mnt/work/projects/tobi_v1/logs/myrun-%j.log \
       jobs/your_job.sbatch
squeue
tail -f logs/myrun-<JOBID>.log
```

### Slurm resources that work

```
#SBATCH --partition=gpu
#SBATCH --gres=gpu:l4:1
#SBATCH --cpus-per-task=4      # 8 only on the g6.2xlarge
#SBATCH --mem=14G              # 30G only on the g6.2xlarge
#SBATCH --time=00:50:00
```

Never target `l4big` — it is above our allocation and queues forever. The team
budget is 8 logical CPUs of GPU worker **total**, shared with teammates, so
check `squeue` first (a `rfd4test` job from someone else has already appeared).

Set `TMPDIR` to node-local disk in GPU jobs. Pointing it at `/mnt/work` produced
a harmless-but-noisy `OSError: Device or resource busy` during multiprocessing
teardown.

## 3. Interactive poking at RFD4 (no GPU work)

```bash
source /mnt/work/repos/RFD4-Proteina/.venv-cu12/bin/activate
cd /mnt/work/repos/RFD4-Proteina
python -c "import rfd4_proteina, tmol; print('ok')"
deactivate
```

Fine for imports, config composition and input validation. Anything that
touches the GPU goes through Slurm.

## Timing, for budgeting

One RFD4 design at `batch_size=1`, including model load: **~87 s** on an L4.
Model load dominates, so batch aggressively — plan §7's "batch, never one job
per design" applies directly.
