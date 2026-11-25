#!/usr/bin/env bash
#SBATCH --job-name=maxi-placer-arms
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=12:00:00
# PLACER on ALL arms (tsa/sub/prod holo + apo). Calls maxi_eval.run_placer EXPLICITLY (10_placer.sbatch had been pointed at
# run_placer_paul, a TSA-only copy, so it silently skipped sub/prod/apo). Resumable. Env: ARMS=sub,apo  NSAMPLES=20  MAXREP=3  ENTRIES=...
set -o pipefail
cd "${PIPE_ROOT:?set PIPE_ROOT to the suite directory}" && source env.sh
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
$PY_PLACER -m maxi_eval.run_placer --inputs work/inputs --out work/placer -n "${NSAMPLES:-20}" --modes holo,apo \
    ${ARMS:+--arms "$ARMS"} ${ENTRIES:+--entries "$ENTRIES"} --max-replicates "${MAXREP:-3}" --threads 2
