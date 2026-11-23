#!/usr/bin/env bash
#SBATCH --job-name=maxi-placer
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=12:00:00
# PLACER holo+apo ensembles. ~1-3 s/model on an L4 -> 50 models x 2 modes ~ 3-5 min per replicate.
# Resumable (existing CSVs are skipped). Optional: ENTRIES=4BS0,7K4Q  NSAMPLES=100  MAXREP=3
set -o pipefail
cd /mnt/work/projects/Maxi/pipeline && source env.sh
export TMPDIR="/tmp/maxi-$SLURM_JOB_ID"; mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
$PY_PLACER -m maxi_eval.run_placer_paul --inputs work/inputs --out work/placer \
	-n "${NSAMPLES:-50}" --modes holo,apo ${ENTRIES:+--entries "$ENTRIES"} --max-replicates "${MAXREP:-3}"
