#!/usr/bin/env bash
#SBATCH --job-name=maxi-tmol
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=03:00:00
# tmol protocol: raw / relax / theozyme-constrained relax / substrate+product in the same pocket (+ optional apo, ablations).
# Resumable (one JSON per replicate). Optional: ENTRIES=...  MAXREP=3  EXTRA="--apo --ablate"  CST=/path/theozyme.cst
set -o pipefail
cd /mnt/work/projects/Maxi/pipeline && source env.sh
export TMPDIR="/tmp/maxi-$SLURM_JOB_ID"; mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
$PY_TMOL -m maxi_eval.run_tmol_entry --inputs work/inputs --out work/tmol --max-replicates "${MAXREP:-3}" \
    ${ENTRIES:+--entries "$ENTRIES"} ${CST:+--cst "$CST"} ${EXTRA:---apo --ablate}
