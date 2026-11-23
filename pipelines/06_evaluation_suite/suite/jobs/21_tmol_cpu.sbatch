#!/usr/bin/env bash
#SBATCH --job-name=maxi-tmol-cpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --time=06:00:00
# tmol on CPU, walking the entry list from the END; CPU-only: deliberately NO --gres (no GPU reserved). Counts 4 of the 8 team CPUs.
set -o pipefail
cd /mnt/work/projects/Maxi/pipeline && source env.sh
export OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 TMPDIR="/tmp/maxi-$SLURM_JOB_ID"; mkdir -p "$TMPDIR"; trap 'rm -rf "$TMPDIR"' EXIT
ENTS=$($PY_TMOL -m maxi_eval.shard_entries -k 1 2>/dev/null)
$PY_TMOL -m maxi_eval.run_tmol_entry --entries "$ENTS" --out work/tmol --apo --reverse
