#!/bin/bash
# Whole-node job (needs KEMP_PIPELINE exported, which jobs/submit.sh does): one worker per GPU (worker.sh). Partition, GPU count and walltime come from the command line or config/site.env:
#   sbatch -p "$NODE_PARTITION" --gres=gpu:"$NODE_GPUS" -t "$NODE_TIME" --exclusive --mem=0 --export=ALL,TASKS=<task file> jobs/node_worker.sbatch
# Whole nodes were the only thing that started quickly on the shared queue; single-GPU array jobs waited for many hours.
#SBATCH --exclusive
#SBATCH --mem=0
D=${KEMP_PIPELINE:?export KEMP_PIPELINE=<this pipeline directory> before sbatch}/jobs; source "$D/lib.sh"
nvidia-smi -L | head -"${NODE_GPUS:-8}"
for g in $(seq 0 $(( ${NODE_GPUS:-8} - 1 ))); do
  bash "$D/worker.sh" "$g" > "$KEMP_SCRATCH/logs/$CAMP/W_${SLURM_JOB_ID}_g$g.log" 2>&1 & sleep 3
done
# stay alive while any worker (including the fanned-out copies) runs
while pgrep -f "worker.sh" > /dev/null; do sleep 60; done
