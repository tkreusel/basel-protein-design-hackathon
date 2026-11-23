#!/usr/bin/env bash
#SBATCH --job-name=maxi-pyrosetta-cpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --time=06:00:00
# 4 PyRosetta processes walking the entry list from the END (lock files avoid collisions). CPU-only: no --gres, no GPU reserved.
set -o pipefail
cd "${PIPE_ROOT:?set PIPE_ROOT to the suite directory}" && source env.sh
ENTS=$($PY_TMOL -m maxi_eval.shard_entries -k 1 2>/dev/null)
for i in 1 2 3 4; do
  sleep $((i*6))
  env PYTHONPATH=src $PY_ROSETTA -m maxi_eval.run_pyrosetta --entries "$ENTS" --out work/pyrosetta --relax-repeats 1 --nstruct 1 --reverse > work/full_ros_cpu_$i.log 2>&1 &
done
wait
