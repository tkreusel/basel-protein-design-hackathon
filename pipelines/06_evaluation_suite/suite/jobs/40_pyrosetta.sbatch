#!/usr/bin/env bash
#SBATCH --job-name=maxi-pyrosetta
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --time=04:00:00
# PyRosetta is CPU-only: deliberately NO --gres, so no GPU is reserved. (Single-threaded per replicate; ~3.5-10 min each.)
# Optional: ENTRIES=...  MAXREP=3  EXTRA="--ablate --nstruct 2"  CST=/path/theozyme.cst
set -o pipefail
cd "${PIPE_ROOT:?set PIPE_ROOT to the suite directory}" && source env.sh
export PYTHONPATH="$PIPE_ROOT/src"
"$PY_ROSETTA" -m maxi_eval.run_pyrosetta --inputs work/inputs --out work/pyrosetta \
    --max-replicates "${MAXREP:-3}" ${ENTRIES:+--entries "$ENTRIES"} ${CST:+--cst "$CST"} ${EXTRA:---ablate --nstruct 2}
