#!/usr/bin/env bash
#SBATCH --job-name=kemp-s1
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=02:00:00

# Stage 1: RFD4-Proteina backbone + sequence co-design around the HG3.17
# Kemp active site (inputs/kemp_rfd4.yaml).
#
# RFD4 co-designs sequence with structure, so its own output sequences are
# design arm A. Arm B re-designs the same backbones with LigandMPNN (stage 2).
#
# Env override order matters: paths.sh must come first so caches land on
# /mnt/work, then TMPDIR is pointed at node-local disk -- a /mnt/work TMPDIR
# caused a noisy "Device or resource busy" during multiprocessing teardown.

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

RFD4=/mnt/work/repos/RFD4-Proteina
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# Which condition spec to run. Lets one job script serve every theozyme
# variant (full motif, tip atoms, burial-conditioned, ...).
SPEC=${KEMP_SPEC:-/mnt/work/projects/tobi_v1/inputs/kemp_rfd4.yaml}
TAG=${KEMP_TAG:-$(basename "$SPEC" .yaml)}
N=${KEMP_N:-10}
BS=${KEMP_BS:-1}
SEED=${KEMP_SEED:-1}
OUT="$KEMP_RUNS/raw_${TAG}_job${SLURM_JOB_ID}"
mkdir -p "$OUT"

echo "=== node $(hostname) | job $SLURM_JOB_ID | n=$N bs=$BS seed=$SEED ==="
echo "spec: $SPEC"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo "out_dir: $OUT"

cd "$RFD4" || exit 1

# recycle_strategy=medium/midpeak is the documented schedule for enzyme /
# small-molecule (MCSA) tasks -- NOT the thinking_* aliases, which are
# front-loaded rather than task-aware (repo AGENTS.md).
#
# skip_existing=true makes the run resumable after a kill (plan 0.3).
# save_chain_of_interest=true is required for any later foldseek clustering.
START=$(date +%s)
VIRTUAL_ENV="$RFD4/.venv-cu12" PATH="$RFD4/.venv-cu12/bin:$PATH" \
  "$RFD4/.venv-cu12/bin/rfd4-generate" \
      inputs="$SPEC" \
      prevalidate_inputs=true \
      num_replicates="$N" batch_size="$BS" seed="$SEED" \
      recycle_strategy=medium/midpeak \
      skip_existing=true \
      save_chain_of_interest=true \
      out_dir="$OUT" 2>&1 | tail -60
RC=$?
END=$(date +%s)

echo
echo "=== rfd4-generate exit=$RC wall=$((END-START))s for n=$N ==="
echo "=== outputs ==="
find "$OUT" -name "*.cif" | head -20
echo "count: $(find "$OUT" -name '*.cif' | wc -l) cif files"
find "$OUT" -name "*.csv" | head
echo "=== done ==="
