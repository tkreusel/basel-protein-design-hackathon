#!/usr/bin/env bash
#SBATCH --job-name=kemp-refold
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=04:00:00

# Boltz-2 refold survival test.
#
# Refolds ONLY the stage-1 success designs of one RFD4 run, from sequence
# alone, with 6NT bound, then re-gates on the same crystal-calibrated active
# site criteria. The first non-circular validation in the pipeline.
#
# Required: KEMP_RUN  = a finalized stage-1 run folder (has manifest.jsonl)
# Optional: KEMP_TAG  = short label for the output folder
#
# No `set -e`: each step reports and the next decides, so one bad prediction
# never kills the batch (plan 0.2).

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

RUN="${KEMP_RUN:?set KEMP_RUN to a finalized stage-1 run folder}"
TAG="${KEMP_TAG:-refold}"
PY="$KEMP_ROOT/env/kemp-orch/bin/python"
WORK="$KEMP_RUNS/raw_${TAG}_job${SLURM_JOB_ID}"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"
mkdir -p "$TMPDIR" "$WORK"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== node $(hostname) | job $SLURM_JOB_ID ==="
echo "source run: $RUN"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

echo
echo "=== STEP 1/3: prepare Boltz inputs (stage-1 successes only) ==="
"$PY" -I "$KEMP_ROOT/src/prepare_refold.py" "$RUN" "$WORK/boltz_in"
echo "prepare rc=$?"

echo
echo "=== STEP 2/3: Boltz-2 prediction ==="
BC="$WORK/boltz-cache"
mkdir -p "$BC"
for ASSET in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do
    [ -e "$BC/$ASSET" ] || ln -s "$BOLTZ_ASSETS/$ASSET" "$BC/$ASSET"
done
[ -e "$BC/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BC/mols"

START=$(date +%s)
# Flags copied from the organisers' verified example; --no_kernels kept
# deliberately. Single seed for now -- plan 6.6's 4-seed consistency arm is a
# separate, more expensive run.
"$ENV_ROOT/env-boltz2/bin/boltz" predict "$WORK/boltz_in" \
    --model boltz2 --cache "$BC" --out_dir "$WORK/boltz" \
    --accelerator gpu --devices 1 --no_kernels \
    --sampling_steps 200 --recycling_steps 3 --diffusion_samples 1 \
    --num_workers 0 --preprocessing-threads 1 \
    --output_format mmcif --seed 1 2>&1 | tail -25
RC=$?
END=$(date +%s)
echo "boltz rc=$RC wall=$((END-START))s"
echo "predictions: $(find "$WORK/boltz" -name '*_model_0.cif' 2>/dev/null | wc -l)"

echo
echo "=== STEP 3/3: score refolds and lay out the run folder ==="
OUTDIR="$KEMP_RUNS/$(date +%Y-%m-%d)_refold-${TAG}_job${SLURM_JOB_ID}"
"$PY" -I "$KEMP_ROOT/src/finalize_refold.py" "$KEMP_ROOT" \
    "$WORK/boltz" "$WORK/designs.json" "$OUTDIR"
echo "finalize rc=$?"
echo "=== done ==="
