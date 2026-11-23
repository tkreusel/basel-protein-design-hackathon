#!/usr/bin/env bash
#SBATCH --job-name=kemp-rfshard
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=04:00:00

# Pin to the SMALL nodes only. This is not a preference, it is the only way
# to get two GPUs at once.
#
# The team quota is 8 AWS vCPUs. l4-1 is a g6.2xlarge = 8 vCPUs, i.e. the
# entire quota, so any job landing there blocks BOTH l4small nodes (and any
# teammate). The l4small pair is 4 + 4 = 8, which fits exactly.
#
# Slurm does not model the AWS quota, so left to itself it picks l4-1 first --
# precisely the choice that forfeits parallelism. Hence the explicit exclude.
# Memory is NOT the constraint: l4small has RealMemory=15564M and we ask for
# 14G (14336M), which fits.
#SBATCH --exclude=gpu-dy-l4big-1,gpu-dy-l4-1

# Refold one SHARD of a prepared boltz_in directory.
#
# Splitting lets both l4small nodes work the same batch in parallel (4 CPUs
# each = the team's 8-CPU budget). Each shard writes to its own out_dir; the
# analysis step reads them all together.
#
# Required: KEMP_WORK  = dir containing boltz_in/ and designs.json
#           KEMP_SHARD = this shard index, 0-based
#           KEMP_NSHARD= total shards

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

W="${KEMP_WORK:?set KEMP_WORK}"
SHARD="${KEMP_SHARD:?set KEMP_SHARD}"
NSHARD="${KEMP_NSHARD:?set KEMP_NSHARD}"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"
SH="$W/shard_${SHARD}"
mkdir -p "$TMPDIR" "$SH/in"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== node $(hostname) | job $SLURM_JOB_ID | shard $SHARD of $NSHARD ==="

# Deterministic round-robin split, so shards never overlap and together cover
# every input exactly once.
i=0
for f in $(ls "$W/boltz_in"/*.yaml | sort); do
    if [ $((i % NSHARD)) -eq "$SHARD" ]; then cp "$f" "$SH/in/"; fi
    i=$((i+1))
done
echo "this shard: $(ls "$SH/in"/*.yaml 2>/dev/null | wc -l) of $i inputs"

BC="$W/boltz-cache"; mkdir -p "$BC"
for ASSET in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do
    [ -e "$BC/$ASSET" ] || ln -s "$BOLTZ_ASSETS/$ASSET" "$BC/$ASSET"
done
[ -e "$BC/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BC/mols"

START=$(date +%s)
"$ENV_ROOT/env-boltz2/bin/boltz" predict "$SH/in" \
    --model boltz2 --cache "$BC" --out_dir "$SH/out" \
    --accelerator gpu --devices 1 --no_kernels \
    --sampling_steps 200 --recycling_steps 3 --diffusion_samples 1 \
    --num_workers 0 --preprocessing-threads 1 \
    --output_format mmcif --seed 1 2>&1 | tail -15
RC=$?
END=$(date +%s)
echo "boltz rc=$RC wall=$((END-START))s"
echo "predictions: $(find "$SH/out" -name '*_model_0.cif' 2>/dev/null | wc -l)"
echo "=== done ==="
