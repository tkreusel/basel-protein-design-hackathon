#!/usr/bin/env bash
#SBATCH --job-name=kemp-s25
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=03:00:00

# Stages 2 + 5: sequence design (arm B) and Boltz-2 refolding of BOTH arms.
#
#   arm A = RFD4's own co-designed sequences (from stage 1)
#   arm B = LigandMPNN re-design of the same RFD4 backbones
#
# Both refolded with 6NT so the comparison is like-for-like. Pass the stage-1
# run directory as KEMP_STAGE1.
#
# No `set -e`: each step reports and the next decides whether it can proceed,
# so one bad design never kills the batch (plan 0.2).

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

STAGE1="${KEMP_STAGE1:?set KEMP_STAGE1 to the stage-1 run dir}"
NSEQ="${KEMP_NSEQ:-4}"
STAGE="$KEMP_RUNS/designs_$SLURM_JOB_ID"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"
mkdir -p "$TMPDIR" "$STAGE"
trap 'rm -rf "$TMPDIR"' EXIT

PY="$KEMP_ROOT/env/kemp-orch/bin/python"
echo "=== node $(hostname) | job $SLURM_JOB_ID ==="
echo "stage1: $STAGE1"
echo "stage:  $STAGE"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

echo
echo "=== STEP 1/4: extract arm A + backbone PDBs from RFD4 output ==="
"$PY" -I "$KEMP_ROOT/src/extract_designs.py" "$KEMP_ROOT" "$STAGE1" "$STAGE"
echo "extract rc=$?"

echo
echo "=== STEP 2/4: LigandMPNN re-design (arm B), catalytic residues fixed ==="
"$PY" -I "$KEMP_ROOT/src/run_stage2_ligandmpnn.py" "$KEMP_ROOT" "$STAGE" "$NSEQ"
echo "stage2 rc=$?"

echo
echo "=== STEP 3/4: prepare Boltz-2 inputs for both arms ==="
"$PY" -I "$KEMP_ROOT/src/prepare_boltz_both.py" "$STAGE" "$STAGE/boltz_in"
echo "prepare rc=$?"

echo
echo "=== STEP 4/4: Boltz-2 refolding, both arms ==="
BOLTZ_CACHE="$STAGE/boltz-cache"
mkdir -p "$BOLTZ_CACHE"
# Boltz needs its assets visible in the cache dir; symlink rather than copy.
for ASSET in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do
    [ -e "$BOLTZ_CACHE/$ASSET" ] || ln -s "$BOLTZ_ASSETS/$ASSET" "$BOLTZ_CACHE/$ASSET"
done
[ -e "$BOLTZ_CACHE/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BOLTZ_CACHE/mols"

# Flags copied from the organisers' verified example; --no_kernels is kept
# deliberately. 4 seeds would be plan 6.6's seed-consistency arm -- here we
# take 1 seed first to establish wall-clock, then scale.
"$ENV_ROOT/env-boltz2/bin/boltz" predict "$STAGE/boltz_in" \
    --model boltz2 --cache "$BOLTZ_CACHE" --out_dir "$STAGE/boltz" \
    --accelerator gpu --devices 1 --no_kernels \
    --sampling_steps 200 --recycling_steps 3 --diffusion_samples 1 \
    --num_workers 0 --preprocessing-threads 1 \
    --output_format mmcif --seed 1 2>&1 | tail -40
echo "boltz rc=$?"

echo
echo "=== outputs ==="
echo "arm A: $(grep -c '^>' "$STAGE/armA_rfd4_sequences.fasta" 2>/dev/null || echo 0) sequences"
echo "arm B: $(grep -c '^>' "$STAGE/armB_ligandmpnn_sequences.fasta" 2>/dev/null || echo 0) sequences"
echo "boltz predictions: $(find "$STAGE/boltz" -name '*.cif' 2>/dev/null | wc -l)"
echo "stage dir: $STAGE"
echo "=== done ==="
