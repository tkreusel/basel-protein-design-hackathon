#!/usr/bin/env bash
#SBATCH --job-name=kemp-overnight
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=12:00:00

# Unsteered RFD4 control for Maxi's steered run, end to end on ONE GPU:
#
#   1 GEN     RFD4, Maxi's with-Trp spec (conditioning only, no search)
#   2 SELECT  measure every design with Maxi's core criteria; take top-K
#   3 MPNN    LigandMPNN, motif fixed, KRHDENQSTYC forbidden near the base
#   4 BOLTZ   refold RFD4 + 4 MPNN sequences per design, 6NT bound, no MSA
#   5 SCORE   fold / site / contacts / base environment / access / ipSAE
#
# Every stage is idempotent (RFD4 skip_existing, Boltz skips finished
# predictions), so a requeue resumes instead of restarting. Stage wall times go
# to $RUN/STAGES.log.
#
# Settings that were NOT confirmed with the owner of the steered run are marked
# [ASSUMED] and are overridable by environment variable (see overnight/README.md).
#
# Env:  KEMP_N=250  KEMP_K=100  KEMP_MODE=topk|all|core  KEMP_SEED=201
#       KEMP_BS=4   KEMP_RECYCLE=<empty = repo default, as Maxi's run_r5.py>
#       KEMP_RUN=<existing run dir to RESUME>

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

N=${KEMP_N:-250}; K=${KEMP_K:-100}; MODE=${KEMP_MODE:-topk}
SEED=${KEMP_SEED:-201}                   # [ASSUMED] his seeds were 7, 11, 12
BS=${KEMP_BS:-4}                         # his run_r5.py uses batch_size=4
RECYCLE=${KEMP_RECYCLE:-}                # [ASSUMED] his run leaves this at the repo default
BOLTZ_SEED=${KEMP_BOLTZ_SEED:-0}         # [ASSUMED] his run_boltz_chunk.sbatch uses --seed 0
RUN=${KEMP_RUN:-$KEMP_RUNS/$(date +%Y-%m-%d)_unsteered-control_job${SLURM_JOB_ID}_${N}designs_seed${SEED}}
RFD4=/mnt/work/repos/RFD4-Proteina
PY="$KEMP_ROOT/env/kemp-orch/bin/python"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"
mkdir -p "$TMPDIR" "$RUN"; trap 'rm -rf "$TMPDIR"' EXIT
LOG="$RUN/STAGES.log"
stage() { echo "$(date +%H:%M:%S) $*" | tee -a "$LOG"; }
t0=$(date +%s); T=$t0
lap() { local n=$(date +%s); stage "  -> $1 took $(( (n-T)/60 )) min"; T=$n; }

stage "=== job $SLURM_JOB_ID on $(hostname) | N=$N K=$K mode=$MODE seed=$SEED bs=$BS recycle='${RECYCLE:-default}'"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | tee -a "$LOG"
echo "run dir: $RUN" | tee -a "$LOG"

# ---------------------------------------------------------------- 1 GEN
stage "STAGE 1 GEN"
count() { find "$RUN/raw_gen" -name '*_generated.cif' 2>/dev/null | wc -l; }
gen() {
  ( cd "$RFD4" && VIRTUAL_ENV="$RFD4/.venv-cu12" PATH="$RFD4/.venv-cu12/bin:$PATH" \
    "$RFD4/.venv-cu12/bin/rfd4-generate" \
      inputs="$KEMP_ROOT/inputs/maxi/kemp4_hg317-tips_with-trp.yaml" \
      out_dir="$RUN/raw_gen" num_replicates="$N" batch_size="$1" compile_lengths=null \
      seed="$SEED" skip_existing=true '~metrics' ${RECYCLE:+recycle_strategy=$RECYCLE} \
      > "$RUN/gen_bs$1.log" 2>&1 ) ; return $?
}
if [ "$(count)" -lt "$N" ]; then gen "$BS"; fi
if [ "$(count)" -lt "$N" ] && [ "$BS" -gt 2 ]; then
  stage "  only $(count)/$N designs after bs=$BS (OOM?) -> resuming with bs=2"
  grep -iE "out of memory|CUDA error" "$RUN/gen_bs$BS.log" | tail -2 | tee -a "$LOG"
  gen 2
fi
stage "  generated $(count)/$N designs"; lap GEN
[ "$(count)" -ge 1 ] || { stage "FATAL: no designs"; exit 1; }

# -------------------------------------------------------------- 2 SELECT
stage "STAGE 2 SELECT (CPU)"
"$PY" -I "$KEMP_ROOT/overnight/select_designs.py" "$RUN/raw_gen" "$RUN/selection" \
      --mode "$MODE" --k "$K" --procs 4 2>&1 | tee "$RUN/select.log" | tee -a "$LOG" | tail -14
[ -s "$RUN/selection/selection.csv" ] || { stage "FATAL: empty selection"; exit 1; }
lap SELECT

# ----------------------------------------------------------------- 3 MPNN
stage "STAGE 3 MPNN"
if [ ! -s "$RUN/sequences.csv" ]; then
  "$PY" -I "$KEMP_ROOT/overnight/mpnn_batch.py" "$RUN/selection/selection.csv" "$RUN" 2>&1 | tail -6 | tee -a "$LOG"
fi
[ -s "$RUN/sequences.csv" ] || { stage "FATAL: no sequences"; exit 1; }
stage "  sequences: $(( $(wc -l < "$RUN/sequences.csv") - 1 ))"; lap MPNN

# ---------------------------------------------------------------- 4 BOLTZ
stage "STAGE 4 BOLTZ"
"$PY" -I "$KEMP_ROOT/overnight/boltz_batch.py" "$RUN/sequences.csv" "$RUN" 125 | tee -a "$LOG"
BC="$RUN/boltz-cache"; mkdir -p "$BC"
for A in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do [ -e "$BC/$A" ] || ln -s "$BOLTZ_ASSETS/$A" "$BC/$A"; done
[ -e "$BC/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BC/mols"
for D in "$RUN"/boltz_in/chunk_*; do
  C=$(basename "$D"); [ -d "$D" ] || continue
  n_in=$(ls "$D"/*.yaml | wc -l)
  n_out=$(find "$RUN/boltz_out/$C" -name '*_model_0.cif' 2>/dev/null | wc -l)
  if [ "$n_out" -ge "$n_in" ]; then stage "  $C already complete ($n_out)"; continue; fi
  stage "  $C: $n_in inputs ($n_out done)"
  "$ENV_ROOT/env-boltz2/bin/boltz" predict "$D" --model boltz2 --cache "$BC" \
      --out_dir "$RUN/boltz_out/$C" --accelerator gpu --devices 1 --no_kernels \
      --diffusion_samples 1 --recycling_steps 3 --output_format mmcif \
      --num_workers 2 --preprocessing-threads 2 --seed "$BOLTZ_SEED" \
      >> "$RUN/boltz.log" 2>&1
  stage "  $C finished rc=$? ($(find "$RUN/boltz_out/$C" -name '*_model_0.cif' | wc -l) predictions)"
done
lap BOLTZ

# ---------------------------------------------------------------- 5 SCORE
stage "STAGE 5 SCORE (CPU)"
mkdir -p "$RUN/ipsae"
/home/ec2-user/.conda/envs/ipsae2/bin/ipsae "$RUN/boltz_out" --output-file "$RUN/ipsae/ipsae.tsv" \
    > "$RUN/ipsae/ipsae.log" 2>&1; stage "  ipsae rc=$?"
"$PY" -I "$KEMP_ROOT/overnight/score_refolds.py" "$RUN/selection/selection.csv" "$RUN/sequences.csv" \
      "$RUN/boltz_out" "$RUN/refold_scores.csv" --ipsae-dir "$RUN/ipsae" --procs 4 2>&1 | tee "$RUN/score.log" | tee -a "$LOG" | tail -8
lap SCORE
stage "=== ALL DONE in $(( ($(date +%s)-t0)/60 )) min ==="
