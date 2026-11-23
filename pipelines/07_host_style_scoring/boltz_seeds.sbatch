#!/usr/bin/env bash
#SBATCH --job-name=kemp-seeds
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=08:00:00

# Multi-seed Boltz-2 refolds for the host-style "consistency" component.
# The hosts score reproducibility across Boltz-2 seeds as a proxy for
# preorganisation, so each candidate is refolded with SEVERAL independent seeds
# (separate runs, not diffusion samples of one trunk pass, which would share a
# trunk and overstate consistency).
#
# Work is ordered chunk-outer / seed-inner: a candidate gets ALL its seeds before
# the next chunk starts, and chunks are in priority order. If time runs out the
# finished candidates are complete. Boltz skips finished predictions, so a
# requeue resumes.
#
# Required: KEMP_CANDS (candidates.csv)  KEMP_OUT (output dir)
# Optional: KEMP_SEEDS="0 1 2 3 4"  KEMP_CHUNK=25  KEMP_STOP_AT=<epoch seconds>

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh
CANDS="${KEMP_CANDS:?}"; OUT="${KEMP_OUT:?}"; SEEDS="${KEMP_SEEDS:-0 1 2 3 4}"; CH="${KEMP_CHUNK:-25}"
STOP="${KEMP_STOP_AT:-0}"
PY="$KEMP_ROOT/env/kemp-orch/bin/python"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"; mkdir -p "$TMPDIR" "$OUT"; trap 'rm -rf "$TMPDIR"' EXIT
log() { echo "$(date +%H:%M:%S) $*" | tee -a "$OUT/SEEDS.log"; }
log "job $SLURM_JOB_ID on $(hostname) seeds=[$SEEDS] chunk=$CH"
"$PY" -I "$KEMP_ROOT/scoring/boltz_seed_inputs.py" "$CANDS" "$OUT" "$CH" | tee -a "$OUT/SEEDS.log"
BC="$OUT/boltz-cache"; mkdir -p "$BC"
for A in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do [ -e "$BC/$A" ] || ln -s "$BOLTZ_ASSETS/$A" "$BC/$A"; done
[ -e "$BC/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BC/mols"
for D in "$OUT"/boltz_in/chunk_*; do
  C=$(basename "$D")
  if [ "$STOP" -gt 0 ] && [ "$(date +%s)" -ge "$STOP" ]; then log "STOP_AT reached before $C; stopping cleanly"; break; fi
  for S in $SEEDS; do
    n_in=$(ls "$D"/*.yaml | wc -l); n_out=$(find "$OUT/seed_$S/boltz_results_$C" -name '*_model_0.cif' 2>/dev/null | wc -l)
    [ "$n_out" -ge "$n_in" ] && continue
    "$ENV_ROOT/env-boltz2/bin/boltz" predict "$D" --model boltz2 --cache "$BC" \
        --out_dir "$OUT/seed_$S" --accelerator gpu --devices 1 --no_kernels \
        --diffusion_samples 1 --recycling_steps 3 --output_format mmcif \
        --num_workers 2 --preprocessing-threads 2 --seed "$S" >> "$OUT/boltz.log" 2>&1
    log "$C seed $S rc=$? ($(find "$OUT/seed_$S" -name '*_model_0.cif' | wc -l) predictions in seed_$S)"
  done
done
log "multi-seed refolds done; scoring (CPU)"
"$PY" -I "$KEMP_ROOT/scoring/host_score.py" "$CANDS" "$OUT" "$OUT/scores" --top 25 2>&1 | tee -a "$OUT/SEEDS.log" | tail -32
log "ALL DONE"
