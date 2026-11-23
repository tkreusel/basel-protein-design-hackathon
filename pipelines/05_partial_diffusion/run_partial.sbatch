#!/usr/bin/env bash
#SBATCH --job-name=kemp-partial
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=06:00:00

# Partial diffusion around our best designs, then the same MPNN -> Boltz -> host-style scoring.
#
#   0 PARENTS  rank everything scored so far on composite3, take the best backbones
#   1 GEN      RFD4 partial diffusion (motif pinned by the parent's own conditioning), 2 noise levels
#   2 FILTER   CPU: keep children at least as good as their parent on Maxi's core criteria
#              and base burial (hosts' Burial score), top M per parent
#   3 MPNN     LigandMPNN, same settings as the overnight run
#   4 BOLTZ    seed 0, 4 MPNN sequences per child, interleaved so every parent's best
#              children are refolded first; stops launching chunks at KEMP_BOLTZ_STOP
#   5 COMPARE  host-style scores of children vs their parent (composite3)
#   6 SEEDS    multi-seed Boltz for the children that beat their parent (consistency)
#
# Env (all optional except the two source runs):
#   KEMP_JOB61  multi-seed run of the pool        (has scores/ranking.csv)
#   KEMP_JOB60  overnight unsteered run           (has boltz_out, sequences.csv)
#   KEMP_NPAR=12 KEMP_REPS=8 KEMP_LEVELS="0.3 0.5" KEMP_M=8 KEMP_SEED=401
#   KEMP_BOLTZ_STOP / KEMP_MS_STOP  epoch seconds after which no new Boltz chunk starts
#   KEMP_RUN=<existing run dir to RESUME>

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh
J61=${KEMP_JOB61:?}; J60=${KEMP_JOB60:?}
NPAR=${KEMP_NPAR:-12}; REPS=${KEMP_REPS:-8}; LEVELS=${KEMP_LEVELS:-"0.3 0.5"}; M=${KEMP_M:-8}
SEED=${KEMP_SEED:-401}; BS=${KEMP_BS:-4}
BSTOP=${KEMP_BOLTZ_STOP:-0}; MSTOP=${KEMP_MS_STOP:-0}
RUN=${KEMP_RUN:-$KEMP_RUNS/$(date +%Y-%m-%d)_partial-diffusion_job${SLURM_JOB_ID}_${NPAR}parents_${REPS}reps}
RFD4=/mnt/work/repos/RFD4-Proteina
PY="$KEMP_ROOT/env/kemp-orch/bin/python"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"; mkdir -p "$TMPDIR" "$RUN"; trap 'rm -rf "$TMPDIR"' EXIT
LOG="$RUN/STAGES.log"
stage() { echo "$(date +%H:%M:%S) $*" | tee -a "$LOG"; }
t0=$(date +%s); T=$t0
lap() { local n=$(date +%s); stage "  -> $1 took $(( (n-T)/60 )) min"; T=$n; }
stage "=== job $SLURM_JOB_ID on $(hostname) | parents=$NPAR reps=$REPS levels=[$LEVELS] M=$M seed=$SEED"
echo "run dir: $RUN" | tee -a "$LOG"

# ------------------------------------------------------------ 0 PARENTS
stage "STAGE 0 PARENTS (CPU)"
if [ ! -s "$RUN/parents/parents.csv" ]; then
  mkdir -p "$RUN/parent_scoring"
  # single-seed scores of the overnight designs (job 60), on the same footing as the pool
  "$PY" -I "$KEMP_ROOT/scoring/pool_candidates.py" "$RUN/parent_scoring/overnight_pool" --overnight "$J60" >> "$LOG" 2>&1
  "$PY" -I - "$RUN/parent_scoring/overnight_pool/candidates.csv" <<'PYEOF'
import csv, sys
p = sys.argv[1]; rows = [r for r in csv.DictReader(open(p)) if r["source"] == "overnight"]
with open(p.replace("candidates.csv", "overnight_only.csv"), "w") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
PYEOF
  "$PY" -I "$KEMP_ROOT/scoring/host_score.py" "$RUN/parent_scoring/overnight_pool/overnight_only.csv" \
        "$J60/boltz_out" "$RUN/parent_scoring/overnight_scores" --top 5 >> "$LOG" 2>&1
  "$PY" -I "$KEMP_ROOT/partial/merge_rankings.py" "$RUN/parent_scoring/merged_ranking.csv" \
        "$J61/scores/ranking.csv" "$RUN/parent_scoring/overnight_scores/ranking.csv" | tee -a "$LOG"
  "$PY" -I "$KEMP_ROOT/partial/select_parents.py" "$RUN/parent_scoring/merged_ranking.csv" "$RUN/parents" \
        --n "$NPAR" --overnight "$J60" | tee -a "$LOG"
fi
[ -s "$RUN/parents/parents.csv" ] || { stage "FATAL: no parents"; exit 1; }
lap PARENTS

# -------------------------------------------------------------- 1 GEN
stage "STAGE 1 GEN (partial diffusion)"
NP=$(( $(wc -l < "$RUN/parents/parents.csv") - 1 ))
RAWARGS=""
for L in $LEVELS; do
  RAW="$RUN/raw_gen_f$L"; RAWARGS="$RAWARGS $RAW=$L"
  count() { find "$RAW" -name '*_generated.cif' 2>/dev/null | wc -l; }
  gen() {
    ( cd "$RFD4" && VIRTUAL_ENV="$RFD4/.venv-cu12" PATH="$RFD4/.venv-cu12/bin:$PATH" \
      "$RFD4/.venv-cu12/bin/rfd4-generate" inputs="$RUN/parents/partial_spec.yaml" \
        generation.args.partial_trajectory_fraction_steps="$L" \
        out_dir="$RAW" num_replicates="$REPS" batch_size="$1" compile_lengths=null \
        seed="$SEED" skip_existing=true '~metrics' > "$RUN/gen_f${L}_bs$1.log" 2>&1 ) ; return $?
  }
  if [ "$(count)" -lt $((NP*REPS)) ]; then gen "$BS"; fi
  if [ "$(count)" -lt $((NP*REPS)) ] && [ "$BS" -gt 2 ]; then
    stage "  level $L: only $(count)/$((NP*REPS)) after bs=$BS -> resuming with bs=2"; gen 2
  fi
  stage "  level $L: generated $(count)/$((NP*REPS)) children"
  if [ "$(count)" -eq 0 ]; then
    stage "FATAL: partial diffusion produced nothing at level $L; log tail:"; tail -15 "$RUN"/gen_f${L}_bs*.log | tee -a "$LOG"; exit 1
  fi
done
lap GEN

# ------------------------------------------------------------ 2 FILTER
stage "STAGE 2 FILTER (CPU)"
"$PY" -I "$KEMP_ROOT/partial/select_children.py" "$RUN/parents/parents.csv" "$RUN/children" $RAWARGS \
      --m "$M" --procs 4 2>&1 | tee "$RUN/filter.log" | tee -a "$LOG" | tail -30
[ -s "$RUN/children/selection.csv" ] && [ "$(wc -l < "$RUN/children/selection.csv")" -gt 1 ] \
  || { stage "FATAL: no child survived the filter"; exit 1; }
lap FILTER

# --------------------------------------------------------------- 3 MPNN
stage "STAGE 3 MPNN"
if [ ! -s "$RUN/sequences_all.csv" ]; then
  "$PY" -I "$KEMP_ROOT/overnight/mpnn_batch.py" "$RUN/children/selection.csv" "$RUN" 2>&1 | tail -6 | tee -a "$LOG"
  [ -s "$RUN/sequences.csv" ] || { stage "FATAL: no sequences"; exit 1; }
fi
"$PY" -I "$KEMP_ROOT/partial/order_sequences.py" "$RUN" | tee -a "$LOG"; lap MPNN

# -------------------------------------------------------------- 4 BOLTZ
stage "STAGE 4 BOLTZ (seed 0)"
"$PY" -I "$KEMP_ROOT/overnight/boltz_batch.py" "$RUN/sequences.csv" "$RUN" 25 | tee -a "$LOG"
BC="$RUN/boltz-cache"; mkdir -p "$BC"
for A in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do [ -e "$BC/$A" ] || ln -s "$BOLTZ_ASSETS/$A" "$BC/$A"; done
[ -e "$BC/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BC/mols"
for D in "$RUN"/boltz_in/chunk_*; do
  C=$(basename "$D"); [ -d "$D" ] || continue
  if [ "$BSTOP" -gt 0 ] && [ "$(date +%s)" -ge "$BSTOP" ]; then stage "  BOLTZ_STOP reached before $C; later chunks (lower-ranked children) not refolded"; break; fi
  n_in=$(ls "$D"/*.yaml | wc -l); n_out=$(find "$RUN/boltz_out/$C" -name '*_model_0.cif' 2>/dev/null | wc -l)
  [ "$n_out" -ge "$n_in" ] && { stage "  $C already complete"; continue; }
  "$ENV_ROOT/env-boltz2/bin/boltz" predict "$D" --model boltz2 --cache "$BC" --out_dir "$RUN/boltz_out/$C" \
      --accelerator gpu --devices 1 --no_kernels --diffusion_samples 1 --recycling_steps 3 \
      --output_format mmcif --num_workers 2 --preprocessing-threads 2 --seed 0 >> "$RUN/boltz.log" 2>&1
  stage "  $C rc=$? ($(find "$RUN/boltz_out/$C" -name '*_model_0.cif' | wc -l)/$n_in predictions)"
done
lap BOLTZ

# ------------------------------------------------------------- 5 COMPARE
stage "STAGE 5 COMPARE (CPU)"
"$PY" -I "$KEMP_ROOT/scoring/host_score.py" "$RUN/child_candidates.csv" "$RUN/boltz_out" "$RUN/scores_seed0" --top 15 2>&1 \
      | tee "$RUN/score_seed0.log" | tee -a "$LOG" | tail -22
"$PY" -I "$KEMP_ROOT/partial/compare_children.py" "$RUN" 2>&1 | tee -a "$LOG"
lap COMPARE

# --------------------------------------------------------------- 6 SEEDS
stage "STAGE 6 MULTI-SEED for children that beat their parent"
if [ "$(wc -l < "$RUN/multiseed_candidates.csv")" -gt 1 ]; then
  KEMP_CANDS="$RUN/multiseed_candidates.csv" KEMP_OUT="$RUN/multiseed" KEMP_SEEDS="0 1 2" KEMP_CHUNK=12 \
  KEMP_STOP_AT="$MSTOP" bash "$KEMP_ROOT/scoring/boltz_seeds.sbatch" 2>&1 | tail -30 | tee -a "$LOG"
else
  stage "  no child beat its parent; nothing to refold with more seeds"
fi
lap SEEDS
stage "=== ALL DONE in $(( ($(date +%s)-t0)/60 )) min ==="
