#!/bin/bash
# One worker bound to one GPU. Loops over the available work in this priority order until the walltime is too short or nothing is left:
#   P  PLACER for finished stage-B runs (only on GPUs matching $PLACER_GPU_PATTERN)
#   B  evaluation suite for runs whose stage A is done (needs 1.5 h of walltime)
#   A  generation-to-Boltz tasks from the task list $TASKS ("motif arm batch [n_rep]" per line; needs est_A seconds)
#   X  extended evaluation of finished non-sheet runs (only with WITH_EXTENDED=1)
# usage: worker.sh <gpu index>      env: TASKS, WITH_EXTENDED, FAN (extra copies per GPU, default 2; GPUs are CPU-bound, 3 processes per GPU gave about 2x throughput)
D=$(dirname "$0"); source "$D/lib.sh"
G=$1; export CUDA_VISIBLE_DEVICES=$G
if [ -z "${WF:-}" ] && [ -n "${SLURM_JOB_ID:-}" ]; then
  for k in $(seq 1 "${FAN:-2}"); do sleep 25; WF=1 bash "$0" "$G" > "$KEMP_SCRATCH/logs/$CAMP/W_${SLURM_JOB_ID}_g${G}_f$k.log" 2>&1 & done
fi
GPUN=$(nvidia-smi --query-gpu=name --format=csv,noheader -i "$G" | head -1); echo "[worker $G] $GPUN"
CAN_PLACER=0; echo "$GPUN" | grep -q "$PLACER_GPU_PATTERN" && CAN_PLACER=1
idle=0; max_idle=$([ -n "${WITH_EXTENDED:-}" ] && echo 120 || echo 30)

while true; do
  rem=$(remaining_s); [ "$rem" -lt 900 ] && break; did=0
  if [ $CAN_PLACER = 1 ]; then
    for r in $(all_runs); do parse_run "$r"
      for X in "" X; do
        [ -f "$r/STAGE_B${X}_DONE" ] && [ ! -f "$r/STAGE_P${X}_DONE" ] && [ ! -d "$r/.lockP$X" ] || continue
        echo "[worker $G] P$X $m $a $b"; bash "$D/stage_P.sh" "$m" "$a" "$b" ${X:+--extended}; did=1; break 2
      done
    done
  fi
  if [ $did = 0 ] && [ "$rem" -ge 5400 ]; then
    for r in $(all_runs); do
      [ -f "$r/STAGE_A_DONE" ] && [ -s "$r/survivors.txt" ] && [ ! -f "$r/STAGE_B_DONE" ] && [ ! -d "$r/.lockB" ] || continue
      parse_run "$r"; echo "[worker $G] B $m $a $b"; bash "$D/stage_B.sh" "$m" "$a" "$b"; did=1; break
    done
  fi
  if [ $did = 0 ] && [ -n "${TASKS:-}" ]; then
    while read -r m a b rest; do [ -z "$m" ] && continue
      RUN=$(run_dir "$m" "$a" "$b")
      ls "$RUN"/STAGE_A_DONE* >/dev/null 2>&1 && continue; [ -d "$RUN/.lockA" ] && continue
      [ "$(remaining_s)" -lt "$(est_A "$a")" ] && continue
      echo "[worker $G] A $m $a $b $rest (remaining $(remaining_s)s)"; bash "$D/stage_A.sh" "$m" "$a" "$b" $rest; did=1; break
    done < "$TASKS"
  fi
  if [ $did = 0 ] && [ -n "${WITH_EXTENDED:-}" ] && [ "$rem" -ge 5400 ]; then
    for r in $(all_runs | grep -v sheet); do
      [ -f "$r/STAGE_B_DONE" ] && [ ! -f "$r/STAGE_BX_DONE" ] && [ ! -d "$r/.lockX" ] && [ ! -d "$r/.lockBX" ] || continue
      parse_run "$r"; echo "[worker $G] X $m $a $b"; bash "$D/stage_X.sh" "$m" "$a" "$b"; did=1; break
    done
  fi
  if [ $did = 0 ]; then
    idle=$((idle+1)); [ "$idle" -ge "$max_idle" ] && { echo "[worker $G] nothing to do, exiting"; break; }; sleep 120
  else idle=0; fi
done
echo "[worker $G] exit"
