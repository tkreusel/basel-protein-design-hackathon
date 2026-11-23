#!/usr/bin/env bash
# Poll every 2 min until NEW complete Boltz predictions exist under ROOT and no `boltz` Slurm job is queued/running; exit 0 = ready.
# usage: jobs/wait_for_boltz.sh <root> <baseline_complete_models> [max_hours]
ROOT="$1"; BASE="$2"; MAXH="${3:-4}"; START=$(date +%s); PREV=-1; STABLE=0
complete_count() { n=0; while IFS= read -r c; do d=$(dirname "$c"); b=$(basename "$c" .cif); b=${b%_model_*}; k=$(basename "$c" .cif); k=${k##*_model_}
  [ -f "$d/pae_${b}_model_${k}.npz" ] && [ -f "$d/confidence_${b}_model_${k}.json" ] && n=$((n+1)); done < <(find "$ROOT" -name "*_model_*.cif" 2>/dev/null); echo $n; }
while true; do
  TOTAL=$(find "$ROOT" -name "*_model_*.cif" 2>/dev/null | wc -l); DONE=$(complete_count)
  JOBS=$(squeue -h -o "%j|%T" 2>/dev/null | grep -i boltz | tr '\n' ' ')
  echo "$(date +%H:%M:%S) cif=$TOTAL complete=$DONE baseline=$BASE boltz_jobs=[${JOBS:-none}]"
  if [ "$DONE" -gt "$BASE" ] && [ "$DONE" = "$TOTAL" ] && [ -z "$JOBS" ]; then
    if [ "$DONE" = "$PREV" ]; then STABLE=$((STABLE+1)); else STABLE=0; fi
    [ "$STABLE" -ge 1 ] && { echo "READY: $DONE complete models (was $BASE), no boltz job, count stable"; exit 0; }
  else STABLE=0; fi
  PREV=$DONE
  [ $(( $(date +%s) - START )) -gt $(( MAXH * 3600 )) ] && { echo "TIMEOUT after ${MAXH}h"; exit 3; }
  sleep 120
done
