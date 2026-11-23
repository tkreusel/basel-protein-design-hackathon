#!/bin/bash
# Stage P: PLACER for a run whose stage B finished on a GPU that cannot run it.   usage: stage_P.sh <motif> <arm> <batch> [--extended]
set -uo pipefail
source "$(dirname "$0")/lib.sh"
MOTIF=$1 ARM=$2 B=$3; X=""; [ "${4:-}" = "--extended" ] && X="X"
RUN=$(run_dir "$MOTIF" "$ARM" "$B"); E=$RUN/eval${X:+_x}
[ -f "$RUN/STAGE_P${X}_DONE" ] && exit 0
acquire_lock "$RUN/.lockP$X" || exit 0
trap 'rm -rf "$RUN/.lockP$X"' EXIT
cd "$E" || exit 1; export KEMP_WORK=$E/work; source env.sh; export OMP_NUM_THREADS=2
i=0
while read -r L; do i=$((i+1))
  $PY_PLACER -m maxi_eval.run_placer --inputs work/inputs --out work/placer -n "${NSAMPLES:-20}" --modes holo,apo --max-replicates 1 --threads 2 --entries "$L" > "work/placer_$i.log" 2>&1 &
done < <($PY_TMOL -m maxi_eval.shard_entries -k "${NPL:-3}" 2>/dev/null)
wait
$PY_TMOL -m maxi_eval.analyze_placer --placer work/placer --out work/placer_summary.csv --field-models 5 2>&1 | tail -2
$PY_TMOL -m maxi_eval.aggregate 2>&1 | tail -1
date > "$RUN/STAGE_P${X}_DONE"
