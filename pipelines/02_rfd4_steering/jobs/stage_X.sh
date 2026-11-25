#!/bin/bash
# Stage X (extended evaluation): the next-best sequences of a finished run (kemp.extended_set) get the full suite:
# 9x TSA + 10x substrate + 10x product with Boltz-2, then stage_B.sh --extended.     usage: stage_X.sh <motif> <arm> <batch>
set -uo pipefail
source "$(dirname "$0")/lib.sh"; D=$(dirname "$0")
MOTIF=$1 ARM=$2 B=$3; RUN=$(run_dir "$MOTIF" "$ARM" "$B"); T0=$(date +%s)
[ -f "$RUN/STAGE_A_DONE" ] && [ -f "$RUN/boltz1_eval.csv" ] || exit 0
[ -f "$RUN/STAGE_BX_DONE" ] && exit 0
acquire_lock "$RUN/.lockX" || exit 0
trap 'rm -rf "$RUN/.lockX"' EXIT
export OMP_NUM_THREADS=4
[ -f "$RUN/survivors_x.txt" ] || kpy extended_set "$RUN" ${K_EXTRA:+--k "$K_EXTRA"}
[ "$(grep -c . "$RUN/survivors_x.txt")" -gt 0 ] || { date > "$RUN/STAGE_BX_DONE"; exit 0; }
SEQ=$RUN/mpnn/sequences.csv
for spec in "6NT 6NT_more_x 9 21" "H5J H5J_x 10 22" "6VP 6VP_x 10 23"; do
  read -r lig name n seed <<< "$spec"
  [ -d "$RUN/boltz_in/$name" ] || kpy boltz_prep --seq "$SEQ" --ligand "$lig" --ids "$RUN/survivors_x.txt" --out "$RUN/boltz_in/$name" --chunk 100000 | tail -1
  log_stage stageX "Boltz $name x$n"; boltz_run "$RUN/boltz_in/$name/chunk_000" "$RUN/boltz_out/$name" "$n" "$seed"
done
bash "$D/stage_B.sh" "$MOTIF" "$ARM" "$B" --extended
