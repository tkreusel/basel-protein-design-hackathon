#!/bin/bash
# Stage A (1 GPU): RFD4 generation -> design gate -> LigandMPNN -> Boltz-2 1x TSA -> Boltz gate -> survivors get 9x more TSA, 10x substrate, 10x product
# -> curated Boltz forest and per-sequence report.     usage: stage_A.sh <motif> <arm> <batch> [n_replicates]
# Every step is resumable: finished outputs are skipped. Environment knobs: BS (RFD4 batch size, 4), MIN_KEEP_DESIGNS, MIN_KEEP_SEQ.
set -uo pipefail
source "$(dirname "$0")/lib.sh"
MOTIF=$1 ARM=$2 B=$3 NREP=${4:-}
RUN=$(run_dir "$MOTIF" "$ARM" "$B"); mkdir -p "$RUN"
ls "$RUN"/STAGE_A_DONE* >/dev/null 2>&1 && { echo "already done: $RUN"; exit 0; }
acquire_lock "$RUN/.lockA" || { echo "running elsewhere"; exit 0; }
trap 'rm -rf "$RUN/.lockA"' EXIT
nvidia-smi --query-gpu=name --format=csv,noheader | tee -a "$RUN/gpu_A.txt"
export OMP_NUM_THREADS=4; T0=$(date +%s)
FILTER='warn|atomworks|environment variable|\(1\)|\(2\)|\(3\)'

# 1 generation (RFD4 reads its Hydra configs relative to the clone)
cd "$KEMP_PIPELINE/rfd4" || exit 1
log_stage stageA "run dir $RUN"
[ -f "$RUN/gen/metrics.csv" ] || kpy gen --motif "$MOTIF" --arm "$ARM" --batch "$B" ${NREP:+--n-rep "$NREP"} --batch-size "${BS:-4}" --out "$RUN/gen" > "$RUN/gen.log" 2>&1
grep -E "^\[t\]|Traceback|Error|tmol fail" "$RUN/gen.log" | cut -c1-200
[ -f "$RUN/gen/metrics.csv" ] || { log_stage stageA "generation produced no metrics"; exit 2; }

# 2 design gate with fallback
kpy filter_pre --metrics "$RUN/gen/metrics.csv" --designs "$RUN/gen/designs" --out "$RUN/gate.csv" --select "$RUN/sel.csv" --min-keep "${MIN_KEEP_DESIGNS:-30}" 2>&1 | grep -v -i warn

# 3 LigandMPNN, four shards sharing the GPU
if [ ! -f "$RUN/mpnn/sequences.csv" ]; then
  mkdir -p "$RUN/mpnn"
  for k in 0 1 2 3; do kpy mpnn --sel "$RUN/sel.csv" --out "$RUN/mpnn" --shard "$k/4" > "$RUN/mpnn/shard$k.log" 2>&1 & done; wait
  kpy mpnn --merge --out "$RUN/mpnn"
fi
SEQ=$RUN/mpnn/sequences.csv

# 4 Boltz-2, 1x TSA for every sequence
[ -d "$RUN/boltz_in/6NT" ] || kpy boltz_prep --seq "$SEQ" --ligand 6NT --out "$RUN/boltz_in/6NT" --chunk 100000 | tail -1
log_stage stageA "Boltz 1x 6NT"; boltz_run "$RUN/boltz_in/6NT/chunk_000" "$RUN/boltz_out/6NT" 1 0

# 5 refold gate with fallback
kpy boltz_gate --boltz-out "$RUN/boltz_out/6NT" --seq "$SEQ" --designs "$RUN/gate.csv" --out "$RUN/boltz1_eval.csv" --survivors "$RUN/survivors.txt" \
    --min-keep "${MIN_KEEP_SEQ:-10}" 2>&1 | grep -v -i -E "$FILTER"
if [ "$(grep -c . "$RUN/survivors.txt")" -eq 0 ]; then
  log_stage stageA "no survivors"; date > "$RUN/STAGE_A_DONE_NO_SURVIVORS"; kpy batch_summary "$RUN" | tail -2; exit 0
fi

# 6 survivors: 9 more TSA samples (seed 1), 10 substrate (seed 2), 10 product (seed 3)
for spec in "6NT 6NT_more 9 1" "H5J H5J 10 2" "6VP 6VP 10 3"; do
  read -r lig name n seed <<< "$spec"
  [ -d "$RUN/boltz_in/$name" ] || kpy boltz_prep --seq "$SEQ" --ligand "$lig" --ids "$RUN/survivors.txt" --out "$RUN/boltz_in/$name" --chunk 100000 | tail -1
  log_stage stageA "Boltz $name x$n"; boltz_run "$RUN/boltz_in/$name/chunk_000" "$RUN/boltz_out/$name" "$n" "$seed"
done

# 7 curate the forest (rank samples by pocket RMSD to the design) and write the per-sequence report
kpy select_models --boltz-dirs "$RUN/boltz_out/6NT" "$RUN/boltz_out/6NT_more" "$RUN/boltz_out/H5J" "$RUN/boltz_out/6VP" --ids "$RUN/survivors.txt" \
    --out "$RUN/eval_forest" --seq "$SEQ" --designs "$RUN/gate.csv" 2>&1 | grep -v -i -E "$FILTER" | tail -3
log_stage stageA "STAGE A DONE"; date > "$RUN/STAGE_A_DONE"
kpy batch_summary "$RUN" | tail -3
