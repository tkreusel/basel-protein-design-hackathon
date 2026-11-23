#!/bin/bash
# Write the task list of a campaign and (optionally) submit whole-node workers for it. DRY RUN unless --submit is given.
# usage: submit.sh --motifs "hgnt hgwt" --batches "1 2 3" [--arms "unsteered steer_tmol steer_tmol_ch18"] [--sheet 1] [--nodes 3] [--extended 1] [--submit]
# One task = (motif, arm, batch) = 300 designs (steered: 100 replicates x 3 beam members; unsteered: 300 replicates).
# --sheet 1 adds a small beta-sheet-conditioned diversity run per (motif, arm, batch) with batch id 50+batch (steered 10 replicates, unsteered 30 designs).
source "$(dirname "$0")/lib.sh"
MOT="hgnt hgwt"; BAT="1"; ARMS="unsteered steer_tmol steer_tmol_ch18"; SHEET=0; NODES=1; EXT=0; SUB=0
while [ $# -gt 0 ]; do case $1 in
  --motifs) MOT=$2; shift 2;; --batches) BAT=$2; shift 2;; --arms) ARMS=$2; shift 2;; --sheet) SHEET=$2; shift 2;;
  --nodes) NODES=$2; shift 2;; --extended) EXT=$2; shift 2;; --submit) SUB=1; shift;; *) echo "unknown option $1"; exit 2;; esac; done
T=$KEMP_SCRATCH/runs/$CAMP/tasks_$(date +%Y%m%d_%H%M%S).txt; mkdir -p "$(dirname "$T")"; : > "$T"
for b in $BAT; do for m in $MOT; do for a in $ARMS; do
  echo "$m $a $b" >> "$T"
  if [ "$SHEET" = 1 ]; then n=10; [ "$a" = unsteered ] && n=30; echo "$m ${a}_sheet $((50+b)) $n" >> "$T"; fi
done; done; done
echo "campaign $CAMP: $(wc -l < "$T") tasks in $T"; head -40 "$T"
[ $SUB = 0 ] && { echo "DRY RUN, nothing submitted."; exit 0; }
for i in $(seq 1 "$NODES"); do
  sbatch --parsable -J kemp_node -p "$NODE_PARTITION" --gres=gpu:"$NODE_GPUS" -t "$NODE_TIME" -o "$KEMP_SCRATCH/logs/$CAMP/node_%j.log" \
         --export=ALL,TASKS="$T",WITH_EXTENDED="$([ "$EXT" = 1 ] && echo 1)" "$KEMP_PIPELINE/jobs/node_worker.sbatch"
done
