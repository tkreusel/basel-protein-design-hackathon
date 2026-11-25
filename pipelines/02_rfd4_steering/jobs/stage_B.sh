#!/bin/bash
# Stage B: evaluation suite on the curated Boltz forest. ipSAE2 (all samples, three ligands) -> PLACER (holo x3 ligands + apo, 20 models, best sample)
# -> tmol protocol and PyRosetta (best TSA sample plus independent substrate/product refolds) -> aggregate.
# usage: stage_B.sh <motif> <arm> <batch> [--extended]       Resumable.
# --extended evaluates the extended set (survivors_x.txt; eval_x/, eval_forest_x/, flags STAGE_BX_DONE / STAGE_PX_DONE).
# PLACER runs inline on GPUs matching $PLACER_GPU_PATTERN; otherwise stage_P.sh finishes it later. Knobs: NSAMPLES (20), N_TMOL (3), N_ROS (6), NPL (3).
set -uo pipefail
source "$(dirname "$0")/lib.sh"
MOTIF=$1 ARM=$2 B=$3; X=""; [ "${4:-}" = "--extended" ] && X="X"
SUF=${X:+_x}; RUN=$(run_dir "$MOTIF" "$ARM" "$B")
E=$RUN/eval$SUF; FOREST=$RUN/eval_forest$SUF; IDS=$RUN/survivors$SUF.txt; DONE=$RUN/STAGE_B${X}_DONE; PDONE=$RUN/STAGE_P${X}_DONE
NS=${NSAMPLES:-20}; NT=${N_TMOL:-3}; NR=${N_ROS:-6}; NPL=${NPL:-3}; T0=$(date +%s)
[ -f "$IDS" ] || { echo "no $IDS"; exit 2; }
[ -f "$DONE" ] && { echo "already done"; exit 0; }
acquire_lock "$RUN/.lockB$X" || { echo "running elsewhere"; exit 0; }
trap 'rm -rf "$RUN/.lockB$X"' EXIT
nvidia-smi --query-gpu=name --format=csv,noheader | tee -a "$RUN/gpu_B$X.txt"

# the suite works in a run-local folder with symlinks to its code and a config pointing at the curated forest
mkdir -p "$E/work" "$E/config"; cd "$E" || exit 1
for x in src ligands jobs; do ln -sfn "$KEMP_SUITE/$x" "$x"; done; ln -sfn "$KEMP_SUITE/env.sh" env.sh
cp "$KEMP_SUITE"/config/*.yaml "$KEMP_SUITE"/config/*.csv config/; cp config/entries.design_run.yaml config/entries.yaml
sed -i "s|root: .*|root: $FOREST|" config/entries.yaml
export KEMP_WORK=$E/work; source env.sh; export OMP_NUM_THREADS=2 MKL_NUM_THREADS=2

if [ ! -f "$FOREST/stage5_report.csv" ]; then
  kpy select_models --boltz-dirs "$RUN/boltz_out/6NT" "$RUN/boltz_out/6NT_more$SUF" "$RUN/boltz_out/H5J$SUF" "$RUN/boltz_out/6VP$SUF" --ids "$IDS" \
      --out "$FOREST" --seq "$RUN/mpnn/sequences.csv" --designs "$RUN/gate.csv" | tail -1
fi
$PY_TMOL -m maxi_eval.check_boltz_ready --root "$FOREST" --require-pae 2>&1 | grep -E "proteins found|incomplete|lacking" | head -5
log_stage stageB "ligands, input preparation, ipSAE2"
$PY_TMOL -m maxi_eval.build_ligands >/dev/null 2>&1; PYTHONPATH=$PIPE_ROOT/src $PY_ROSETTA -m maxi_eval.build_rosetta_params >/dev/null 2>&1
$PY_TMOL -m maxi_eval.prepare_inputs 2>&1 | grep -E "registered|warn|skip" | head -5
$PY_TMOL -m maxi_eval.run_ipsae 2>&1 | tail -1

# PLACER in the background where the GPU allows it; tmol and PyRosetta shards run on the CPU meanwhile
PL_INLINE=0; nvidia-smi --query-gpu=name --format=csv,noheader | head -1 | grep -q "$PLACER_GPU_PATTERN" && PL_INLINE=1
if [ $PL_INLINE = 1 ]; then
  log_stage stageB "PLACER inline ($NPL shards)"
  ( i=0; while read -r L; do i=$((i+1))
      $PY_PLACER -m maxi_eval.run_placer --inputs work/inputs --out work/placer -n "$NS" --modes holo,apo --max-replicates 1 --threads 2 --entries "$L" > "work/placer_$i.log" 2>&1 &
    done < <($PY_TMOL -m maxi_eval.shard_entries -k "$NPL" 2>/dev/null); wait; echo placer_done > work/placer.log ) &
fi
ENT=$($PY_TMOL -m maxi_eval.shard_entries -k 1 2>/dev/null)
log_stage stageB "tmol x$NT, PyRosetta x$NR"; mkdir -p work/tmol work/pyrosetta
for i in $(seq 1 "$NT"); do sleep $((i*5)); R=""; [ $((i%2)) = 0 ] && R=--reverse
  $PY_TMOL -m maxi_eval.run_tmol_entry --entries "$ENT" --out work/tmol --apo --max-replicates 1 $R > "work/tmol_$i.log" 2>&1 &
done
for i in $(seq 1 "$NR"); do sleep $((i*3)); R=""; [ $((i%2)) = 0 ] && R=--reverse
  env PYTHONPATH=src $PY_ROSETTA -m maxi_eval.run_pyrosetta --entries "$ENT" --out work/pyrosetta --relax-repeats 1 --nstruct 1 --max-replicates 1 $R > "work/ros_$i.log" 2>&1 &
done
wait; log_stage stageB "analysis"
[ $PL_INLINE = 1 ] && $PY_TMOL -m maxi_eval.analyze_placer --placer work/placer --out work/placer_summary.csv --field-models 5 2>&1 | tail -2
$PY_TMOL -m maxi_eval.aggregate 2>&1 | tail -1
log_stage stageB "STAGE B$X DONE"; date > "$DONE"; [ $PL_INLINE = 1 ] && date > "$PDONE"
exit 0
