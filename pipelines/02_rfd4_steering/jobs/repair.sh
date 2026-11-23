#!/bin/bash
# Re-run parts of stage B for finished runs, then re-aggregate. Idempotent.
#   repair.sh placer     PLACER analysis for runs whose placer_summary.csv is missing
#   repair.sh pyrosetta  PyRosetta for replicates whose JSON carries an interface error (an Sc failure used to discard dG_sep)
source "$(dirname "$0")/lib.sh"
case "${1:-}" in
placer)
  for r in $(all_runs); do
    [ -f "$r/STAGE_P_DONE" ] && [ -d "$r/eval/work/placer" ] && [ ! -f "$r/eval/work/placer_summary.csv" ] || continue
    ( cd "$r/eval" && export KEMP_WORK=$PWD/work && source env.sh && export OMP_NUM_THREADS=2
      $PY_TMOL -m maxi_eval.analyze_placer --placer work/placer --out work/placer_summary.csv --field-models 5 > work/analyze_placer.log 2>&1
      $PY_TMOL -m maxi_eval.aggregate > work/aggregate.log 2>&1 )
    echo "repaired ${r#$KEMP_SCRATCH/runs/$CAMP/} $([ -f "$r/eval/work/placer_summary.csv" ] && echo ok || echo STILL_MISSING)"
  done;;
pyrosetta)
  for r in $(all_runs); do
    [ -f "$r/STAGE_B_DONE" ] && [ -d "$r/eval/work/pyrosetta" ] || continue
    bad=$(grep -l '"ros.error\.\(native_sub\|native_prod\|[a-z_.]*interface\)"' "$r"/eval/work/pyrosetta/*.json 2>/dev/null); [ -z "$bad" ] && continue
    ( cd "$r/eval" && export KEMP_WORK=$PWD/work && source env.sh && export OMP_NUM_THREADS=2
      rm -f $bad; ENT=$($PY_TMOL -m maxi_eval.shard_entries -k 1 2>/dev/null)
      for i in 1 2 3 4; do sleep "$i"; R=""; [ $((i%2)) = 0 ] && R=--reverse
        env PYTHONPATH=src $PY_ROSETTA -m maxi_eval.run_pyrosetta --entries "$ENT" --out work/pyrosetta --relax-repeats 1 --nstruct 1 --max-replicates 1 $R > "work/ros_repair_$i.log" 2>&1 &
      done; wait
      $PY_TMOL -m maxi_eval.aggregate > work/aggregate.log 2>&1 )
    echo "repaired pyrosetta ${r#$KEMP_SCRATCH/runs/$CAMP/} ($(echo "$bad" | wc -w) replicates)"
  done;;
*) echo "usage: repair.sh placer|pyrosetta"; exit 2;;
esac
