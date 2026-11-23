#!/bin/bash
# Pack the key results of every finished run (stage P done, or stage A without survivors) into one tar.gz per run under $KEMP_RESULTS/run_archives.
# Scratch is not backed up and home directories often cap the file count, so one archive per run is the durable form. Idempotent.
source "$(dirname "$0")/lib.sh"
A=$KEMP_RESULTS/run_archives; mkdir -p "$A"
for r in $(all_runs); do
  [ -f "$r/STAGE_P_DONE" ] || [ -f "$r/STAGE_A_DONE_NO_SURVIVORS" ] || continue
  n=$(echo "$r" | sed "s|$KEMP_SCRATCH/runs/$CAMP/||;s|/|__|g"); out=$A/$n.tar.gz
  [ -f "$out" ] && [ "$out" -nt "$r/STAGE_B_DONE" ] && continue
  T=$(mktemp -d "$KEMP_SCRATCH/arch_tmp.XXXX"); d=$T/$n; mkdir -p "$d/best_boltz" "$d/designs" "$d/eval/tmol" "$d/eval/pyrosetta"
  cp "$r/gen/metrics.csv" "$d/design_metrics.csv" 2>/dev/null
  cp "$r"/gate*.csv "$r/sel.csv" "$r/batch_summary.csv" "$r/survivors.txt" "$r"/boltz1_eval*.csv "$r/mpnn/sequences.csv" "$d/" 2>/dev/null
  cp "$r"/eval_forest/{ranking,boltz_model_summary,stage5_report}.csv "$d/" 2>/dev/null
  cp -r "$r/eval/work/report" "$r/eval/work/placer_summary.csv" "$r/eval/work/inputs/manifest.json" "$d/eval/" 2>/dev/null
  cp "$r"/eval/work/tmol/*.json "$d/eval/tmol/" 2>/dev/null; cp "$r"/eval/work/pyrosetta/*.json "$d/eval/pyrosetta/" 2>/dev/null
  for f in "$r"/eval_forest/*/*_model_0.cif; do [ -e "$f" ] && cp -L "$f" "$d/best_boltz/"; done 2>/dev/null
  for s in $(cut -d, -f1 "$r/sel.csv" | tail -n +2); do cp "$r/gen/designs/$s.cif" "$d/designs/" 2>/dev/null; done
  tar -C "$T" -czf "$out" "$n" && echo "archived $n ($(du -h "$out" | cut -f1))"; rm -rf "$T"
done
