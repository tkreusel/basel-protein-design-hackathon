#!/bin/bash
# Progress of the campaign: queue, stage flags, steered generation progress, tracebacks in worker and stage-B logs.
source "$(dirname "$0")/lib.sh"
R=$KEMP_SCRATCH/runs/$CAMP; L=$KEMP_SCRATCH/logs/$CAMP
echo "== $(date +%H:%M)  jobs: $(squeue -u "$USER" -r -h -o '%j %T' | grep -E 'kemp' | sort | uniq -c | tr '\n' ';')"
echo "== stage flags A / B / P done: $(ls "$R"/*/*/b*/STAGE_A_DONE* 2>/dev/null | wc -l) / $(ls "$R"/*/*/b*/STAGE_B_DONE 2>/dev/null | wc -l) / $(ls "$R"/*/*/b*/STAGE_P_DONE 2>/dev/null | wc -l) of $(ls -d "$R"/*/*/b* 2>/dev/null | wc -l) run directories"
echo "== steered generation (designs / 300):"
for d in "$R"/*/steer_tmol*/b0*; do [ -d "$d/gen/chunks" ] && echo "  ${d#$R/} $(ls -d "$d"/gen/chunks/c*/inference_samples/* 2>/dev/null | wc -l)/300 $([ -f "$d/STAGE_A_DONE" ] && echo A_DONE)"; done | head -30
echo "== tracebacks:"; grep -l "Traceback" "$L"/W_*.log "$R"/*/*/b*/eval/work/*.log 2>/dev/null | head -5
