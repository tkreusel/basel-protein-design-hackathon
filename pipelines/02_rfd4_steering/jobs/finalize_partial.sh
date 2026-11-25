#!/bin/bash
# Deadline tool: end the generation of a still-running run now and keep the chunks already finished (each chunk = 25 replicates).
# usage: finalize_partial.sh <motif> <arm> <batch>
# Merges the finished chunk tables into gen/metrics.csv and kills that run's gen process; the stage-A task then continues with the gate,
# LigandMPNN and Boltz on the designs that exist. Writes gen/PARTIAL.
source "$(dirname "$0")/lib.sh"; RUN=$(run_dir "$1" "$2" "$3"); G=$RUN/gen
[ -d "$G/chunks" ] || { echo "no chunks in $G"; exit 1; }
[ -f "$G/metrics.csv" ] && { echo "already has metrics.csv"; exit 0; }
"$PY_RFD4" -s - <<PY
import glob, os, pandas as pd
done = [x for x in sorted(glob.glob("$G/chunks/c*/metrics.csv")) if os.path.exists(os.path.dirname(x) + "/DONE")]
assert done, "no finished chunk"
d = pd.concat([pd.read_csv(x) for x in done], ignore_index=True); d.to_csv("$G/metrics.csv", index=False)
open("$G/PARTIAL", "w").write(f"{len(done)} chunks, {len(d)} designs\n"); print("merged", len(done), "chunks,", len(d), "designs")
PY
for j in $(squeue -u "$USER" -h -t R -o "%A %j" | awk '$2=="kemp_node"{print $1}'); do
  srun --jobid="$j" --overlap -N1 -n1 bash -c "pkill -f -- '--out $G' && echo killed on job $j" 2>/dev/null
done; true
