#!/usr/bin/env bash
# FULL SUITE on a Boltz output folder with 6NT / H5J / 6VP per protein:
#   prepare -> ipSAE2 (all arms) -> PLACER holo on enzyme-TSA, enzyme-substrate AND enzyme-product complexes (+ apo) ->
#   tmol + PyRosetta (head-node CPU, sharded) -> PLACER analysis -> aggregate -> validate / discriminative power.
#
#   bash jobs/run_full_suite.sh --root /mnt/work/projects/Maxi/boltz_outputs            # DRY RUN: readiness, plan, cost (default)
#   bash jobs/run_full_suite.sh --root ... --submit                                   # actually runs (1 GPU job + head-node CPU)
# Options: --nsamples 20  --shards-tmol 2  --shards-ros 3  --ablate (tmol+PyRosetta ablations, slower)  --force (run despite incomplete arms)
# Rules: never submits without --submit; ONE GPU job (PLACER); check `squeue` for teammates first (8-CPU team budget).
set -o pipefail
cd "$(dirname "$0")/.." && source env.sh
ROOT=""; SUBMIT=0; NS=20; KT=2; KR=3; ABL=""; FORCE=0
while [ $# -gt 0 ]; do case "$1" in
  --root) ROOT="$2"; shift 2;; --submit) SUBMIT=1; shift;; --nsamples) NS="$2"; shift 2;;
  --shards-tmol) KT="$2"; shift 2;; --shards-ros) KR="$2"; shift 2;; --ablate) ABL="--ablate"; shift;; --force) FORCE=1; shift;;
  *) echo "unknown option $1"; exit 2;; esac; done
[ -z "$ROOT" ] && { echo "need --root <boltz output folder>"; exit 2; }
echo "=== 1. readiness"; $PY_TMOL -m maxi_eval.check_boltz_ready --root "$ROOT" --nsamples "$NS" --require-pae 2>&1 | grep -v -E "MIRROR|\(1\)|\(2\)|\(3\)|CUDA runtime"
READY=${PIPESTATUS[0]}
if [ "$READY" != 0 ] && [ $FORCE = 0 ]; then echo "NOT every protein has all three arms; fix or pass --force"; [ $SUBMIT = 1 ] && exit 1; fi
echo; echo "=== 2. queue (other users' jobs count against the 8-CPU team budget)"; squeue
echo; echo "=== 3. plan"
echo "  prepare + ipSAE2 (CPU, minutes)"; echo "  PLACER: sbatch jobs/10_placer_allarms.sbatch NSAMPLES=$NS  (1 GPU, all arms holo + one apo per protein)"
echo "  tmol x$KT and PyRosetta x$KR processes on the head node ${ABL:+(with ablations)}"; echo "  then analyze_placer, aggregate, validate, validate_dataset"
[ $SUBMIT = 0 ] && { echo; echo "DRY RUN - nothing started. Re-run with --submit after sign-off."; exit 0; }

echo; echo "=== 4. prepare"; BOLTZ_ROOT="$ROOT" bash jobs/00_prepare_and_ipsae.sh 2>&1 | grep -E "registered|processed|FAILED|warn|skip"
mkdir -p logs work/tmol work/pyrosetta
JID=$(sbatch --parsable --export=ALL,NSAMPLES=$NS,MAXREP=3 --chdir=$PWD --output=$PWD/logs/placer-full-%j.log jobs/10_placer_allarms.sbatch); echo "PLACER job: $JID"
i=0; while read -r L; do i=$((i+1)); nohup $PY_TMOL -m maxi_eval.run_tmol_entry --entries "$L" --out work/tmol --apo $ABL > work/full_tmol_$i.log 2>&1 &
done < <($PY_TMOL -m maxi_eval.shard_entries -k "$KT" 2>/dev/null)
i=0; while read -r L; do i=$((i+1)); nohup env PYTHONPATH=src $PY_ROSETTA -m maxi_eval.run_pyrosetta --entries "$L" --out work/pyrosetta --relax-repeats 1 $ABL --nstruct 1 > work/full_ros_$i.log 2>&1 &
done < <($PY_TMOL -m maxi_eval.shard_entries -k "$KR" 2>/dev/null)
# finisher: waits for the GPU job and all CPU processes, then analyses and writes work/report/FULL_SUITE_DONE
nohup bash -c "while squeue -j $JID -h 2>/dev/null | grep -q . || pgrep -f 'maxi_eval.run_(tmol_entry|pyrosetta)' >/dev/null; do sleep 30; done;
  cd $PWD && source env.sh && \$PY_TMOL -m maxi_eval.analyze_placer --placer work/placer --out work/placer_summary.csv --field-models 5 &&
  \$PY_TMOL -m maxi_eval.apo_holo && \$PY_TMOL -m maxi_eval.aggregate && \$PY_TMOL -m maxi_eval.validate && \$PY_TMOL -m maxi_eval.validate_dataset --group boltz --min-proteins 10;
  date > work/report/FULL_SUITE_DONE" > work/full_finish.log 2>&1 &
echo "started. progress: squeue; ls work/tmol work/pyrosetta | wc -l; tail work/full_*.log; done-marker: work/report/FULL_SUITE_DONE"
