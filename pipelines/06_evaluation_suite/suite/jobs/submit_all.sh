#!/usr/bin/env bash
# DRY RUN by default: prints the sbatch chain. Pass --submit to actually submit (the user has asked NOT to submit yet).
# Team budget is 8 logical CPUs of GPU worker in total -> the three GPU jobs are chained, not parallel.
cd "$(dirname "$0")/.." && source env.sh
SUBMIT=0; [ "$1" = "--submit" ] && SUBMIT=1
run() { if [ $SUBMIT = 1 ]; then "$@"; else echo "[dry-run] $*"; fi; }
run bash jobs/00_prepare_and_ipsae.sh
J1=$( [ $SUBMIT = 1 ] && sbatch --parsable --chdir=$PWD --output=$PWD/logs/placer-%j.log jobs/10_placer.sbatch || echo "<jobid-placer>" )
echo "placer job: $J1"
J2=$( [ $SUBMIT = 1 ] && sbatch --parsable --dependency=afterok:$J1 --chdir=$PWD --output=$PWD/logs/tmol-%j.log jobs/20_tmol.sbatch || echo "<jobid-tmol>" )
echo "tmol job: $J2 (after placer)"
run sbatch --dependency=afterok:$J2 --chdir=$PWD --output=$PWD/logs/placer-an-%j.log jobs/30_placer_analysis.sbatch
run sbatch --dependency=afterok:$J2 --chdir=$PWD --output=$PWD/logs/pyrosetta-%j.log jobs/40_pyrosetta.sbatch   # CPU-only, no GPU reserved
echo "then (CPU): $PY_TMOL -m maxi_eval.aggregate"
