#!/bin/bash
# Checks of the shell helpers in jobs/lib.sh (no GPU or Slurm needed).   usage: bash tests/test_jobs.sh
set -u
P="$(cd "$(dirname "$0")/.." && pwd)"; T=$(mktemp -d); trap 'rm -rf "$T"' EXIT
cp -r "$P/config" "$T/config"; cp -r "$P/jobs" "$T/jobs"
cat > "$T/config/site.env" <<ENV
export CAMP=testcamp KEMP_SCRATCH=$T/scratch KEMP_RESULTS=$T/results PY_RFD4=python3 PY_BOLTZ_BIN=/none BOLTZ_CACHE=/none PLACER_GPU_PATTERN=H200
ENV
export KEMP_PIPELINE=$T; source "$T/jobs/lib.sh" || exit 1
fail=0; check() { [ "$2" = "$3" ] || { echo "FAIL $1: got '$2', expected '$3'"; fail=1; }; }

check run_dir "$(run_dir hgnt steer_tmol 3)" "$T/scratch/runs/testcamp/hgnt/steer_tmol/b03"
check run_dir_sheet "$(run_dir hgwt unsteered_sheet 51)" "$T/scratch/runs/testcamp/hgwt/unsteered_sheet/b51"
parse_run "$T/scratch/runs/testcamp/hgwt/steer_tmol_ch18/b07"; check parse_run "$m $a $b" "hgwt steer_tmol_ch18 7"
check est_A_steered "$(est_A steer_tmol_ch18)" 17000; check est_A_unsteered "$(est_A unsteered)" 8000; check est_A_sheet "$(est_A unsteered_sheet)" 5400
check remaining_outside_slurm "$(remaining_s)" 999999

mkdir -p "$T/lock"
acquire_lock "$T/lock/a"; check lock_first $? 0
acquire_lock "$T/lock/a"; first_retry=$?                      # owner is "manual", not in the queue: a stale lock is taken over
check lock_stale_takeover $first_retry 0

[ "$fail" = 0 ] && echo "jobs/lib.sh: all checks passed"; exit $fail
