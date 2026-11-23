# Shared helpers. Source from every job script:  source "$(dirname "$0")/lib.sh"
# shellcheck shell=bash

KEMP_PIPELINE="${KEMP_PIPELINE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
export KEMP_PIPELINE
if [ ! -f "$KEMP_PIPELINE/config/site.env" ]; then
  echo "config/site.env not found; copy config/site.env.example and edit it" >&2
  return 1 2>/dev/null || exit 1
fi
# shellcheck disable=SC1091
source "$KEMP_PIPELINE/config/site.env"
export CLUSTER=local HF_HUB_OFFLINE=1 TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1
export KEMP_SUITE="${KEMP_SUITE:-$KEMP_PIPELINE/../06_evaluation_suite/suite}"
CAMP="${CAMP:?CAMP is not set}"
export CAMP
mkdir -p "$KEMP_SCRATCH/logs/$CAMP"

# python -m kemp.<module> ...   (stage scripts live in src/kemp)
kpy() { PYTHONPATH="$KEMP_PIPELINE/src${PYTHONPATH:+:$PYTHONPATH}" "$PY_RFD4" -s -m "kemp.$1" "${@:2}"; }

# run_dir <motif> <arm> <batch>  ->  $KEMP_SCRATCH/runs/<camp>/<motif>/<arm>/bNN
run_dir() { printf '%s/runs/%s/%s/%s/b%02d\n' "$KEMP_SCRATCH" "$CAMP" "$1" "$2" "$((10#$3))"; }
all_runs() { ls -d "$KEMP_SCRATCH/runs/$CAMP"/*/*/b* 2>/dev/null; }
# parse_run <run dir>  ->  sets m (motif), a (arm), b (batch number)
parse_run() { m=$(basename "$(dirname "$(dirname "$1")")"); a=$(basename "$(dirname "$1")"); b=$(basename "$1"); b=$((10#${b#b})); }

# boltz_run <input dir> <output dir> <diffusion samples> <seed>: one `boltz predict` over all YAMLs in the input directory (finished predictions are skipped)
boltz_run() {
  local n; n=$(ls "$1"/*.yaml 2>/dev/null | wc -l); [ "$n" = 0 ] && return 0
  "$PY_BOLTZ_BIN/boltz" predict "$1" --out_dir "$2" --accelerator gpu --devices 1 --diffusion_samples "$3" --recycling_steps 3 \
      --output_format mmcif --write_full_pae --no_kernels --num_workers 4 --seed "$4" --cache "$BOLTZ_CACHE" 2>&1 | grep -v -i -E "warn|it/s|^$" | tail -3
}

# acquire_lock <dir>: atomic mkdir lock. A lock whose Slurm job is no longer in the queue is taken over.
acquire_lock() {
  local d=$1 owner
  if ! mkdir "$d" 2>/dev/null; then
    owner=$(cat "$d/jobid" 2>/dev/null)
    if [ -n "$owner" ] && squeue -h -j "$owner" 2>/dev/null | grep -q .; then return 1; fi
    rm -rf "$d"; mkdir "$d" 2>/dev/null || return 1
  fi
  echo "${SLURM_JOB_ID:-manual}" > "$d/jobid"
}

# remaining_s: seconds left in this Slurm job (999999 outside Slurm)
remaining_s() {
  local l; l=$(squeue -h -j "${SLURM_JOB_ID:-0}" -o %L 2>/dev/null | head -1)
  [ -z "$l" ] && { echo 999999; return; }
  python3 - "$l" <<'PY'
import sys
s, d = sys.argv[1], 0
if '-' in s:
    d, s = s.split('-'); d = int(d)
p = [int(x) for x in s.split(':')]; p = [0] * (3 - len(p)) + p
print(d * 86400 + p[0] * 3600 + p[1] * 60 + p[2])
PY
}
# seconds a worker needs to start stage A of an arm (steered runs take about 4.5 h including Boltz)
est_A() { case "$1" in *_sheet) echo 5400;; unsteered) echo 8000;; *) echo 17000;; esac; }
log_stage() { echo "[$1 $(date +%H:%M:%S) +$(( $(date +%s) - ${T0:-$(date +%s)} ))s] ${*:2}"; }
