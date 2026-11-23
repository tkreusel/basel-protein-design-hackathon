#!/usr/bin/env bash
# Prepare the RFD4-Proteina clone this pipeline expects at <pipeline>/rfd4:
#   1. clone RFD4-Proteina at the pinned commit
#   2. drop bionemo-ir from the extras (its wheel needs GLIBC 2.34, which some clusters lack)   rfd4_patch/pyproject.diff
#   3. copy the Kemp reward metrics and the metric registry (both are gitignored upstream, so a plain clone does not have them)
#   4. create the virtual environment (skipped with --no-venv)
#   5. write the experiment configs of the steered arms (src/kemp/make_configs.py)
# usage: install_rfd4.sh [--no-venv]        env: CUDA=cu13 (default) | cu12, RFD4_REPO (clone URL), UV_CACHE_DIR, PY (interpreter with PyYAML when using --no-venv)
set -euo pipefail
P="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="${RFD4_REPO:-https://github.com/RosettaCommons/RFD4-Proteina}"
COMMIT=6306a210dd17dcae8e4e0f381b5b74a719d9b30e         # 2026-11-18, the commit used for every design of the campaign
CUDA="${CUDA:-cu13}"

[ -d "$P/rfd4/.git" ] || git clone "$REPO" "$P/rfd4"
git -C "$P/rfd4" checkout -q "$COMMIT"
if git -C "$P/rfd4" apply --check "$P/rfd4_patch/pyproject.diff" 2>/dev/null; then git -C "$P/rfd4" apply "$P/rfd4_patch/pyproject.diff"; fi
mkdir -p "$P/rfd4/src/rfproteina/metrics/custom_metrics" "$P/rfd4/configs/metrics/custom_metrics"
cp "$P"/rfd4_patch/custom_metrics/*.py "$P/rfd4/src/rfproteina/metrics/custom_metrics/"
cp "$P"/rfd4_patch/configs/metrics/all_custom_metrics.yaml "$P/rfd4/configs/metrics/custom_metrics/"
if [ "${1:-}" != "--no-venv" ]; then
  cd "$P/rfd4"
  uv venv --python 3.12 .venv
  BC=$(.venv/bin/python -c "import tomllib; print(next(d for d in tomllib.load(open('pyproject.toml','rb'))['project']['dependencies'] if d.startswith('torch')))")
  uv pip install --python .venv/bin/python --build-constraint <(printf '%s\n' "$BC") --torch-backend "$CUDA" -e ".[all-$( [ "$CUDA" = cu12 ] && echo cuda12 || echo cuda13 )]"
  uv pip check --python .venv/bin/python
  PY="$P/rfd4/.venv/bin/python"
fi
# the experiment configs need PyYAML: use the new venv, or PY=<interpreter> with --no-venv
PYTHONPATH="$P/src" "${PY:-python3}" -m kemp.make_configs --out "$P/rfd4/configs/inference_experiment"
echo "RFD4 ready at $P/rfd4 (commit $COMMIT, $CUDA)"
