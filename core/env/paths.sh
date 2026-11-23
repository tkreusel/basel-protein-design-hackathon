#!/usr/bin/env bash
# Kemp eliminase pipeline -- central path and cache configuration.
# Source this at the top of EVERY job script (plan section 7).
#
# Design rule: inherit the organisers' known-good tool paths from
# /mnt/ref/hackathon/scripts/reference-paths.sh rather than restating them, so a
# reference update does not silently desync us. Everything we add is either a
# project path or a cache redirect.
#
# Plan section 2.1: no model weight cache may live in $HOME. $HOME is on / with
# ~69 GB free and is shared with the OS; /mnt/work has ~459 GB. Every cache below
# is therefore pinned under the project.

KEMP_ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)
export KEMP_ROOT

# --- provided tools, weights and environments (do not duplicate; inherit) -----
# Supplies: PDH_REF_ROOT PDH_TOOLKIT ENV_ROOT CONDA_ROOT LIGANDMPNN_ROOT
#           MMSEQS MMSEQS_DB RFD3_CHECKPOINT LIGANDMPNN_CHECKPOINT BOLTZ_ASSETS
source /mnt/ref/hackathon/scripts/reference-paths.sh

# --- our own tool locations ---------------------------------------------------
export RFD4_ROOT="/mnt/work/repos/RFD4-Proteina"
export RFD4_PYTHON="$RFD4_ROOT/.venv/bin/python"

# --- project layout ----------------------------------------------------------
export KEMP_CONFIGS="$KEMP_ROOT/configs"
export KEMP_LIGANDS="$KEMP_ROOT/ligands"
export KEMP_RUNS="$KEMP_ROOT/runs"
export KEMP_BENCHMARK="$KEMP_ROOT/benchmark"
export PYTHONPATH="$KEMP_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

# --- caches: all on /mnt/work, never $HOME ----------------------------------
KEMP_CACHE="$KEMP_ROOT/cache"
export XDG_CACHE_HOME="$KEMP_CACHE"
export HF_HOME="$KEMP_CACHE/huggingface"
export TORCH_HOME="$KEMP_CACHE/torch"
export TORCH_EXTENSIONS_DIR="$KEMP_CACHE/torch-extensions"
export TRITON_CACHE_DIR="$KEMP_CACHE/triton"
export CUDA_CACHE_PATH="$KEMP_CACHE/cuda"
export PIP_CACHE_DIR="$KEMP_CACHE/pip"
export CONDA_PKGS_DIRS="$KEMP_CACHE/conda"
export NUMBA_CACHE_DIR="$KEMP_CACHE/numba"
export MPLCONFIGDIR="$KEMP_CACHE/matplotlib"
export TMPDIR="$KEMP_ROOT/tmp"
export PYTHONDONTWRITEBYTECODE=1

mkdir -p "$HF_HOME" "$TORCH_HOME" "$TORCH_EXTENSIONS_DIR" "$TRITON_CACHE_DIR" \
         "$CUDA_CACHE_PATH" "$PIP_CACHE_DIR" "$CONDA_PKGS_DIRS" "$NUMBA_CACHE_DIR" \
         "$MPLCONFIGDIR" "$TMPDIR" "$KEMP_RUNS" "$KEMP_ROOT/logs"

# --- tool discovery ----------------------------------------------------------
export PATH="$(dirname -- "$MMSEQS"):$PATH"
source "$CONDA_ROOT/etc/profile.d/conda.sh"
