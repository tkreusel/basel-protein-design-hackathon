#!/usr/bin/env bash
# Evaluation suite environment. Source before running anything; PIPE_ROOT is the folder that contains this file.
# One interpreter per tool (they cannot be merged). Override any of them by exporting it before sourcing,
# e.g. from the campaign's config/site.env (pipelines/02_rfd4_steering).
_HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
export PIPE_ROOT="${PIPE_ROOT:-$_HERE}"
export PYTHONPATH="$PIPE_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

# tmol and analysis: python 3.12, tmol 0.2.0, torch cu126, biotite, rdkit, scipy, pandas. Also runs aggregate and validate; CPU or GPU.
export PY_TMOL="${PY_TMOL:?set PY_TMOL (the RFD4 venv works)}"
# PLACER: python 3.10, torch 2.3.1 cu12, DGL 2.4.0. GPU stage.
export PY_PLACER="${PY_PLACER:?set PY_PLACER}"
export PLACER_ROOT="${PLACER_ROOT:?set PLACER_ROOT}"
# ipSAE2 command line tool (CPU)
export IPSAE_BIN="${IPSAE_BIN:?set IPSAE_BIN}"
# PyRosetta 2026.47 (CPU only; academic licence)
export PY_ROSETTA="${PY_ROSETTA:?set PY_ROSETTA}"

# caches inside the work folder, not $HOME
export XDG_CACHE_HOME="${KEMP_WORK:-$PIPE_ROOT/work}/.cache"
export CLUSTER=local
export TORCH_EXTENSIONS_DIR="$XDG_CACHE_HOME/torch-extensions"
export MPLCONFIGDIR="$XDG_CACHE_HOME/matplotlib"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
