#!/usr/bin/env bash
# Maxi Kemp-eliminase evaluation pipeline -- environment. Source before running anything.
# Three environments, one per tool; stages shell out between them (never mix them).
export MAXI_ROOT="/mnt/work/projects/Maxi"
export PIPE_ROOT="$MAXI_ROOT/pipeline"
export PYTHONPATH="$PIPE_ROOT/src${PYTHONPATH:+:$PYTHONPATH}"

# tmol + analysis (python 3.12, tmol 0.2.0, torch cu126, biotite, rdkit, scipy, pandas). GPU for relax.
export PY_TMOL="/mnt/work/repos/RFD4-Proteina/.venv-cu12/bin/python"
# PLACER (python 3.10, torch 2.3.1 cu12, dgl 2.4.0). GPU.
export PY_PLACER="/home/ec2-user/.conda/envs/placer_env/bin/python"
export PLACER_ROOT="/mnt/work/repos/placer"
# ipSAE2 (CPU only)
export IPSAE_BIN="/home/ec2-user/.conda/envs/ipsae2/bin/ipsae"

# caches off $HOME / (small, shared)
export XDG_CACHE_HOME="$PIPE_ROOT/work/.cache"
export TORCH_EXTENSIONS_DIR="$XDG_CACHE_HOME/torch-extensions"
export MPLCONFIGDIR="$XDG_CACHE_HOME/matplotlib"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$XDG_CACHE_HOME" "$MPLCONFIGDIR"
export PY_ROSETTA="/mnt/work/envs/pyrosetta/bin/python"   # PyRosetta 2026.47, CPU only
