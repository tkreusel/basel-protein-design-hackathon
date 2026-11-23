#!/usr/bin/env bash
# CPU-only, seconds-to-minutes: runs fine on the head node (no GPU, no sbatch needed).
# Standardises inputs (crystals + Boltz refolds) and scores ipSAE2. Safe to re-run.
set -o pipefail
cd "$(dirname "$0")/.." && source env.sh
$PY_TMOL -m maxi_eval.build_ligands                       # once; idempotent: 6NT / H5J / 6VP from the CCD (tmol params)
PYTHONPATH="$PIPE_ROOT/src" $PY_ROSETTA -m maxi_eval.build_rosetta_params   # Rosetta params for the same ligands
$PY_TMOL -m maxi_eval.prepare_inputs ${BOLTZ_ROOT:+--boltz-root "$BOLTZ_ROOT"} ${NAME_MAP:+--name-map "$NAME_MAP"}
$PY_TMOL -m maxi_eval.run_ipsae                           # needs Boltz refolds with PAE npz; no-op otherwise
