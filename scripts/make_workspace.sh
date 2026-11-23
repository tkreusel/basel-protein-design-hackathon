#!/usr/bin/env bash
# Build a working directory with the folder layout the scripts were written for
# ($KEMP_ROOT/{src,configs,inputs,ligands,jobs,scoring,overnight,partial,eligibility,env,runs,logs,cache}),
# using symlinks into this repository. Runs, logs and caches are created as real (empty) directories.
#
#   scripts/make_workspace.sh /path/to/workdir
#   cd /path/to/workdir && source env/paths.sh
#
# Then follow docs/03_REPRODUCING.md. env/paths.sh sources the organisers' reference-paths.sh (cluster specific);
# edit it for your site. The kemp-orch venv is NOT created here: see docs/02_ENVIRONMENTS_AND_CLUSTER.md.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
W="${1:?usage: make_workspace.sh <workdir>}"
mkdir -p "$W"/{runs,logs,cache,env}
ln -sfn "$REPO/core/src"                                   "$W/src"
ln -sfn "$REPO/core/configs"                               "$W/configs"
ln -sfn "$REPO/core/tests"                                 "$W/tests"
ln -sfn "$REPO/core/ligands"                               "$W/ligands"
ln -sfn "$REPO/pipelines/01_backbones_sequence_refold/inputs" "$W/inputs"
ln -sfn "$REPO/pipelines/01_backbones_sequence_refold/jobs"   "$W/jobs"
ln -sfn "$REPO/pipelines/04_overnight_unsteered_control"   "$W/overnight"
ln -sfn "$REPO/pipelines/05_partial_diffusion"             "$W/partial"
ln -sfn "$REPO/pipelines/07_host_style_scoring"            "$W/scoring"
ln -sfn "$REPO/pipelines/08_eligibility"                   "$W/eligibility"
ln -sfn "$REPO/pipelines/11_selection_and_submission"      "$W/submission"
ln -sfn "$REPO/pipelines/09_bioemu"                        "$W/bioemu"
ln -sfn "$REPO/pipelines/06_evaluation_suite"              "$W/maxi_eval"
ln -sfn "$REPO/pipelines/02_rfd4_steering"                 "$W/steering"
cp -n "$REPO/core/env/paths.sh" "$W/env/paths.sh"
echo "workspace ready: $W"
echo "next: scripts/fetch_benchmark.sh $W ; create the kemp-orch venv at $W/env/kemp-orch (docs/02) ; cd $W && source env/paths.sh && export KEMP_ROOT=$W"
