#!/bin/bash
# Smoke test of the whole chain at tiny scale (one steered channel-arm sheet batch with 2 replicates, relaxed gates); run before every large launch.
#   export KEMP_PIPELINE=$PWD; sbatch -p "$SINGLE_PARTITION" --gres=gpu:1 -t 01:00:00 --cpus-per-task=12 --mem=150G jobs/smoke.sbatch
D=${KEMP_PIPELINE:?export KEMP_PIPELINE=<this pipeline directory> before sbatch}/jobs; source "$D/lib.sh"
export CAMP=${CAMP}_smoke MIN_KEEP_DESIGNS=3 MIN_KEEP_SEQ=3 N_TMOL=3 N_ROS=3
mkdir -p "$KEMP_SCRATCH/logs/$CAMP"
bash "$D/stage_A.sh" hgnt steer_tmol_ch18_sheet 97 2 && bash "$D/stage_B.sh" hgnt steer_tmol_ch18_sheet 97
