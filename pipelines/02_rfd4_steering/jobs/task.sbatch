#!/bin/bash
# One stage on one GPU:  sbatch -p "$SINGLE_PARTITION" --gres=gpu:1 -t 08:00:00 --cpus-per-task=12 --mem=150G jobs/task.sbatch A hgnt steer_tmol_ch18 1
# arguments: stage (A|B|P|X) motif arm batch [n_replicates | --extended]
D=${KEMP_PIPELINE:?export KEMP_PIPELINE=<this pipeline directory> before sbatch}/jobs; source "$D/lib.sh"
stage=$1; shift
bash "$D/stage_$stage.sh" "$@"
