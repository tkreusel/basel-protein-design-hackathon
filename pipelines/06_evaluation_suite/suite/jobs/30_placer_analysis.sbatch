#!/usr/bin/env bash
#SBATCH --job-name=maxi-placer-an
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=01:00:00
# PLACER ensemble statistics + tmol field per model (top FIELD_MODELS by prmsd).
set -o pipefail
cd /mnt/work/projects/Maxi/pipeline && source env.sh
$PY_TMOL -m maxi_eval.analyze_placer --placer work/placer --out work/placer_summary.csv --field-models "${FIELD_MODELS:-20}"
