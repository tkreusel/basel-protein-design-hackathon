#!/usr/bin/env bash
#SBATCH --job-name=bioemu_4Z08
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=01:00:00


INPUT="$1"
NUM_SAMPLES="$2"

source ~/.bashrc
echo "activated"
echo nvidia-smi
which python
conda activate bioemu
cd /mnt/work/repos/bioemu2
bioemu2 sample-equilibrium --system "$INPUT" --out /mnt/work/projects/Wasabi/outputs --num-samples "$NUM_SAMPLES"
