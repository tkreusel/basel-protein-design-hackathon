#!/bin/bash
# usage: sbatch -J name --export=ALL,N_PART=12 run_rfd3_arms.sbatch <key> "arm1 arm2 ..."
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH -o /home/hd/hd_hd/hd_la299/HACKATHON_ETH/logs/%x_%j.out
source activate fkguid_boltz
echo start $(date +%s)
python ~/HACKATHON_ETH/rfd3/run_arms_rfd3.py $1 $2
echo end $(date +%s)
