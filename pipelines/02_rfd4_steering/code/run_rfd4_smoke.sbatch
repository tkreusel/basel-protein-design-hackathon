#!/bin/bash
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH -J 20261121_rfd4_smoke
#SBATCH -o /home/hd/hd_hd/hd_la299/HACKATHON_ETH/logs/%x_%j.out
export CLUSTER=local HF_HUB_OFFLINE=1 PROJECT_ROOT=$HOME/HACKATHON_ETH/rfd4 DATA_PATH=$HOME/HACKATHON_ETH/rfd4/data_path
cd $HOME/HACKATHON_ETH/rfd4; source .venv/bin/activate
echo start $(date +%s)
rfd4-generate inputs=../rfd4_runs/inputs/kemp_rfd4.yaml num_replicates=4 batch_size=4 seed=1 out_dir=../rfd4_runs/out/20261121_smoke '~metrics' 
echo end $(date +%s)
