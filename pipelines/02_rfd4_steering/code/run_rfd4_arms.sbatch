#!/bin/bash
# usage: sbatch run_rfd4_arms.sbatch "arm1 arm2" jobname
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH -o /home/hd/hd_hd/hd_la299/HACKATHON_ETH/logs/%x_%j.out
export HF_HUB_OFFLINE=1 PROJECT_ROOT=$HOME/HACKATHON_ETH/rfd4 DATA_PATH=$HOME/HACKATHON_ETH/rfd4/data_path CLUSTER=local
L=/tmp/rfd4env_$SLURM_JOB_ID; mkdir -p $L; s=$(date +%s)
tar -C $L -xf $HOME/HACKATHON_ETH/rfd4_venv.tar; echo "copy $(( $(date +%s)-s ))s"
cd $PROJECT_ROOT   # config dir is found from the package location
$L/.venv/bin/python ~/HACKATHON_ETH/rfd4_runs/run_arms.py $1
rm -rf $L
