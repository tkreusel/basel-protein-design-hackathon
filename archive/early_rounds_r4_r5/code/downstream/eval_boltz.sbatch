#!/bin/bash
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
export HF_HUB_OFFLINE=1 PROJECT_ROOT=$HOME/HACKATHON_ETH/rfd4 DATA_PATH=$HOME/HACKATHON_ETH/rfd4/data_path CLUSTER=local
L=/dev/shm/rfd4env_$SLURM_JOB_ID; mkdir -p $L; tar -C $L -xf $HOME/HACKATHON_ETH/rfd4_venv.tar; cd $PROJECT_ROOT
$L/.venv/bin/python -I $HOME/HACKATHON_ETH/downstream/eval_boltz.py; rm -rf $L
