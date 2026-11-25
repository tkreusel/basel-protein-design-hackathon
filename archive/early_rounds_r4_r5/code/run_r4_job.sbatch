#!/bin/bash
# usage: sbatch -J <name> -o <log> --export=ALL,N_REP=3 run_r4_job.sbatch <campaign> <arm> [<arm>...]
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
export HF_HUB_OFFLINE=1 PROJECT_ROOT=$HOME/HACKATHON_ETH/rfd4 DATA_PATH=$HOME/HACKATHON_ETH/rfd4/data_path CLUSTER=local
L=/tmp/rfd4env_$SLURM_JOB_ID; mkdir -p $L; s=$(date +%s)
tar -C $L -xf $HOME/HACKATHON_ETH/rfd4_venv.tar; echo "env copy $(( $(date +%s)-s ))s"
cd $PROJECT_ROOT
$L/.venv/bin/python $HOME/HACKATHON_ETH/run_r4.py "$@"
rm -rf $L
