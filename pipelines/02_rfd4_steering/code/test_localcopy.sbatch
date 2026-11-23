#!/bin/bash
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:15:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH -J 20261121_rfd4_localcopy
#SBATCH -o /home/hd/hd_hd/hd_la299/HACKATHON_ETH/logs/%x_%j.out
echo TMPDIR=$TMPDIR; df -h /tmp ${TMPDIR:-/tmp} | tail -2
L=${TMPDIR:-/tmp}/rfd4env; mkdir -p $L
s=$(date +%s); tar -C ~/HACKATHON_ETH/rfd4 -cf - .venv | tar -C $L -xf -; echo copy $(( $(date +%s)-s ))s; du -sh $L
s=$(date +%s); $L/.venv/bin/python -c "import torch,rfproteina,boto3,tmol;print(torch.__version__,torch.cuda.is_available())" 2>&1 | tail -1; echo import $(( $(date +%s)-s ))s
