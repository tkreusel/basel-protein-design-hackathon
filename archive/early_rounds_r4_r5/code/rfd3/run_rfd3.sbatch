#!/bin/bash
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -J 20261121_rfd3_kemp_test
#SBATCH -o /home/hd/hd_hd/hd_la299/HACKATHON_ETH/logs/%x_%j.out
source activate fkguid_boltz
cd ~/HACKATHON_ETH/rfd3
echo "start $(date +%s)"; hostname
rfd3 design out_dir=out/20261121_test inputs=inputs/kemp_rfd3.json ckpt_path=/home/hd/hd_hd/hd_la299/rfd3_latest.ckpt diffusion_batch_size=2 n_batches=1 prevalidate_inputs=True
echo "end $(date +%s)"
