#!/bin/bash
# usage: sbatch -J <name> -o <log> run_boltz_chunk.sbatch <chunk_dir_name e.g. chunk_00>
#SBATCH -A <SLURM_ACCOUNT>
#SBATCH -p devel
#SBATCH --gres=gpu:1
#SBATCH -t 00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
D=$HOME/HACKATHON_ETH/runs/20261121_r5_downstream; C=$1
L=/dev/shm/boltzenv_$SLURM_JOB_ID; mkdir -p $L; s=$(date +%s); tar -C $L -xf $HOME/HACKATHON_ETH/boltz_env.tar; echo "env copy $(( $(date +%s)-s ))s"
export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1 PYTHONNOUSERSITE=1
n=$(ls $D/boltz_in/$C/*.yaml | wc -l); echo "chunk $C: $n inputs"; s=$(date +%s)
$L/boltz/bin/python -m boltz.main predict $D/boltz_in/$C --out_dir $D/boltz_out/$C --accelerator gpu --devices 1 --diffusion_samples 1 --recycling_steps 3 \
   --output_format mmcif --no_kernels --num_workers 2 --seed 0 --cache $HOME/.boltz
echo "predict done $(( $(date +%s)-s ))s"; rm -rf $L
