#!/usr/bin/env bash
#SBATCH --job-name=kemp-arms2
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=01:00:00
# substrate (H5J) + product (6VP) Boltz-2 refolds for the extra pool designs (same settings as arms_then_placer.sbatch)
set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh
cd /mnt/work/projects/tobi_v1/maxi_eval && source env.sh
BC="$PWD/boltz-cache"
"$ENV_ROOT/env-boltz2/bin/boltz" predict boltz_arms_in2 --model boltz2 --cache "$BC" --out_dir boltz_arms_out2 \
    --accelerator gpu --devices 1 --no_kernels --diffusion_samples 1 --recycling_steps 3 --output_format mmcif \
    --num_workers 2 --preprocessing-threads 2 --seed 0 >> logs/boltz_arms2.log 2>&1
echo "arm refolds: $(find boltz_arms_out2 -name '*_model_0.cif' | wc -l)/20"
