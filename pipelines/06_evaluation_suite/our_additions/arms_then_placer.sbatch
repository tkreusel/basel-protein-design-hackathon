#!/usr/bin/env bash
#SBATCH --job-name=kemp-arms-placer
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=04:00:00
# 1 Boltz-2 refolds of every finalist with the substrate (H5J) and product (6VP) ligand (same settings as our 6NT refolds:
#   msa empty, seed 0, 3 recycles, 1 sample) -> copied into boltz_root as <id>__H5J / <id>__6VP
# 2 wait for the TSA-only CPU shards (they read work/inputs), then re-run prepare + ipSAE2 so all three arms are registered
#   -> work/ARMS_PREPARED marker (the full tmol / PyRosetta runs on the head node start from it)
# 3 PLACER on all arms (holo x3 + apo), 20 samples, <= 3 replicates; resumable
set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh
cd /mnt/work/projects/tobi_v1/maxi_eval && source env.sh
LOG=logs/arms_then_placer.log; say() { echo "$(date +%H:%M:%S) $*" | tee -a $LOG; }
say "job $SLURM_JOB_ID on $(hostname)"
BC="$PWD/boltz-cache"; mkdir -p "$BC"
for A in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do [ -e "$BC/$A" ] || ln -s "$BOLTZ_ASSETS/$A" "$BC/$A"; done
[ -e "$BC/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BC/mols"
n_in=$(ls boltz_arms_in/*.yaml | wc -l); n_out=$(find boltz_arms_out -name '*_model_0.cif' 2>/dev/null | wc -l)
if [ "$n_out" -lt "$n_in" ]; then
  "$ENV_ROOT/env-boltz2/bin/boltz" predict boltz_arms_in --model boltz2 --cache "$BC" --out_dir boltz_arms_out \
      --accelerator gpu --devices 1 --no_kernels --diffusion_samples 1 --recycling_steps 3 --output_format mmcif \
      --num_workers 2 --preprocessing-threads 2 --seed 0 >> logs/boltz_arms.log 2>&1
fi
say "arm refolds: $(find boltz_arms_out -name '*_model_0.cif' | wc -l)/$n_in"
for d in $(find boltz_arms_out -type d -path '*predictions/*'); do
  id=$(basename "$d"); mkdir -p "boltz_root/$id"; cp -f "$d"/* "boltz_root/$id/"
done
for i in $(seq 1 40); do pgrep -f "maxi_eval.run_(tmol_entry|pyrosetta)" >/dev/null || break; sleep 30; done
say "CPU shards finished or waited 20 min; running prepare + ipSAE2"
BOLTZ_ROOT="$PWD/boltz_root" bash jobs/00_prepare_and_ipsae.sh > logs/prepare_ipsae_arms.log 2>&1; say "prepare rc=$?"
date > work/ARMS_PREPARED
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | tee -a $LOG
$PY_PLACER -m maxi_eval.run_placer --inputs work/inputs --out work/placer -n 20 --modes holo,apo --max-replicates 3 --threads 2 >> logs/placer.log 2>&1
say "PLACER rc=$?"; date > work/PLACER_DONE
