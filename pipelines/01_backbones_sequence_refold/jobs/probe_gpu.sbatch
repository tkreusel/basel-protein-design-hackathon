#!/usr/bin/env bash
#SBATCH --job-name=kemp-probe
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=00:15:00

# Milestone 0 probe. Answers three questions that gate the whole pipeline:
#   1. What GPU driver / CUDA version do the compute nodes actually have?
#   2. Do the compute nodes have outbound internet (plan 2.0: the 2am killer)?
#   3. Does the RFD4-Proteina venv (torch cu132) see the GPU, or only the
#      provided cu126 envs?
# Deliberately does not use `set -e`: every probe must run even if one fails.

set -o pipefail

echo "=== node: $(hostname) ==="
nvidia-smi
echo
echo "=== driver version ==="
modinfo nvidia 2>/dev/null | grep -i '^version'
echo
echo "=== compute-node internet egress ==="
timeout 20 curl -sS -o /dev/null -w "pypi http=%{http_code}\n" https://pypi.org/simple/tmol/ \
    || echo "NO EGRESS from compute node"
timeout 20 curl -sS -o /dev/null -w "rcsb http=%{http_code}\n" https://files.rcsb.org/download/7K4Z.pdb \
    || echo "NO EGRESS to RCSB from compute node"
echo
echo "=== provided env-rfd3 (torch cu126) ==="
/mnt/ref/software/envs/env-rfd3/bin/python -c "
import torch
ok = torch.cuda.is_available()
print('cuda_available', ok, 'build', torch.version.cuda, torch.cuda.get_device_name(0) if ok else '')
" 2>&1 | tail -3
echo
echo "=== RFD4-Proteina venv (torch cu132) ==="
/mnt/work/repos/RFD4-Proteina/.venv/bin/python -c "
import torch
ok = torch.cuda.is_available()
print('cuda_available', ok, 'build', torch.version.cuda, torch.cuda.get_device_name(0) if ok else '')
" 2>&1 | tail -3
echo
echo "=== done ==="
