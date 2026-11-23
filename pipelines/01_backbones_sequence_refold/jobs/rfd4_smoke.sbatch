#!/usr/bin/env bash
#SBATCH --job-name=rfd4-smoke
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=00:50:00

# Does the rebuilt CUDA-12.6 RFD4 environment actually work on an L4?
#
# Three gates, in order. Later ones only mean something if the earlier pass:
#   1. torch sees the GPU (the original cu132 venv failed here)
#   2. RFD4's own preflight passes
#   3. one real generation completes end to end at batch_size=1
#
# No `set -e`: every gate must report, even after an earlier failure.

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

RFD4=/mnt/work/repos/RFD4-Proteina
PY="$RFD4/.venv-cu12/bin/python"
OUT="${KEMP_RUNS}/rfd4_smoke_${SLURM_JOB_ID}"
mkdir -p "$OUT"

echo "=== node $(hostname) | job $SLURM_JOB_ID ==="
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader

echo
echo "=== GATE 1: torch sees the GPU (cu126 vs the old cu132) ==="
"$PY" -c "
import torch
ok = torch.cuda.is_available()
print('torch', torch.__version__, '| cuda build', torch.version.cuda)
print('cuda_available', ok)
if ok:
    print('device', torch.cuda.get_device_name(0))
    x = torch.randn(2048, 2048, device='cuda')
    print('matmul on GPU ok, trace =', float((x @ x.T).diagonal().sum()))
raise SystemExit(0 if ok else 1)
"
G1=$?
echo "GATE 1: $([ $G1 -eq 0 ] && echo PASS || echo FAIL)"

echo
echo "=== GATE 1b: tmol imports and sees CUDA (stage 3 physics scorer) ==="
"$PY" -c "
import tmol, torch
print('tmol imported from', tmol.__file__)
print('torch cuda for tmol:', torch.cuda.is_available())
" 2>&1 | tail -5
echo "GATE 1b exit: $?"

echo
echo "=== GATE 2: RFD4 preflight ==="
cd "$RFD4" || exit 1
VIRTUAL_ENV="$RFD4/.venv-cu12" PATH="$RFD4/.venv-cu12/bin:$PATH" \
    bash skills/rfd4-proteina-setup/scripts/preflight.sh 2>&1 | tail -20
echo "GATE 2 exit: $?"

echo
echo "=== GATE 3: one real generation, batch_size=1 ==="
# The documented representative smoke example from the repo's AGENTS.md.
VIRTUAL_ENV="$RFD4/.venv-cu12" PATH="$RFD4/.venv-cu12/bin:$PATH" \
    "$RFD4/.venv-cu12/bin/rfd4-generate" \
        inputs=data/benchmarks/monomer/monomer-short.cif.zst \
        num_replicates=1 batch_size=1 \
        out_dir="$OUT/monomer" 2>&1 | tail -40
echo "GATE 3 exit: $?"

echo
echo "=== outputs produced ==="
find "$OUT" -type f \( -name "*.cif*" -o -name "*.csv" \) | head -20
echo "=== done ==="
