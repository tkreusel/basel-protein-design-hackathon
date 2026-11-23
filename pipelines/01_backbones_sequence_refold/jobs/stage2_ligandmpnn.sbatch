#!/usr/bin/env bash
#SBATCH --job-name=kemp-mpnn
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=01:00:00

# Stage 2: LigandMPNN re-design of stage-1 success backbones.
#
# Why this matters here: RFD4's own co-designed sequences fold correctly
# (25/31 at TM>=0.85) but pin the oxyanion donor rotamer in only 2/31 refolds.
# LigandMPNN is explicitly ligand-aware and packs side chains against the
# ligand, which is exactly the failure mode. This is also the plan's actual
# sequence-design stage -- RFD4's sequence was always a bonus, not the plan.
#
# Catalytic residues are FIXED and the output is VERIFIED to confirm they were
# not redesigned (plan 6.3: "assert this, do not trust it").
#
# Required: KEMP_RUN = finalized stage-1 run folder
# Optional: KEMP_NSEQ (sequences per backbone, default 4)

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

RUN="${KEMP_RUN:?set KEMP_RUN}"
NSEQ="${KEMP_NSEQ:-4}"
TAG="${KEMP_TAG:-mpnn}"
PY="$KEMP_ROOT/env/kemp-orch/bin/python"
W="$KEMP_RUNS/raw_${TAG}_job${SLURM_JOB_ID}"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"
mkdir -p "$TMPDIR" "$W"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== node $(hostname) | job $SLURM_JOB_ID ==="
echo "source run: $RUN"
echo "sequences per backbone: $NSEQ"
nvidia-smi --query-gpu=name --format=csv,noheader

"$PY" -I "$KEMP_ROOT/src/run_ligandmpnn.py" "$KEMP_ROOT" "$RUN" "$W" "$NSEQ"
RC=$?
echo "ligandmpnn rc=$RC"

echo
echo "=== prepared for refold ==="
echo "boltz inputs: $(ls "$W/boltz_in"/*.yaml 2>/dev/null | wc -l)"
echo "work dir: $W"
echo "=== done ==="
