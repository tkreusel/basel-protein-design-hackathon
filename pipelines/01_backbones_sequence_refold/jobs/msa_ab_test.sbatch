#!/usr/bin/env bash
#SBATCH --job-name=kemp-msatest
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=00:40:00

# A/B test: does a real ColabFold MSA change the Boltz-2 refold?
#
# Uses Boltz's OWN MSA path (--use_msa_server, default server
# https://api.colabfold.com), not a hand-rolled client: boltz/data/msa/mmseqs2.py
# does the submit / poll / download round-trip itself. This is also what the
# organisers' evaluation will do, so matching it matters.
#
#   arm "with_msa": yaml has NO msa field  + --use_msa_server
#                   -> Boltz fetches the alignment from ColabFold
#   arm "no_msa"  : yaml has msa: empty    (what our refold runs used)
#
# Everything else identical: same sequence, same seed, same flags.
#
# Expectation to falsify: a direct ColabFold query for this design returned
# depth=1, nonquery_rows=0 (query-only), matching the organisers' own example
# for their de novo design. If so the two arms carry identical information and
# should be indistinguishable. Measure it rather than assume it.
#
# Required: KEMP_DESIGN (design id)

set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh

DID="${KEMP_DESIGN:?set KEMP_DESIGN}"
SRC_RUN="${KEMP_SRC_RUN:-$KEMP_ROOT/runs/2026-11-21_rfd4-backbones-BURIED_job21_32designs_seed201}"
PY="$KEMP_ROOT/env/kemp-orch/bin/python"
W="$KEMP_RUNS/raw_msa-abtest_job${SLURM_JOB_ID}"
export TMPDIR="/tmp/kemp-$SLURM_JOB_ID"
mkdir -p "$TMPDIR" "$W/with_msa" "$W/no_msa"
trap 'rm -rf "$TMPDIR"' EXIT

echo "=== node $(hostname) | job $SLURM_JOB_ID | design $DID ==="
SEQ=$("$PY" -I -c "
import sys; sys.path.insert(0,'$KEMP_ROOT/src')
from kemp.record import read_manifest
for r in read_manifest('$SRC_RUN/manifest.jsonl'):
    if r.design_id=='$DID': print(r.sequence)
")
echo "sequence length: ${#SEQ}"

# with_msa: omit the msa key entirely so Boltz generates one from the server.
cat > "$W/with_msa/${DID}.yaml" <<EOF
version: 1
sequences:
- protein:
    id: A
    sequence: $SEQ
- ligand:
    id: B
    ccd: 6NT
EOF

# no_msa: single-sequence mode, exactly what the refold runs used.
cat > "$W/no_msa/${DID}.yaml" <<EOF
version: 1
sequences:
- protein:
    id: A
    sequence: $SEQ
    msa: empty
- ligand:
    id: B
    ccd: 6NT
EOF

BC="$W/boltz-cache"; mkdir -p "$BC"
for ASSET in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do
    [ -e "$BC/$ASSET" ] || ln -s "$BOLTZ_ASSETS/$ASSET" "$BC/$ASSET"
done
[ -e "$BC/mols" ] || ln -s "$BOLTZ_ASSETS/mols" "$BC/mols"

COMMON=(--model boltz2 --cache "$BC" --accelerator gpu --devices 1 --no_kernels
        --sampling_steps 200 --recycling_steps 3 --diffusion_samples 1
        --num_workers 0 --preprocessing-threads 1 --output_format mmcif --seed 1)

echo
echo "=== arm with_msa (Boltz fetches from ColabFold) ==="
"$ENV_ROOT/env-boltz2/bin/boltz" predict "$W/with_msa" \
    --out_dir "$W/out_with_msa" "${COMMON[@]}" \
    --use_msa_server --msa_server_url https://api.colabfold.com 2>&1 | tail -10
echo "rc=$?"

echo
echo "=== arm no_msa (single sequence) ==="
"$ENV_ROOT/env-boltz2/bin/boltz" predict "$W/no_msa" \
    --out_dir "$W/out_no_msa" "${COMMON[@]}" 2>&1 | tail -6
echo "rc=$?"

echo
echo "=== the alignment Boltz actually fetched ==="
find "$W" -name "*.a3m" | head -3 | while read -r f; do
    echo "  $f"
    echo "    rows: $(grep -c '^>' "$f")"
done

echo
echo "=== COMPARISON ==="
KEMP_SRC_RUN="$SRC_RUN" "$PY" -I "$KEMP_ROOT/src/compare_msa_arms.py" \
    "$KEMP_ROOT" "$W" "$DID"
echo "=== done ==="
