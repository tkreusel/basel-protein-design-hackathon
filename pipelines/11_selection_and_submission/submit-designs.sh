#!/usr/bin/env bash
# Submit your team's designs to the hackathon.
#
#   ./submit-designs.sh designs.fasta
#
# Validates the file first and refuses to upload anything that would be
# rejected by the scoring pipeline. Safe to run as many times as you like:
# each upload replaces the previous one, and only your last upload is scored.
set -euo pipefail

BUCKET="perlin-pdh-2026-paris"
REGION="eu-west-3"
HERE="$(cd "$(dirname "$0")" && pwd)"
# Find a Python 3 that actually runs. On Windows, "python3" can resolve to a
# Microsoft Store alias stub that exits non-zero, so test each candidate.
PY_BIN=""
for _cand in python3 python py; do
    if command -v "$_cand" >/dev/null 2>&1        && "$_cand" -c 'import sys; sys.exit(0 if sys.version_info[0]==3 else 1)' >/dev/null 2>&1; then
        PY_BIN="$_cand"; break
    fi
done
if [[ -z "$PY_BIN" ]]; then
    echo "Python 3 is required but no working python3 was found on PATH." >&2
    exit 2
fi


FASTA="${1:-}"
if [[ -z "$FASTA" ]]; then
    echo "usage: $0 <designs.fasta>" >&2
    exit 2
fi
if [[ ! -f "$FASTA" ]]; then
    echo "No such file: $FASTA" >&2
    exit 2
fi

echo "==> Validating $FASTA"
"$PY_BIN" "$HERE/validate-submission.py" "$FASTA" || {
    echo
    echo "Not uploaded. Fix the problems above and run this again." >&2
    exit 1
}

echo
echo "==> Identifying your team account"
ACCOUNT=$(aws sts get-caller-identity --query Account --output text | tr -d '\r')
if [[ -z "$ACCOUNT" || "$ACCOUNT" == "None" ]]; then
    echo "Could not read your AWS account id. Are you signed in?" >&2
    exit 1
fi
KEY="submissions/${ACCOUNT}/designs.fasta"
echo "    account ${ACCOUNT}"
echo "    target  s3://${BUCKET}/${KEY}"

echo
echo "==> Uploading"
aws s3api put-object \
    --bucket "$BUCKET" --region "$REGION" \
    --key "$KEY" --body "$FASTA" \
    --content-type "text/plain" \
    --query 'VersionId' --output text | tr -d '\r' | sed 's/^/    version /'

echo
echo "SUBMITTED. Your designs are in. Re-run this any time before the deadline"
echo "to replace them; the file that is present at the deadline is the one scored."
