#!/usr/bin/env bash
#SBATCH --job-name=kemp-uniref-id
#SBATCH --partition=gpu
#SBATCH --exclude=gpu-dy-l4big-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --gres=gpu:l4:1
#SBATCH --time=01:30:00
# Sequence identity of the candidates against UniRef90 (2022-05) with the GPU MMseqs2 DB (needs ~28 GB RAM -> l4-1).
# Same verified invocation as the example pipeline; output TSV is read by check_eligibility.py
set -o pipefail
source /mnt/work/projects/tobi_v1/env/paths.sh
cd /mnt/work/projects/tobi_v1/eligibility; mkdir -p results tmp/uniref
python3 - <<'PY'
import csv
rows=list(csv.DictReader(open("candidates.csv")))
open("tmp/uniref/q.fasta","w").write("".join(f">{r['cand_id']}\n{r['sequence']}\n" for r in rows))
PY
"$MMSEQS" createdb tmp/uniref/q.fasta tmp/uniref/queryDB --dbtype 1 --shuffle 0
"$MMSEQS" search tmp/uniref/queryDB "$MMSEQS_DB.idx" tmp/uniref/resultDB tmp/uniref/tmp --db-load-mode 2 --gpu 1 \
    --max-seqs 500 -e 1e-3 --threads 4
"$MMSEQS" convertalis tmp/uniref/queryDB "$MMSEQS_DB.idx" tmp/uniref/resultDB results/uniref90_hits.tsv \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov" --threads 4
echo "done: $(wc -l < results/uniref90_hits.tsv) hits"
