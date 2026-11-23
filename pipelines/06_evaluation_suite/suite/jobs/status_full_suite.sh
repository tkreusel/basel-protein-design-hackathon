#!/usr/bin/env bash
# Progress of the full-suite run: bash jobs/status_full_suite.sh
cd "$(dirname "$0")/.."
NT=$(ls work/tmol/*_boltz*_tsa_*.json 2>/dev/null | wc -l); NR=$(ls work/pyrosetta/*_boltz*_tsa_*.json 2>/dev/null | wc -l)
NP=$(python3 - <<'PY'
import json, glob
m = json.load(open('work/inputs/manifest.json'))
want = {e for e in m['entries'] if e.endswith('_boltz')}
print(sum(1 for f in glob.glob('work/placer/*/*.csv') if f.split('/')[-2] in want and not f.endswith('_models.csv')))
PY
)
TOT=$(python3 -c "import json;m=json.load(open('work/inputs/manifest.json'));print(sum(1 for e,v in m['entries'].items() if e.endswith('_boltz') for r in v['replicates'] if r['arm']=='tsa'))")
echo "$(date +%H:%M)  tmol $NT/$TOT TSA replicates | PyRosetta $NR/$TOT | PLACER ensembles $NP/100 (holo x3 arms + apo, Boltz entries)"
squeue -h -o "  job %i %j %T %M" ; [ -f work/report/FULL_SUITE_DONE ] && echo "  FULL SUITE DONE $(cat work/report/FULL_SUITE_DONE) -> work/report/" || echo "  (finisher waiting; writes work/report/FULL_SUITE_DONE)"
