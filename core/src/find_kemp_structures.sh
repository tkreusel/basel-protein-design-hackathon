#!/usr/bin/env bash
# Reproduce the RCSB full-text search behind docs/LITERATURE_STRUCTURES.md.
# Prints every PDB entry matching "Kemp eliminase", then each entry's title and
# resolution. 6NT presence must be verified on the downloaded coordinates
# (grep -c 6NT), NOT from the title -- several entries use a different analogue.
set -euo pipefail
OUT=${1:-/tmp/kemp_search}
mkdir -p "$OUT"
cat > "$OUT/q.json" <<'JSON'
{"query":{"type":"terminal","service":"full_text","parameters":{"value":"Kemp eliminase"}},
 "return_type":"entry","request_options":{"paginate":{"start":0,"rows":200},"results_verbosity":"compact"}}
JSON
curl -sS -X POST "https://search.rcsb.org/rcsbsearch/v2/query" \
     -H "Content-Type: application/json" --data-binary @"$OUT/q.json" -o "$OUT/hits.json"
python3 -c "
import json,sys
d=json.load(open('$OUT/hits.json'))
ids=d.get('result_set',[])
print('total',d.get('total_count'),'returned',len(ids))
open('$OUT/ids.txt','w').write('\n'.join(ids))
print(' '.join(ids))
"
while read -r id; do
  curl -sS "https://data.rcsb.org/rest/v1/core/entry/$id" -o "$OUT/e_$id.json" || continue
  python3 -c "
import json
d=json.load(open('$OUT/e_$id.json'))
r=(d.get('rcsb_entry_info',{}).get('resolution_combined') or [None])[0]
print('$id', r, d.get('struct',{}).get('title','?')[:90])
"
done < "$OUT/ids.txt"
