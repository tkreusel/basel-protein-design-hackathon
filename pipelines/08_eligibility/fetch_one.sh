#!/usr/bin/env bash
# download one PDB entry (mmCIF, gz); skip if present; 3 tries
id=$1; out=/mnt/work/projects/tobi_v1/eligibility/db/pdb30_cif/$id.cif.gz
[ -s "$out" ] && exit 0
for i in 1 2 3; do
  code=$(curl -s -m 60 -o "$out.part" -w "%{http_code}" "https://files.rcsb.org/download/$id.cif.gz")
  if [ "$code" = 200 ] && [ -s "$out.part" ]; then mv "$out.part" "$out"; exit 0; fi
  sleep 2
done
rm -f "$out.part"; echo "$id failed" >> /mnt/work/projects/tobi_v1/eligibility/db/fetch_failed.txt
