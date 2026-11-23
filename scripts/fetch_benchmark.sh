#!/usr/bin/env bash
# Download the literature Kemp-eliminase crystal structures used as the calibration/benchmark set (public PDB entries, mmCIF).
# Only PDB ids go out; nothing from the project is sent anywhere.
#   scripts/fetch_benchmark.sh <workdir>      # fills <workdir>/benchmark/{pdb,apo}
set -euo pipefail
W="${1:?usage: fetch_benchmark.sh <workdir>}"; mkdir -p "$W/benchmark/pdb" "$W/benchmark/apo"
# 6NT-bound (holo) entries used by tests, thresholds and the host-style reference (docs/reference/LITERATURE_STRUCTURES.md)
HOLO="4BS0 6NW4 6SRW 6SS1 6SS3 6TU6 7K4Q 7K4U 7K4X 7K4Z 8RD5 9TZD"
APO="7K4T 7K4W 9U01 "
for id in $HOLO; do [ -s "$W/benchmark/pdb/$id.cif" ] || curl -sf -m 120 -o "$W/benchmark/pdb/$id.cif" "https://files.rcsb.org/download/$id.cif" || echo "failed: $id"; done
for id in $APO; do [ "$id" = "__APO__" ] && break; [ -s "$W/benchmark/apo/$id.cif" ] || curl -sf -m 120 -o "$W/benchmark/apo/$id.cif" "https://files.rcsb.org/download/$id.cif" || echo "failed: $id"; done
ls "$W/benchmark/pdb" | wc -l
