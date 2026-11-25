#!/bin/bash
# Regenerate the campaign figures from $KEMP_RESULTS/tables (run `python -m kemp.export_results` first).
# usage: KEMP_RESULTS=/path/to/results PY=/path/to/python bash figures/make_all.sh
set -e
cd "$(dirname "$0")"; PY=${PY:-python3}; : "${KEMP_RESULTS:?set KEMP_RESULTS}"
for s in f01_funnel f02_boltz_confidence_hist f03_tsa_vs_substrate_product_hist f04_energy_placer_hist f05_ipsae_vs_dGsep_scatter \
         f06_placer_holo_vs_apo_scatter f07_channel_vs_tmol_designs f08_boltz_gate_criteria f09_scaffold_population d_design_metric_hists e_steering_objectives; do
  $PY "$s.py" 2>&1 | grep -v -i warn | tail -1
done
