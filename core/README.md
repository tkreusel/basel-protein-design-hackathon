# core — shared library, calibration and tests

| Path | What |
|---|---|
| `src/kemp/` | geometry library (`geometry.py`: roles from the atom map, distances, burial, ligand-copy discovery) and the `DesignRecord` data class |
| `src/derive_atom_map.py`, `check_motif.py`, `check_burial.py`, `milestone1.py`, `explore_benchmark.py` | derive the 6NT atom map from crystals, score a motif/design, check ligand burial, reproduce the Milestone 1 validation table |
| `src/motif_to_cif.py`, `motif_locate.py`, `make_his_motif.py` | build RFD4 inputs from the theozyme; locate unindexed motif residues in RFD4 output; His-base variant |
| `src/extract_designs.py`, `finalize_run.py`, `finalize_refold.py`, `his_gate.py`, `collect_passing.py`, `make_shortlist.py` | turn raw RFD4/Boltz output into the human-readable run-folder layout and gate/collect successes |
| `src/run_ligandmpnn.py`, `run_stage2_ligandmpnn.py`, `prepare_refold.py`, `prepare_boltz_both.py`, `analyze_refold.py`, `confidence_report.py`, `compare_*.py`, `colabfold_msa.py` | sequence design, refold inputs, refold analysis and arm/MSA comparisons |
| `src/steering_criteria.py`, `apply_steering_criteria.py`, `base_environment.py`, `substrate_access.py` | our reimplementation of the steering campaign's base-environment criteria (validated: 198/198 of its designs reproduce), substrate access probe |
| `src/host_geometry.py` | host-style approach geometry and base burial (used by `pipelines/07`) |
| `configs/atom_map.yaml` | which 6NT atom plays which catalytic role (derived from 11 ligand copies) |
| `configs/thresholds.yaml` | every cut-off, annotated with what it was calibrated against |
| `configs/host_reference.yaml` | crystal-derived ideal geometry for the stand-in score |
| `ligands/6NT_ccd.cif` | 6NT definition from the CCD |
| `env/paths.sh` | central path and cache configuration (source first in every job) |
| `tests/test_geometry.py` | 30 tests encoding the real observed values (donor separates active/inactive; base distance anti-correlates with activity) |

Run: see `docs/03_REPRODUCING.md` §0–1. Tests need the benchmark crystals: `scripts/fetch_benchmark.sh <workdir>`, then
`KEMP_ROOT=<workdir> PYTHONPATH=src python -m pytest tests -q`.
