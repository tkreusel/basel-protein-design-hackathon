# 09 — BioEmu-2 ensembles and literature correlation

| File | Role |
|---|---|
| `run_bioemu.slurm` | `bioemu2 sample-equilibrium --system <system.yaml> --num-samples N` on one L4 (system = protein + 6NT) |
| `analyze_ensembles.py` | per run: frames.csv (CA RMSD, Rg, ligand pose RMSD/COM shift/tilt, host-style approach geometry, burial, SASA), residues.csv (core-aligned RMSF, contact fraction), summary.json; samples are independent (no time axis) and are superposed on an iteratively chosen rigid core |
| `plot_ensembles.py`, `ligand_rmsf.py`, `ligand_density_grid.py` | figures: ensembles, ligand RMSF, ligand density |
| `ENSEMBLE_REPORT.md` | findings on three campaign designs: protein rigid, ligand bimodal (15–23 % bound), base end swings 2.5–2.8 Å |
| `literature/ensemble_metrics.py`, `correlate.py`, `correlations.csv`, `REPORT.md` | ensemble metrics for literature enzymes vs kcat/KM; a **selection-aware permutation test** (labels shuffled 20,000×, max |ρ| over the metric family) — no metric survives (best raw p 0.011 → 0.18) |

Caveats: BioEmu samples are equilibrium ensembles without a time axis; the 1A53-core ensembles were invalid for these metrics (Val at the base position) and were excluded.
