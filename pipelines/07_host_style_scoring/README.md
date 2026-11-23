# Host-style scoring pipeline (2026-11-22 submission, noon deadline)

Stand-in for the hosts' internal ranking, built from their scoring slide.
**The hosts' formulas, weights and cut-offs are not public.** Every definition
below is ours, calibrated on HG3-series crystals; component values are always
reported so the ranking can be re-weighted the moment more is known.

| Component | Slide wording | Our definition |
| --- | --- | --- |
| iPDE | Boltz-2 metric, protein-ligand token distance error; "is the ligand confidently placed?" | Boltz's own `complex_ipde`; score `clip((4 - ipde)/3.5, 0, 1)` (our 425 refolds: 0.42-5.19, median 1.20) |
| Geometry | distance/angle/dihedral vs the ideal approach geometry for the base and H-bond donor | 10 measurements (5 base O->N3, 5 Gln NE2->N1); ideal = 8 crystal copies with Asp127+Gln50; per-metric `exp(-0.5 z^2)`, sigma = max(sd, floor); mean. Crystals score 0.89-0.95 |
| Burial | number of CA atoms within 12 A of the catalytic base | CA count around the base CG (also around the OD midpoint, reported); score `min(1, n / 45.9)` |
| Consistency | metric reproducibility across Boltz-2 seeds; proxy for preorganisation | needs >= 2 seeds: median pairwise ligand RMSD (after CA superposition, weight 2), plus spread of base distance, donor distance, burial, iPDE |

`composite_abs` = mean of the available component scores.
`composite_rank` = mean percentile rank within the pool.

## Files

| File | Role |
| --- | --- |
| `src/host_geometry.py` | geometry + burial measurements |
| `scoring/calibrate_host_reference.py` -> `configs/host_reference.yaml` | the ideal, measured from crystals |
| `scoring/pool_candidates.py` | gather candidates (ours, Maxi's steered, overnight) |
| `scoring/boltz_seeds.sbatch` | multi-seed Boltz-2 (separate seeds, chunk-outer so finished candidates are complete) |
| `scoring/host_score.py` | scores, ranks, flags |

## Tested

- crystals through the geometry scorer: 0.95 (7K4Z), 0.89 (7K4X)
- 124 existing refolds: scored in 2 s; iPDE correlates with geometry (rho -0.47)
- consistency, on constructed pseudo-seeds: rigid motion -> 1.00; ligand shifted 1.5 A -> pose RMSD 1.50
- real multi-seed Boltz run (job 61): output layout `seed_<k>/boltz_results_chunk_<n>/predictions/<id>/` is found by the scorer
- **bug found and fixed 2026-11-22 04:30 UTC**: `geom_scores` returned per-metric scores under the SAME keys as the raw measurements (`b_dist`, `d_dist`, ...), so the raw distances were overwritten by 0-1 scores. Geometry/burial/iPDE scores were never affected, but the consistency terms for base/donor distance used the spread of scores divided by 0.25 A (inflating consistency) and the reported `b_dist`/`d_dist` columns were scores. Per-metric scores are now `sc_<metric>`. Any ranking made before this fix should be recomputed (the final job-61 scoring runs the fixed file).

## Known limits

- Asp base / Gln donor only. His, Zn or GDM Glu/Ser designs cannot be scored (no crystal reference).
- Our "ideal" comes from one scaffold family; the hosts' may differ.
- Burial counts depend on protein size: the crystal scaffold is 300 residues, our designs 150-260.
- Consistency rewards reproducibly bad designs too; it is combined with the other three, not used alone.
