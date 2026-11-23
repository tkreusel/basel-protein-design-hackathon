# Partial diffusion around our best designs (2026-11-22)

Idea (Tobias): re-noise good designs with RFD4 and regenerate the tail of the
trajectory, so good designs get better for the hosts' ranking and the steered and
unsteered designs alike get extra refinement steps.

## What was decided, and why

- **Best-of-N with the real metric, not steering on a proxy.** The hosts score the
  Boltz-2 *refold* (iPDE, geometry, consistency). Only burial is a backbone property.
  So generation is unsteered; the only "steering" is a CPU rejection filter
  (`select_children.py`) on backbone-level properties before any GPU time is spent:
  motif tips and ligand intact, at least the parent's number of Maxi's core criteria,
  base burial (CA within 12 A) not below the parent's. Survivors then get the whole
  MPNN -> Boltz -> `scoring/host_score.py` chain, and only children that beat their
  parent are kept.
- **Parents** = annotated `*_generated.cif` written by RFD4 (carries the conditioning
  columns). Boltz refolds cannot be parents. Parents are ranked on composite3 =
  mean(iPDE, geometry, burial scores) so single-seed (overnight) and multi-seed
  (pool) designs are comparable.
- **Our own job 21/23 outputs are excluded**: their generated CIFs have no C_CRD
  (fixed-atom) column, so the motif would not stay pinned. Listed with the reason in
  `parents_excluded.csv`. (Maxi's designs and the overnight run have it.)
- **Noise levels** 0.3 and 0.5 (`partial_trajectory_fraction_steps`: 0.3 = rerun the
  last 30 % of the trajectory = close to the parent). No evidence yet which is better;
  both are run and the results are reported per level.

## Files

| File | Role |
| --- | --- |
| `run_partial.sbatch` | the whole chain on one GPU, idempotent, stage times in STAGES.log |
| `merge_rankings.py` | puts pool (multi-seed) and overnight (single-seed) scores on composite3 |
| `select_parents.py` | best backbone per design; finds its annotated RFD4 CIF; writes the RFD4 spec |
| `select_children.py` | CPU rejection filter, top M per parent |
| `order_sequences.py` | MPNN sequences only, ordered so every parent's best children refold first |
| `compare_children.py` | child vs parent on composite3; picks the multi-seed set |

## Tested / not tested

- tested on CPU: parent resolution for all 106 pool candidates (22 usable parents, 10
  excluded for missing C_CRD), child filter on 13 pseudo-children, sbatch syntax
- **not tested**: the RFD4 partial-diffusion call itself (GPU was busy), whether the
  motif and ligand survive it, the multi-example spec naming (`p01_<replicate>`).
  The sbatch fails fast with the log tail if level 0.3 produces nothing; check the
  filter's rejection counts first (`filter.log`): "motif-tips-moved" for most children
  would mean conditioning was lost.

## Launch (after jobs 60 and 61)

    sbatch --dependency=afterany:61 --exclude=gpu-dy-l4big-1,gpu-dy-l4small-1,gpu-dy-l4small-2 \
      --export=ALL,KEMP_JOB61=<job 61 run dir>,KEMP_JOB60=<job 60 run dir>,\
    KEMP_BOLTZ_STOP=$(date -d '2026-11-22 09:20 UTC' +%s),KEMP_MS_STOP=$(date -d '2026-11-22 10:30 UTC' +%s) \
      partial/run_partial.sbatch
