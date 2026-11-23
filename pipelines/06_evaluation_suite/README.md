# 06 — Evaluation suite (ipSAE2 · PLACER · tmol · PyRosetta)

`suite/` is the suite owner's code, unmodified (code, config, jobs, tests; no `work/`). Documents: `PIPELINE.md` (stages, tmol protocol per replicate, what is native vs proxy,
findings along the way), `PROPOSAL.md`, `RESULTS.md`, `PLACER_RESULTS.md`.

```
Boltz-2 refolds (6NT / H5J / 6VP per protein)
   └─ jobs/00_prepare_and_ipsae.sh   standardise inputs (protein chain A, ligand chain B), ipSAE2 for every arm        CPU, minutes
        ├─ jobs/10_placer*.sbatch    PLACER holo (each ligand) + apo ensembles, 20–50 models                           GPU
        ├─ jobs/20_tmol / 21_tmol_cpu  tmol protocol per replicate (raw, rlx, rlxcst, tsa_pk/sub/prod, …)             head-node CPU, sharded
        ├─ jobs/40_pyrosetta / 41     relax, interface dG_sep, Sc (unreliable), packstat, BUNS, ΔΔG                      CPU
        ├─ jobs/30_placer_analysis    ensemble statistics + tmol electric field                                         GPU/CPU
        └─ aggregate → work/report/{summary_entries.csv, summary_replicates.csv, report.md}
jobs/run_full_suite.sh --root <boltz folder> [--submit]    dry-run by default; prints readiness, queue, plan
```

`our_additions/` — what we added to run it on our finalists:

| File | Role |
|---|---|
| `stage_finalists.py` | copies the 6NT refolds of the finalists into `boltz_root/<id>/` and writes Boltz inputs for the substrate (H5J) and product (6VP) arms |
| `stage_extra.py` | the same for the extra pool designs never evaluated |
| `arms_then_placer.sbatch`, `arms2.sbatch` | Boltz arm refolds → prepare + ipSAE2 → PLACER; second batch of arms only |
| `geometry_wall_ring_fallback.patch` | `geometry.py` returns `{"geo.wall_found": 0.0}` when the wall ring has < 3 points (PLACER crops lack wall-ring atoms and the analysis crashed) |
| `plot_best.patch`, `env.sh.as-used`, `finalists.csv`, `finalists_report.md` | plotting tweak, environment as used, the 12 finalists, the suite report for them |

Caveats the owner states: composite/rank are uncalibrated heuristics; tmol relax is stochastic (~5 REU); Sc unreliable; PLACER is a weak, noisy feature
(HG3.17 crystal holo ≈ 0.05; within-backbone variance ≈ 50 %); 6NT is a neutral TSA analogue and 6VP in the CCD is the neutral phenol.
