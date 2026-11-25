# Kemp motif scaffolding: RFD3 / RFD4 trial runs (2026-11-21)
Motif: HG3.17 W44 ring, Q50 side chain, D127 carboxylate (unindexed, fixed tip atoms + sequence) + 6NT ligand (fixed). Scaffold 150-260 aa.

## Layout
- rfd3/: inputs (kemp_rfd3_150_260.json = 1 key per loopy/non-loopy), configs/, run_arms_rfd3.py, run_rfd3_arms.sbatch, eval_kemp.py, out/20261121_{loopy,nonloopy}_{unsteered,pocket,pocket_loops}/, out/20261121_eval.csv
- rfd4/ (git clone, .venv, weights in HF cache; configs/inference_experiment/kemp_*.yaml; custom metrics kemp_pocket.py, kemp_tmol.py)
- rfd4_runs/: inputs/ (kemp_rfd4.yaml, kemp_motif.cif, 6NT.tmol), run_arms.py, run_rfd4_arms.sbatch, eval_rfd4.py, out/20261121_{unsteered,pocket,pocket_loops,tmol}/, out/20261121_eval_rfd4*.csv
- logs/: slurm logs. RFD4 install deviation: bionemo-ir dropped from the pip extras (needs GLIBC 2.34, Helix has 2.28) -> native backend.

## Results (n small: 4 replicates RFD4, 12 particles x 1 trajectory RFD3; indicative only)
RFD4 (same lengths per replicate across arms; steered arms return 2 beam members per replicate = 4 lineages)
| arm | n | volume_ratio | ratio error | loop frac in 8 A shell | TMol ligand ddG |
| unsteered | 4 | 0.86+-1.03 | 0.82 | 0.27 | 0.2+-23.5 |
| pocket | 8 | 0.95+-0.08 | 0.07 | 0.26 | -4.1+-6.0 |
| pocket+loops | 8 | 0.95+-0.19 | 0.13 | 0.50 | -6.8+-4.7 |
| tmol(+0.5 pocket) | 8 | 0.98+-0.13 | 0.11 | 0.26 | -13.1+-2.7 |
RFD3 (fork FK steering, backbone-only proxy): steering made pockets too tight (ratio 0.53-0.70 vs 1.04-1.11 unsteered); ESS collapse to 1-2 of 12 particles; needs proxy calibration + tau calibration.
Caveat: RFD4 reward and evaluation use the same metric (circular); pocket proxy = grid free volume (3 A clearance, ray enclosure), not fpocket.

## Round 2: theozyme GDM_KE_1483 (Glu base, Ser partner, 6NT; all file atoms fixed), burial-conditioned (2026-11-21, RFD4 only)
Specs: rfd4_runs/inputs/kemp2_rfd4_{plain,buried}.yaml (buried: C_RAS BURIED on motif atoms + whole ligand, C_HPHOB on base CG/CD and ligand ring C).
Steered: configs/inference_experiment/kemp2_{tmol,chem}.yaml (chem = TMol dG + pocket + base SASA + base apolar env + nitro H-bonds + base H-bond partners).
Metrics: custom_metrics/kemp_chem.py. Evaluation: rfd4_runs/eval_k2.py -> out/20261121_eval_k2.csv. n = 8 designs per arm (4 lineages for steered).
| arm | base SASA rel | lig SASA | base apolar frac | nitro H-bonds | nitro O buried+unsat | vol ratio | TMol dG |
| unsteered, no burial | 0.32 | 0.25 | 0.64 | 0.6 | 0.25 | 0.38 | -9.7 |
| unsteered, burial | 0.07 | 0.02 | 0.66 | 0.9 | 1.1 | 0.89 | -6.7 |
| steered TMol | 0.07 | 0.00 | 0.67 | 1.4 | 0.6 | 0.74 | -17.9 |
| steered chem | 0.03 | 0.01 | 0.80 | 1.0 | 1.0 | 0.73 | -18.3 |
