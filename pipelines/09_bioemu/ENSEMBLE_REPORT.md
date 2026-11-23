# BioEmu-2 ensembles of three Maxi designs + 6NT (2026-11-22)

Runs analysed (1000 samples each, `bioemu2 sample-equilibrium`, produced on Maxi's cluster):
`hgnt__steer_base_tmol__no_trp_0_beam3_5__mpnn1`, `..._beam3_5__mpnn4`, `..._beam3_6__mpnn1`
(no-Trp HG3.17 motif, 127 aa, base Asp77, donor Gln50; 5_mpnn1 and 6_mpnn1 differ at a few
positions, 5_mpnn4 is another LigandMPNN sequence on the same backbone).
`attempt24_ce941bbd5c43_SUCCESS` has an empty `partial/` folder (run not finished): not analysed.

Samples are independent draws, not a trajectory: "motion" = spread of the ensemble. Everything is in the
frame of the RFD4 design (CA superposition on an iteratively chosen rigid core, 101-109 of 127 residues).
Script `analyze_ensembles.py`, figures `plot_ensembles.py`; per-run `frames.csv`, `residues.csv`,
`summary.json`, `medoid.pdb`, `ensemble_20frames.pdb` (open in PyMOL/ChimeraX) in `analysis/<run>/`.
All 3000 samples passed BioEmu's physicality check.

## Findings

| | 5_mpnn1 | 5_mpnn4 | 6_mpnn1 |
| --- | --- | --- | --- |
| CA RMSD to design, core (median / p90) | 0.83 / 1.09 | 0.83 / 1.23 | 0.81 / 1.09 |
| CA RMSD to design, all residues | 1.51 / 2.07 | 1.28 / 2.07 | 1.52 / 2.18 |
| median CA RMSF (max, at C-terminus) | 0.70 (6.6) | 0.74 (6.0) | 0.73 (7.6) |
| helix fraction | 0.79 | 0.81 | 0.78 |
| ligand in designed pose (RMSD <= 3 A) | 15.9 % | 15.2 % | 22.9 % |
| bound samples: base O->N3, median (design 2.67) | 7.5 A | 5.6 A | 7.5 A |
| bound samples: Gln NE2->N1, median (design 3.22) | 3.75 A | 3.97 A | 3.48 A |
| bound samples with base <= 3.5 A / both <= 3.5 A | 4 % / 0 % | 13 % / 8 % | 3 % / 2 % |
| bound geometry score, median (design 0.87) | 0.20 | 0.28 | 0.21 |
| bound: median N3 / Asp O / N1 / NE2 shift from design | 2.8 / 2.5 / 1.5 / 1.3 | 2.6 / 0.9 / 1.6 / 1.4 | 2.5 / 2.7 / 1.5 / 1.2 |
| base burial, CA within 12 A (design 41) | 32 | 35 | 32 |

1. **The fold is rigid.** Helical core stays within ~0.8 A of the design; only the termini move (RMSF up to 6-8 A at the
   C-terminus, 5-6 A at the N-terminus). Motif residues are among the stiffest (RMSF 0.45-0.8 A).
2. **The ligand is bimodal.** 15-23 % of samples keep it in the designed pocket (RMSD 2 A); the rest sit 12-20 A away
   (on the surface, closest protein atom still ~3 A), with a long tail to 50 A. The occupancy is the key
   number: the pocket is only partially preferred by the model.
3. **Even when bound, the catalytic geometry is lost.** The N1 / Gln NE2 end stays near the design (1.2-1.6 A shift,
   donor distance 3.5-4.0 A, up to half the samples <= 3.5 A), but the N3 end swings ~2.5-2.8 A away and, for 5_mpnn1 and
   6_mpnn1, the Asp carboxylate moves ~2.5 A too, so base->N3 is 7.5 A instead of 2.7 A. 5_mpnn4 keeps Asp
   in place (0.9 A shift) and is clearly the best of the three (13 % of bound samples with base <= 3.5 A, 8 % with both).
4. **Pocket residues** (bound): Gln50, Leu70, Ser57, Leu53, Ile27, Phe17, Phe73, Ile13 (>= 0.7 of bound samples), Asp77
   itself 0.6-0.7. Bound ligand SASA is ~0 (sealed pocket, as in the design).
5. **Backbone modes** (PCA, 23-27 % in PC1) are dominated by the termini and correlate with base burial
   (PC1 r = -0.4 to -0.6); they do not correlate with the ligand pose (|r| < 0.1).

## Caveats

- No reference system (a crystal-validated Kemp enzyme, or a non-binder) was run, so 15-23 % occupancy and the geometry
  loss cannot yet be called good or bad in absolute terms; the three runs are only comparable with each other.
- BioEmu side chains and the ligand placement are model output, not physics: a 2-3 A displacement of the N3 end is
  within what such a model may get wrong. The comparison with the Boltz refold of the same sequences has not been done.
- "Bound" (pose RMSD <= 3 A) is our cut-off; the conditional statistics change with it.
- n = 3 sequences; the 5_mpnn4 advantage is one sequence on one backbone.
