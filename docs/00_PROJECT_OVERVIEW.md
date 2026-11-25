# 00 Project overview

## The challenge

The hackathon asked for de novo Kemp eliminases. The Kemp elimination (5-nitrobenzisoxazole to 2-cyano-4-nitrophenolate) needs a base that abstracts the C3 proton and,
ideally, a donor that stabilises the developing negative charge. The hosts supplied the transition-state analogue 6NT and the HG3.17 theozyme; each team could submit up to 20
sequences in one FASTA file, scored at the deadline with re-uploads allowed. The full rules are in [reference/SUBMISSION_hosts_spec.md](reference/SUBMISSION_hosts_spec.md).

| Rule | Value |
|---|---|
| Designs, length | at most 20, 50 to 260 residues, standard amino acids |
| Header | `>name base=<n> oxy=<n>` (1-based; base is Asp or Glu, oxy is Gln, Asn, Ser, Thr, Tyr, Lys or Arg; His is not allowed as donor) |
| Eligibility | below 70 % sequence identity to published enzymes and TM-score of at most 0.7 to published structures |
| Scoring | Boltz-2 co-fold with 6NT and MSAs; catalytic contact distances and approach geometry (terms multiply and dominate), iPDE, base enclosure, pose consistency over repeated predictions |

## The active site

One theozyme, taken unchanged from the hosts' `HG3_17_kemp_motif.pdb` (the 4BS0 pose of HG3.17):

| Residue | Role | Distance to 6NT |
|---|---|---|
| Asp127 | catalytic base | OD2 to N3 2.53 A |
| Gln50 | oxyanion donor | NE2 to N1 3.10 A |
| Trp44 | aromatic wall (optional, the "with Trp" motif) | ring centroids 5.53 A |

The atom map was derived from coordinates over 11 ligand copies and agrees with the hosts' validation code: the abstraction centre is N3 and the oxyanion centre is N1.
The nitro group is not the oxyanion contact (it lies 7.6 to 8.7 A from the donor tip); see `core/configs/atom_map.yaml`.

One result shaped the scoring. Across the HG3 series the catalytic base distance anti-correlates with activity (HG3, the weakest variant, has the shortest, 2.23 A), while the donor
distance discriminates (active Gln50 2.96 to 3.10 A, inactive 6 to 7 A). A score that rewards a shorter base distance therefore ranks backwards;
[reference/MILESTONE1.md](reference/MILESTONE1.md) documents this and a regression test in `core/tests` pins it.

## Workstreams

| Workstream | Product | Location |
|---|---|---|
| Foundation and early generation | theozyme conditioning recipe, geometry and burial scoring with calibrated thresholds, stage 1 to 5 jobs, unsteered control, partial diffusion, host-style scoring, eligibility gate | `core/`, `pipelines/01, 04, 05, 07, 08` |
| Campaign | reward-steered and unsteered RFD4 generation, gated LigandMPNN and Boltz-2 funnel, 10,500 backbones over 6 motif/arm combinations | `pipelines/02_rfd4_steering` |
| Evaluation suite | ipSAE2, PLACER (holo and apo), tmol, PyRosetta, extra features, literature validation | `pipelines/06_evaluation_suite` |
| Ensembles and literature | BioEmu-2 ensembles, literature correlation | `pipelines/09_bioemu` |

## Timeline (2026)

| Date | Milestone |
|---|---|
| 11-21 | cluster discovery, tool installs, theozyme conditioning recipe working, milestone 1 scoring validation, first 32-design RFD4 batches, evaluation suite built |
| 11-21 to 22 | overnight: the unsteered control (250 backbones), multi-seed Boltz scoring and partial diffusion on the hackathon cluster; the steered and unsteered campaign on the large cluster |
| 11-22 morning | full evaluation of finalists, BioEmu analysis, literature correlation (no signal after selection-aware correction), eligibility infrastructure |
| 11-22 noon | final selection (merged top 80, eligibility, ranking, top 20), submitted before 12:00 UTC |

## Outcome

The 20 submitted sequences lie on 10 distinct backbones and won the hackathon; the list with the generation method of each is in [06_RESULTS_AND_FUNNELS.md](06_RESULTS_AND_FUNNELS.md).
The ranking rests on an in-silico stand-in for the hosts' score and a weak literature calibration, and no design has been tested experimentally.

## Pitfalls worth knowing first

1. The default `.venv` of a fresh RFD4 clone may carry a torch build newer than the cluster driver; the GPU is then invisible. Check before generating (see [02](02_ENVIRONMENTS_AND_CLUSTER.md)).
2. RFD4 motif conditioning fails silently when mis-specified; use the recipe in [reference/TOOL_NOTES.md](reference/TOOL_NOTES.md).
3. Boltz-2 without an MSA often puts the ligand somewhere other than the designed pose; this was the largest single loss in every funnel.
4. Count lineages, not designs, for steered runs: the three beam members of a replicate are near-identical (CA RMSD about 0.03 A).
5. Pass criteria that equal the steering rewards are circular, so passing them is not evidence of catalytic competence.
6. tmol relaxes are stochastic (about 5 REU) and the PyRosetta Sc term is unreliable.
7. Design identifiers were truncated to 50 characters in one eligibility table and collided; join by sequence.
8. The UniRef90 search needs 28 to 30 GB of RAM and can take an hour; start it early.
9. Project sequences and structures were never sent to external services; all eligibility checks use local databases.
