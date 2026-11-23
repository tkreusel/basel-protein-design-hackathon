# 00 — Project overview

## 1. The challenge

A de novo **Kemp eliminase** design hackathon. Kemp elimination (5-nitrobenzisoxazole → 2-cyano-4-nitrophenolate) needs
a catalytic base that abstracts a proton and, ideally, an oxyanion-hole donor that stabilises the developing negative
charge. The hosts supplied the transition-state analogue **6NT** and the HG3.17 theozyme; teams submit up to 20
sequences in one FASTA. The file at the deadline is scored (re-uploads allowed). Full rules:
[reference/SUBMISSION_hosts_spec.md](reference/SUBMISSION_hosts_spec.md).

| Rule | Value |
|---|---|
| Designs / length | ≤ 20; 50–260 aa; 20 standard amino acids |
| Header | `>name base=<n> oxy=<n>` (1-based; base = Asp/Glu, oxy = Gln/Asn/Ser/Thr/Tyr/Lys/Arg; His not allowed as donor) |
| Eligibility | < 70 % identity to published enzymes **and** TM-score ≤ 0.7 to published structures |
| Scoring | Boltz-2 co-fold with 6NT, **with MSAs**; catalytic contact distances + approach geometry (terms multiply, dominate), iPDE, base enclosure, pose consistency across repeated predictions |

## 2. What the active site is

One theozyme, unmodified from the hosts' `HG3_17_kemp_motif.pdb` (pose from PDB 4BS0 / HG3.17):

| Residue | Role | Distance to 6NT |
|---|---|---|
| Asp127 | catalytic base | OD2 → N3 = 2.53 Å |
| Gln50 | oxyanion donor | NE2 → N1 = 3.10 Å |
| Trp44 | aromatic wall (optional "with-Trp" motif) | ring centroids 5.53 Å |

**Atom map** (derived from coordinates over 11 ligand copies, confirmed by the hosts' own validation code):
`abstraction_centre = N3` (the base's target), `oxyanion_centre = N1`. The nitro group is **not** the oxyanion contact
(it sits 7.6–8.7 Å from the donor tip). See `core/configs/atom_map.yaml` and [reference/TOOL_NOTES.md](reference/TOOL_NOTES.md).

Counter-intuitive finding that shaped scoring: over the HG3 series the **catalytic base distance anti-correlates with activity**
(HG3, the weakest variant, has the shortest base distance, 2.23 Å); the *donor* distance is the discriminating metric
(active Gln50 2.96–3.10 Å vs inactive 6–7 Å). A score that rewards "shorter base distance" ranks backwards
([reference/MILESTONE1.md](reference/MILESTONE1.md); pinned by a regression test in `core/tests`).

## 3. How the work was organised

Three workstreams ran in parallel and met in the final selection:

| Workstream | What it produced | Repo location |
|---|---|---|
| **Foundation + our own generation** (this project folder) | theozyme conditioning recipe, geometry/burial scoring, calibrated thresholds, stage-1…5 jobs, unsteered overnight control, partial diffusion, host-style scoring, eligibility gate | `core/`, `pipelines/01,04,05,07,08` |
| **Large-cluster steering campaign** (teammate) | RFD4 reward-steered beam search, custom reward metrics, a 10,500-backbone campaign over 6 motif/arm combinations, full evaluation tables | `pipelines/02_rfd4_steering/` |
| **Evaluation suite** (teammate) | ipSAE2, PLACER (holo + apo), tmol, PyRosetta, extra features, literature validation | `pipelines/06_evaluation_suite/` |
| BioEmu-2 ensembles (teammate + us) | ensemble trajectories and the analysis / literature correlation | `pipelines/09_bioemu/` |

## 4. Timeline (all 2026)

| When | Milestone |
|---|---|
| 11-21 | cluster discovery, tool installs (RFD4 CUDA-12 venv), theozyme conditioning recipe working; Milestone 1 scoring validation; 32-design RFD4 batches; steering handoff received; evaluation suite built |
| 11-21 → 22 | overnight: our unsteered control (250 backbones), multi-seed Boltz scoring, partial diffusion; the large cluster ran the 10,500-backbone steered/unsteered campaign |
| 11-22 morning | full evaluation of finalists, BioEmu analysis, literature correlation study (null result after selection-aware correction), eligibility infrastructure |
| 11-22 noon | final selection (merged top 80 → eligibility → ranking → top 20), submission before 12:00 UTC |

## 5. Outcome

The 20 submitted sequences — on 10 distinct backbones — won the hackathon. The final list, with generation method for each,
is in [06_RESULTS_AND_FUNNELS.md](06_RESULTS_AND_FUNNELS.md). Honest framing for reuse: the ranking is built on an
in-silico stand-in for the hosts' score plus a weak literature calibration; **no design was experimentally tested here**.

## 6. Ten things worth knowing before you touch anything

1. The default `.venv` of the RFD4 clone is broken (CUDA 13 torch vs a CUDA 12.4 driver); use `.venv-cu12`.
2. RFD4 spec conditioning is easy to get silently wrong — follow the recipe in `reference/TOOL_NOTES.md`.
3. Boltz-2 *without* an MSA often puts the ligand somewhere other than the designed pose; this "ligand pose mismatch" was the largest single loss in every funnel (see 06).
4. Count **lineages**, not designs, for steered runs: the 3 beam members of a lineage are near-identical (Cα RMSD ≈ 0.03 Å).
5. Pass criteria that equal the steering rewards are circular — passing is not evidence of catalytic competence.
6. tmol relaxes are stochastic (≈ 5 REU noise); PyRosetta Sc is unreliable.
7. Design ids were truncated to 50 characters in an early eligibility table and silently collided — join by **sequence**.
8. A long-running UniRef90 GPU search (needs ≈ 28–30 GB RAM → the 32 GB node) can take an hour; start it early.
9. The team's GPU quota is an AWS vCPU quota: one 8-vCPU node **or** two 4-vCPU nodes, never mixed.
10. Never send project sequences or structures to external services; all eligibility checks here use local databases.
