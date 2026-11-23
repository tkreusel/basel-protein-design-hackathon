# 04 Generation methods: conditioning, steering, partial diffusion

All three use RFD4-Proteina (`nvidia/RFD4-Proteina`) on the same active-site problem. They differ in what is fixed in the input, what is selected during sampling, and what starts the trajectory.

## Conditioning (every arm, including "unsteered")

Conditioning is part of the input specification. The model samples backbones that respect it and nothing checks them afterwards.

| Conditioning | Meaning |
|---|---|
| Motif tips fixed and unindexed | Asp127 carboxylate and Gln50 side chain (plus the Trp44 ring in the with-Trp motif) at fixed coordinates relative to 6NT; the model decides where in the chain they sit |
| Ligand 6NT fixed | HG3.17 pose |
| `BURIED` (relative SASA) on motif and ligand | on its own it lowered the base SASA from 0.32 to 0.07 in round 2 |
| Hydrophobic contact (`C_HPHOB`) | on the base CB and CG and the ligand ring carbons |
| H-bond acceptor (`C_HB`) | on the nitro oxygens |
| `C_SEC: SHEET` | only in the `*_sheet` diversity variants |

The specs are `pipelines/02_rfd4_steering/inputs/kemp4_*.yaml` (campaign) and `pipelines/01_backbones_sequence_refold/inputs/` (early rounds). The recipe, and the failures that led to it,
are in [reference/TOOL_NOTES.md](reference/TOOL_NOTES.md).

Unsteered generation is conditioning plus a post-hoc gate. In the campaign 1.2 to 2.8 % of unsteered backbones passed the full core gate, against 13 to 25 % of steered ones.

## Steering: inference-time beam search on rewards

Same model and conditioning, but samples are selected during denoising with external scores. RFD4's `search:` block (`rfd4/src/rfproteina/search/beam_search.py`) denoises
`n_branch = 4` branches for several steps from each beam member, scores the predicted clean structure of each branch by a weighted sum of metrics, and keeps the best `beam_width = 3`.
Checkpoints are at denoising steps 30, 42 and 52 (scoring the predicted clean structure) and 60 (the final sample). A replicate therefore returns three members of one lineage, near-identical backbones
(CA RMSD about 0.03 A). It is resampling only: no gradients reach the model and no weights change.

The reward terms (`pipelines/02_rfd4_steering/config/arms.yaml`; metric code in `rfd4_patch/custom_metrics/`):

| Metric module | Terms |
|---|---|
| `kemp_base` | base burial (relative SASA), apolar fraction of the environment, side-chain polar atoms within 4.5 A of the base O (weight -4), Lys or Arg within 6 A (-3), nitro H-bond donors far from or near the base, ligand and protein clashes |
| `kemp_pocket` (`kemp_pocket3`) | grid pocket volume over ligand volume with a band of 0.5 +- 0.1, ligand SASA fraction |
| `kemp_chem` (`kemp_chem_hg`) | locates the unindexed motif residues by superposing the theozyme on the ligand; base and partner RMSD as motif rigidity |
| `kemp_channel` (`kemp_channel18`) | probe-radius connectivity from pocket to bulk; reward is the bottleneck radius within a band (1.7 to 1.9 A) |
| `kemp_tmol` | TMol protein-ligand interface energy (needs `inputs/6NT.tmol`; weight 0.5 per 20 REU) |

Campaign arms: `unsteered`, `steer_tmol` (all terms except the channel band) and `steer_tmol_ch18` (adds the band and halves the ligand-SASA weight), each for `hgnt` and `hgwt`.
The earlier rounds 4 and 5 used `steer_base`, `steer_base_tmol` and `steer_base_channel` (code in `archive/early_rounds_r4_r5`).

Design gate on the output (`filter_pre.py`, thresholds in `config/gates.yaml`):

| Flag | Criteria |
|---|---|
| `pass_geom` | ligand SASA fraction and base SASA at most 0.10, pocket ratio 0.4 to 1.4, no clashes, base and partner motif RMSD at most 0.3 A |
| `pass_seq` | no side-chain polar atoms within 4.5 A of the base O (one allowed for with-Trp: Trp44 NE1), no Lys or Arg within 6 A, nitro H-bond away from the base, TMol dG below 0 |
| `pass_core` | both |

Steered arms are decided on `pass_core` and unsteered arms on `pass_geom`, because unsteered backbones almost never satisfy the base-environment rules; both flag sets are stored for every
design. If fewer than 30 designs of a run pass, the best remaining ones are added (`gate_tier` 2 if they pass `pass_geom`, 3 otherwise).

Caveats:

- Circularity: the pass criteria equal the reward terms, so passing is not evidence of catalytic competence. TMol dG is the only partly independent term and is itself a reward in the steered arms.
- The effective sample size is the number of lineages, not designs, and no confidence intervals were computed.
- Pocket volume and channel metrics are heavy-atom grid proxies.

## Partial diffusion

Re-noising an existing RFD4 design to a chosen level and regenerating the tail of the trajectory (`generation.args.partial_trajectory_fraction_steps`; 0.3 re-runs the last 30 %, staying close to the parent).
Generation stays unsteered; the only selection is a CPU rejection filter on backbone properties before any refold (motif tips and ligand intact, at least as many core criteria as the parent, base burial at least the parent's).
Parents must be annotated RFD4 `*_generated.cif` files, which carry the `C_CRD` column; Boltz refolds cannot serve as parents.
Run: 12 parents x 8 replicates at two noise levels gave 192 children, 112 passed the filter, 64 were selected, 256 sequences were designed and 16 beat their parent. None of them entered the final 20
(the best child ranked 22 of the merged pool); details in `pipelines/05_partial_diffusion/README.md`.

## Where the funnel loses

Steering spends compute on the search and therefore many more steered backbones pass the design gate. The bottleneck is downstream: Boltz-2 without an MSA often does not reproduce the designed ligand pose
([06](06_RESULTS_AND_FUNNELS.md)).
