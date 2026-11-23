# 04 — Generation methods: conditioning, steering, partial diffusion

All three use the same model (**RFD4-Proteina**, `nvidia/RFD4-Proteina`) and the same active-site problem. They differ in what is
fixed at input, what is selected during sampling, and what starts the trajectory.

## Conditioning (every arm, incl. "unsteered")

Part of the **input specification**. We tell RFD4 the problem; the model samples backbones that respect it, with no further check.

| Conditioning | Meaning |
|---|---|
| Motif tips fixed, **unindexed** | Asp127 carboxylate + Gln50 side chain (+ Trp44 ring in the "with-Trp" motif) at fixed absolute coordinates relative to 6NT; the model decides where in the chain they sit |
| Ligand 6NT fixed | HG3.17 pose |
| `BURIED` (rel. SASA) on motif + ligand | raises base burial dramatically alone (base SASA 0.32 → 0.07 in round 2) |
| Hydrophobic contact (`C_HPHOB`) | on the base CB/CG and the ligand ring carbons |
| H-bond acceptor (`C_HB`) | on the nitro oxygens |

Specs: `pipelines/01_backbones_sequence_refold/inputs/*.yaml` and `pipelines/02_rfd4_steering/inputs/kemp4_*.yaml`
(`hg317-tips_with-trp`, `hg317-tips_no-trp`, `gdm1483-glu-ser`). Recipe and the failures that led to it: [reference/TOOL_NOTES.md](reference/TOOL_NOTES.md).

**Unsteered** = conditioning only. Quality comes from the model's own distribution plus a post-hoc design gate. In the campaign
only 1–3 % of unsteered backbones passed the strict `pass_core` gate, against 16–30 % in the steered arms.

## Steering (inference-time reward-guided beam search)

Same model, same conditioning; the **samples are selected during denoising** using external scores. RFD4's built-in `search:` block
(`rfd4/src/rfproteina/search/beam_search.py`):

- `n_branch = 4` branches are denoised K steps from each beam member; each branch's *predicted clean structure* is scored by a weighted
  sum of metrics (rewards); the top `beam_width = 3` are kept.
- Checkpoints at denoising steps **30, 42, 52** (predicted-clean scoring, no rollout) and **60** (final).
- Each replicate returns 3 beam members of one **lineage** — near-identical backbones (Cα RMSD ≈ 0.03 Å). Count lineages, not designs.
- It is **non-differentiable resampling only**; no gradients flow to the model and no weights change.

Custom reward metrics (`pipelines/02_rfd4_steering/rfd4_patch/custom_metrics/`; configs reference `kemp_base`, `kemp_pocket3`, `kemp_chem`,
`kemp_chem_hg`, `kemp_channel`, `tmol_lig_ddg`):

| Metric | Rewards / penalises |
|---|---|
| `kemp_base.py` | base burial (rel. SASA), apolar environment fraction, side-chain polar atoms ≤ 4.5 Å of the base O (weight −4), Lys/Arg ≤ 6 Å (−3), nitro H-bond donors far/near the base, ligand and protein clashes |
| `kemp_pocket.py` | grid pocket volume / ligand volume (band 0.4–1.4), ligand SASA fraction, loop fraction in the 8 Å shell |
| `kemp_chem.py` | locates the unindexed motif residues by superposing the reference theozyme on the ligand; base/partner RMSD (motif rigidity) |
| `kemp_channel.py` | probe-radius connectivity pocket → bulk; reward = bottleneck radius in a band (1.4–1.9 Å; 1.7–1.9 Å in `steer_tmol_ch18`) |
| `kemp_tmol.py` | **TMol** protein–ligand interface ΔG (needs a prebuilt `6NT.tmol`; weight 0.5/20) |

Arms: `unsteered` | `steer_base` (geometry rewards) | `steer_base_tmol` (+ TMol dG) | `steer_base_channel` (+ channel band) — round 4/5 —
and in the large overnight campaign `steer_tmol` | `steer_tmol_ch18` | `unsteered`, each for the `hgnt` (no Trp44) and `hgwt` (with Trp44) motifs.

**Design gate on the output** (`pass_geom`, `pass_seq`, `pass_core`): ligand/base SASA ≤ 0.10, pocket ratio 0.4–1.4, no clashes, motif RMSD ≤ 0.3 Å;
no polar atoms within 4.5 Å of the base O, no Lys/Arg within 6 Å, nitro H-bond away from the base, TMol dG < 0. Steered arms need `pass_core`; unsteered needs
`pass_geom` (the unsteered arm essentially never passes the base-environment rules).

**Caveats stated by the authors** (keep them when you report):
- *Circularity:* pass criteria equal the reward terms, so passing is not evidence of catalytic competence; TMol dG is the only partly independent term and is itself a reward in some arms.
- Effective sample size = lineages (3 beam members are not independent). No confidence intervals.
- Pocket volume and channel metrics are heavy-atom grid proxies.

## Partial diffusion (refinement of good designs)

Re-noise an existing RFD4 design to a chosen level and regenerate the tail of the trajectory (`generation.args.partial_trajectory_fraction_steps`;
0.3 = re-run the last 30 %, i.e. close to the parent). Our design (`pipelines/05_partial_diffusion/README.md`):

- Generation stays **unsteered**; the only selection is a CPU rejection filter on backbone properties (motif tips and ligand intact, ≥ the parent's number of core
  criteria, base burial ≥ parent) before any GPU refold. Survivors go through MPNN → Boltz → host-style scoring and only children that beat their parent are kept.
- Parents must be annotated RFD4 `*_generated.cif` files (they carry the `C_CRD` fixed-atom column); Boltz refolds cannot be parents.
- Run: 12 parents × 8 replicates at levels 0.3 and 0.5 → 192 children → 112 passed the filter → 64 selected → 256 sequences → 16 beat their parent.
- None of the children entered the final 20 (best child ranked 22 of the merged pool).

## Why this matters for the funnel
Steering spends compute on search, which is why far more steered backbones pass the design gate — but the real bottleneck is downstream
(see [06](06_RESULTS_AND_FUNNELS.md)): Boltz-2 without an MSA often does not reproduce the designed ligand pose.
