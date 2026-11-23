# Proposed run on the single-chain Boltz refolds: for sign-off (nothing submitted)

## What is already done (CPU, head node, no job, no GPU)
* Boltz refolds in `boltz_outputs/single_chain/` registered as `<ID>_boltz` (29SB, 4BS0, 5RGA, 7K4Q, 8RD5; 1 sample each) next to the 5 crystal/NMR entries.
* ipSAE2 + site confidence computed for all 5 refolds (`work/ipsae/`).
* Catalytic-geometry quick look: all 5 refolds reproduce the site (base–N3 2.5–2.8 Å; Gln 2.8–3.0 Å in 4BS0/8RD5 refolds; no donor in HG3).

## Replicates to run: 16
crystals 2+2+2+2 copies, 29SB NMR capped at 3 of 10 models, plus 5 Boltz refolds. (`MAXREP=3`; say if you want all 10 NMR models.)

## Proposed jobs (all resumable, one at a time because of the 8-CPU team budget)
| # | Job | Where | Estimated time | Needs GPU? |
|---|---|---|---|---|
| 1 | PLACER holo + apo, 50 models each, 16 replicates | `jobs/10_placer.sbatch`, 1x L4, 4 CPUs | ~1 h (GPU 1-3 s/model per PLACER docs; **unmeasured here**) | **yes, the only one** |
| 2 | tmol protocol (relax, theozyme constraints, same-pocket TSA/substrate/product, apo, ablations) | CPU: ~7 min/replicate measured on 8 cores -> ~2 h | ~2 h | no (GPU would only speed relaxation; unmeasured) |
| 3 | PLACER ensemble analysis + tmol field on top-20 models | CPU | ~1 h | no |
| 4 | PyRosetta (relax, interface dG, Sc, packstat, BUNS, independent ddG, ablations) | CPU, 1 core/replicate (~7 min without, ~18 min with ablations) | ~1.5 h with 4 parallel processes | no |
| 5 | aggregate -> `work/report/` | CPU, seconds | | no |

So **GPU use is limited to job 1** (about an hour of one L4). Jobs 2-4 can run on the CPUs of one GPU-partition node without reserving the GPU (`--gres` omitted), or on the head node. Which do you prefer?

## Decisions I need from you
1. Sign off on job 1 (GPU, ~1 h)? When is the GPU free for it?
2. CPU jobs: head node, or sbatch without `--gres`?
3. Include the validation ladder 7K4X (27,000) and 7K4Z (126,000) (adds 4 replicates, ~25 min of CPU, no GPU) so a rank correlation becomes possible? Default: off.
4. Your theozyme CST: still missing; the run uses the crystal-derived constraints until you send it.
5. Replicates: NMR 3 of 10 (default) or all?

## What the results can and cannot tell you (see PIPELINE.md section 6)
* The 5 refolds are native sequences of structures Boltz may have seen: high confidence is not independent evidence. The informative test comes from mutants/designs.
* tmol relax is stochastic: two identical runs differed by several REU (DDE -5.2 vs -3.6). Differences below ~5 REU between entries are noise; the crystal copies (2 per entry) give a first noise estimate.
* No activity calibration exists; features are for ranking only.
