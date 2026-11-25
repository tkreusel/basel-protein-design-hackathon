# 07 Lessons and pitfalls

Each item cost real time. Grouped by where it bites.

## Environment and cluster
| Pitfall | Fix |
|---|---|
| RFD4's default `.venv` has CUDA-13 torch; the cluster driver is CUDA 12.4 → GPU invisible, silently CPU-only | use `.venv-cu12`; never rerun the repo's `setup.sh` |
| Test GPU code on a compute node; the head node has no GPU | everything GPU through Slurm |
| Team quota is an AWS vCPU quota (`VcpuLimitExceeded` → `down#` node, jobs stuck in `CF`) | one `l4-1` or both `l4small` pinned with `--nodelist`; never mix; an idle node still holds its vCPUs for a cooldown |
| `l4big` is outside the allocation (jobs queue forever) | `--exclude=gpu-dy-l4big-1` on every GPU job |
| `sbatch --mem=40G` unsatisfiable on the 32 GB node | UniRef90 search needs `--mem=30G` and node `l4-1` |
| `TMPDIR` on `/mnt/work` in GPU jobs → noisy `Device or resource busy` at multiprocessing teardown | node-local `TMPDIR` |
| A long UniRef90 search (≈ 1 h for ~45 sequences once) | start identity searches early; it is the slowest eligibility step |
| `pkill -f <name>` inside a command whose own text contains `<name>` kills your shell | use PIDs |

## RFD4 / generation
| Pitfall | Fix |
|---|---|
| Hydra syntax is `key=value`; `--key` is rejected; `.pdb` inputs rejected | `.cif`/`.bcif`; run from the repo root |
| Theozyme conditioning failed silently twice before working | follow `reference/TOOL_NOTES.md` recipe exactly |
| 20 of 32 early designs satisfied every catalytic distance with the ligand on the protein surface | always add a burial criterion and verify plausibility (BURIED conditioning fixed it: 12/32 → 31/32 in the full-motif spec) |
| 3 beam members of a steered replicate are one backbone | count lineages |
| Pass criteria identical to rewards are circular | report as screens, not evidence |
| Partial-diffusion parents without the `C_CRD` column lose the pinned motif | use RFD4 `*_generated.cif` only; refold CIFs cannot be parents |
| His-as-base motifs cannot be declared (header has no tautomer field) | design around it; the His-motif batch (job 58) gave 3/32 backbones and was not used for the final designs (see `reference/HIS_ZN_MOTIFS.md`) |

## Sequence design / refolding
| Pitfall | Fix |
|---|---|
| Boltz without MSA often puts the ligand elsewhere (largest funnel loss: 84 % of Cα-OK refolds fail ligand RMSD ≤ 2 Å) | accept it as the gate; the hosts fold *with* MSAs — we checked de novo designs fold identically |
| Boltz ipTM is saturated (0.85–0.94) | do not rank on it alone |
| Catalytic residues must be *verified* fixed after LigandMPNN | the job asserts it |
| `--diffusion_samples` share a trunk → overstate consistency | use independent seeds for consistency |

## Scoring
| Pitfall | Fix |
|---|---|
| Base O…N3 distance anti-correlates with activity in the HG3 series (weakest variant has the shortest) | do not reward "shorter"; the donor distance discriminates; regression test pins this |
| The nitro group is not the oxyanion contact (7.6–8.7 Å from the donor tip) | atom map: N3 = abstraction centre, N1 = oxyanion centre |
| `geom_scores` overwrote raw distances with 0–1 scores → inflated consistency | per-metric scores now `sc_<metric>` (fixed 2026-11-22 04:30 UTC) |
| tmol relax is stochastic (≈ 5 REU) | average replicates; treat < 5 REU as noise |
| PyRosetta `FastRelax` leaves a coordinate-restraint weight on the score function (+5,200 REU on every later score) | relax on a clone (fixed in the suite) |
| PyRosetta Sc unreliable (0.78 unrelaxed vs 0.38 relaxed for a 0.3 Å move) | treat as unvalidated; `dG_sep` is the usable term |
| PLACER analysis crashes when the cropped model lacks wall-ring atoms | `our_additions/geometry_wall_ring_fallback.patch` |
| A PLACER-based "donor median" correlated *positively* (+0.775) with activity, contradicting the criteria | excluded (outlier) |
| Several selection-free correlations looked significant | selection-aware permutation test removed all of them |

## Eligibility
| Pitfall | Fix |
|---|---|
| Foldseek's official PDB download was blocked (HTTP 429) | built our own DB from one representative per 30 % sequence cluster via `files.rcsb.org` |
| TM stage: `Path(cid).stem` truncated dotted ids (TM 0.0); NMR models flooded the top-20 hit list | use the id as is; deduplicate best hit per entry |
| Design ids truncated to 50 characters collided (`mpnn1/2/4`), giving wrong/missing eligibility for 5 of the top 30 | join by sequence; full unique ids |
| Relative structure paths crashed a long run after hours of search | absolute paths in `candidates.csv` |
| A control that is guaranteed to fail is the only proof the gate works | always include published HG3.17 (it must FAIL) |
| Six of the eight eligibility failures were over TM 0.70 by less than 0.01 | treat TM within ~0.02 of the limit as a risk; the hosts' alignment tool may differ slightly from ours |

## Campaign pipeline (large runs)
| Pitfall | Fix |
|---|---|
| The metric registry (`configs/metrics/custom_metrics/`) and the reward code are gitignored in the RFD4 clone, so a plain clone cannot resolve `kemp_pocket3` or `kemp_channel18` | `rfd4_patch/` carries both; `jobs/install_rfd4.sh` installs them |
| `6NT.tmol` (TMol ligand parameters) and the beta-sheet input specs were missing from the first hand-over | now in `pipelines/02_rfd4_steering/inputs/` |
| Single-GPU array jobs waited many hours on a shared queue while whole 8-GPU nodes started quickly | whole-node jobs with one worker per GPU; workers pull tasks and take locks (`jobs/worker.sh`) |
| GPUs were CPU-bound: one process per GPU left the GPU idle | two to three processes per GPU gave about twice the throughput |
| PLACER (torch 2.3, cu121) failed with "no kernel image" on newer GPUs | restrict it to GPUs that work (`PLACER_GPU_PATTERN`); other stage-B jobs hand it to a separate job |
| An Sc failure in PyRosetta discarded dG_sep, packstat and hydrogen-bond terms for that replicate | Sc is computed in a try block; `jobs/repair.sh pyrosetta` re-runs affected replicates |
| One unanalysable PLACER model cost the whole run its PLACER summary | per-model try block in `analyze_placer.py` |
| Unsteered backbones almost never pass the base-environment rules, so a single gate would remove them all | decision flag differs by arm family, both flag sets are stored, fallback tiers fill each run to a minimum and are marked |
| Extended-evaluation PLACER on non-PLACER GPUs was never picked up | `stage_P.sh --extended` and the worker loop cover both evaluation sets |
| Results live on unbacked scratch and the home directory caps the file count | `jobs/archive_runs.sh` packs one archive per run; structures can be regenerated from it |

## Process
- Write the definition of success next to each run folder (`SUCCESS.md`) before reading results.
- Join everything by sequence hash or full id, never by a display name.
- Keep the reproducibility record in the folder: run folders carry their summary, attempts and verdicts.
- Report outcomes faithfully: when a metric failed validation (Milestone 1 criterion 1, Spearman 0.50 at n = 3) say so.
