# Kemp eliminase design pipeline — implementation plan (v0)

Target executor: Claude Code on the hackathon GPU cluster.
Assume nothing about the scheduler, module system or filesystem layout — discover it first (§2.0)
and write down what you find.

---

## 0. Ground rules for the implementing agent

Read this section before writing any code.

1. **Do not invent CLI flags, config keys or API signatures.** Several tools here are weeks old
   (tmol, BioEmu-2, RFD4-Proteina). For every tool: clone it, read the README and `--help`,
   run the bundled example, and **write the real interface into `docs/TOOL_NOTES.md`** before
   wiring it into the pipeline. If the example does not run, record the error and move on —
   do not spend more than the timebox in §7.
2. **Fail loud, never silently.** No `try/except: pass`. A stage that cannot score a design
   writes `status: "failed"` plus the exception into that design's record and continues with
   the rest of the batch.
3. **Everything is resumable.** Each stage checks for its own output and skips designs that
   already have it. We will crash, run out of walltime and get preempted.
4. **No tool imports another tool's environment.** See §2 — one conda env per tool, orchestration
   shells out via `subprocess`. Do not try to make RFdiffusion2, Boltz-2 and tmol coexist.
5. **Milestone 1 (§4) gates everything.** Do not generate a single backbone until the scoring
   stack can rank the HG series correctly. If the score cannot tell HG3 from HG3.17, it cannot
   tell our designs apart either.

---

## 1. What we are building

A filter cascade that turns a transition-state geometry into ranked, compliant Kemp eliminase
designs.

| # | Stage | Tool | Runs on |
|---|-------|------|---------|
| 0 | Theozyme set (6 cells) | hand-built + RDKit | once |
| 1 | Backbone generation | RFdiffusion2 | once |
| 2 | Sequence design | LigandMPNN | per variant |
| 3 | Physics scoring | tmol | **all** variants |
| 4 | Active-site ensemble | PLACER | survivors |
| 5 | Fold + interface | Boltz-2 | top ~10% |
| 6 | Compliance gate | MMseqs2, Foldseek, TM-align | survivors |
| 7 | Apo ensemble | BioEmu-2 | finalists only |
| 8 | Maturation loop | loop over 2–5 | v2 |

**v0 scope (this document's primary goal):** stages 0–3, 5, 6 + the Milestone 1 validation.
**v1:** add stage 4 (PLACER).
**v2:** add stages 7 and 8.

Hard constraints from the hackathon brief, enforced in stage 6:
- design length ≤ 260 aa
- < 70% sequence identity to any published enzyme
- TM-score ≤ 0.7 to any published structure

---

## 2. Environment setup

### 2.0 Discover the machine first — 10 minutes, do not skip

Write the answers into `docs/CLUSTER.md` before installing anything:

- **Scheduler:** `sbatch`/`squeue` (SLURM), `qsub` (PBS/SGE), or none. If none, stages run as plain
  background processes with a worker pool.
- **GPUs:** `nvidia-smi` — how many, which model, how much VRAM. VRAM sets the Boltz-2 and
  BioEmu-2 batch sizes and may rule BioEmu-2 out entirely.
- **Filesystem:** `df -h` plus any quota command. Model weights total tens of GB. Find the big
  shared scratch and put everything there.
- **Modules / containers:** `module avail`, and whether `apptainer`, `singularity` or `docker`
  exists. If published images for these tools exist, prefer them over building envs from source.
- **Conda:** is `conda`/`mamba` present, or does `micromamba` need bootstrapping?
- **Internet from compute nodes.** Many clusters allow outbound traffic only from login nodes. If
  compute nodes are offline, **download every model weight up front from a login node** and point
  the caches at them. This is the single most common way a hackathon pipeline dies at 2am.

### 2.1 Environments

Create **one env per tool**. These have mutually incompatible torch/CUDA pins; do not merge them.

| Env name | Contains | Notes |
|---|---|---|
| `kemp-orch` | python, rdkit, numpy, pandas, biotite, pyarrow | orchestration only, no GPU deps |
| `kemp-rfd2` | RFdiffusion2 | check whether weights need separate download |
| `kemp-mpnn` | LigandMPNN | small, fast to install |
| `kemp-tmol` | `pip install tmol` | pin the version you get and record it |
| `kemp-placer` | github.com/baker-laboratory/PLACER | weights are in the repo per its README |
| `kemp-boltz` | `pip install boltz` | needs a weights cache dir — set it on scratch, not $HOME |
| `kemp-bioemu` | github.com/microsoft/bioemu2 | v2 only; Linux + NVIDIA GPU required |

Put every model weight cache on the shared scratch filesystem, never in `$HOME` — home quotas are
usually small and you will hit them mid-run. Export the cache variables in `env/paths.sh` and
source that file at the top of every job script.

Record for each env: exact package versions, CUDA version, and whether the bundled example ran.

---

## 3. Repo layout and the data contract

```
kemp/
  env/paths.sh                 # scratch dirs, weight caches, module loads
  configs/
    theozyme_cells.yaml        # the 6 cells (§5)
    thresholds.yaml            # all gate cut-offs, one place
    atom_map.yaml              # 6NT <-> TS <-> substrate atom correspondence (§6.1)
  ligands/
    6NT.sdf  substrate.sdf  product.sdf  ts.sdf
  src/kemp/
    record.py                  # the DesignRecord schema, read/write
    stage0_theozyme.py
    stage1_backbone.py
    stage2_sequence.py
    stage3_tmol.py
    stage4_placer.py
    stage5_boltz.py
    stage6_compliance.py
    stage7_bioemu.py
    geometry.py                # all geometric metrics, pure functions, unit-tested
    gate.py                    # applies thresholds.yaml, returns pass/fail + reasons
  jobs/                        # one job script per stage (sbatch or plain shell)
  runs/<run_id>/designs/<design_id>/   # all artifacts for one design
  runs/<run_id>/manifest.jsonl         # one DesignRecord per line
  docs/TOOL_NOTES.md
  tests/
```

### DesignRecord

One JSON object per design, appended to `manifest.jsonl`, also written as
`runs/<run_id>/designs/<design_id>/record.json`. **Every stage reads and writes this and nothing else.**

```python
{
  "design_id": str,            # stable hash of (cell, backbone_seed, seq_seed)
  "run_id": str,
  "cell": {"donor": "a|b|c", "wall": "1|2", "rasa": float},
  "provenance": {"stage1_seed": int, "stage2_seed": int, "parent_id": str|None, "generation": int},
  "paths": {"backbone_pdb": str, "design_pdb": str, "fasta": str, ...},
  "sequence": str,
  "length": int,
  "catalytic_residues": {"base": int, "donor": int|None, "wall": [int, ...]},  # 1-indexed
  "scores": { ... },           # flat float dict, namespaced: "tmol.ddE_ts_sub", "placer.geom_frac"
  "gates": {"stage3": "pass|fail", "stage3_reasons": [str], ...},
  "status": "ok|failed",
  "errors": [str]
}
```

Rules: scores are **append-only** and namespaced by stage. Never overwrite another stage's keys.
Never store large arrays in the record — write them to the design dir and store the path.

---

## 4. MILESTONE 1 — validate the scoring stack before generating anything

This is the most important part of the plan. We have a published series of Kemp eliminases with
measured activities spanning three orders of magnitude, several with crystal structures that
contain the transition-state analogue **6NT**. If our composite score cannot reproduce that
ordering, it is noise.

### 4.1 Benchmark set

Structures (fetch from RCSB, all are the same xylanase scaffold):

| PDB | Variant | kcat/KM (M⁻¹s⁻¹) |
|---|---|---|
| 7K4Q | HG3 | 146 |
| 7K4X | HG3.7 | 27,000 |
| 7K4Z | HG3.17 | 126,000 |
| 8RD5 | HG3.R5 | — |
| 6SS1 | HG3.17 Q50A / E47N / N300D | — (oxyanion ablated) |
| 6SS3 | HG3.17 Q50K / E47N / N300D | — (oxyanion ablated) |

Additional activity points without 6NT structures, usable via the predict-from-sequence path:
HG3.3b 2,200 and HG3.14 52,000. Sequences: start from the HG3 sequence and apply the mutation
lists in Broom et al. 2020 (Nat Commun 11:4808) Table 1. **Verify each reconstructed sequence
length and mutation position against the deposited HG3 sequence before using it.**

Negative controls, constructed in silico from 7K4Z:
- `HG3.17_D127A` — catalytic base ablated. Must score badly.
- `HG3.17_scramble` — pocket residues randomised. Must score badly.

### 4.2 Acceptance criteria

Run stages 3, 5 and (in v1) 4 on the benchmark set and compute:

1. **Spearman correlation** between the composite score and log10(kcat/KM) across HG3, HG3.7,
   HG3.17 (+ HG3.3b, HG3.14 via sequence). Target: ρ ≥ 0.7. Report it honestly; if it is near
   zero, the individual metrics are wrong, not the correlation.
2. **Ablation separation:** HG3.17 must outscore `HG3.17_D127A` and `HG3.17_scramble` on the
   composite, and must outscore 6SS1/6SS3 specifically on the **donor–oxyanion** metric.
3. **Known-geometry recovery:** the geometric metrics computed on 7K4Z must match values measured
   directly from the deposited coordinates. This is a unit test, not a model test.

Write the result to `docs/MILESTONE1.md` as a table plus the scatter plot. If criterion 1 fails,
fix the metric definitions before generating. If criterion 2 fails, the donor metric is broken.

### 4.3 Calibrate the thresholds here

Do not guess cut-offs. Measure the geometric values on 7K4Z and 7K4X and set
`configs/thresholds.yaml` from the observed distributions. The only two cut-offs we take from
literature are RFdiffusion2's: catalytic heavy-atom RMSD < 1.5 Å and no protein–ligand atom pair
closer than 1.5 Å.

---

## 5. Stage 0 — ligands and the theozyme set

### 5.1 Ligands

Build four SDF files with explicit hydrogens and 3D coordinates:

- `6NT.sdf` — **6-nitrobenzotriazole**, PDB chemical component ID `6NT`, formula C6H4N4O2.
  Get the ideal coordinates from the RCSB CCD entry for `6NT`; do not generate from a SMILES
  you typed from memory. Cross-check the formula.
- `substrate.sdf` — 5-nitrobenzisoxazole, the real substrate.
- `product.sdf` — 2-cyano-4-nitrophenolate (anionic form).
- `ts.sdf` — transition-state geometry. v0 approach: take the 6NT pose from 7K4Z, convert the
  ring to the isoxazole connectivity, and lengthen the C3–H and N–O bonds toward the TS.
  **This is a coordinate template, not a chemically valid molecule, and that is fine** — it is
  only ever used to supply atom positions to RFdiffusion2's atomised motif conditioning (§6.2).
  If xtb or ORCA is available on the cluster, do a proper constrained TS optimisation instead and
  note it in TOOL_NOTES.

Sanity check every SDF with RDKit: correct formula, no valence errors (expected for `ts.sdf`
— record it as known), reasonable bond lengths.

### 5.2 The six cells

`configs/theozyme_cells.yaml`. Asp base throughout (not Glu).

| Cell | Oxyanion donor | Pocket wall |
|---|---|---|
| a1 | Gln/Asn amide → phenolate oxygen | two-sided sandwich |
| a2 | Gln/Asn amide → phenolate oxygen | single aromatic face |
| b1 | Tyr/Ser/Thr hydroxyl → nitro group | two-sided sandwich |
| b2 | Tyr/Ser/Thr hydroxyl → nitro group | single aromatic face |
| c1 | none (water-mediated) | two-sided sandwich |
| c2 | none (water-mediated) | single aromatic face |

Cross each cell with a RASA burial sweep over the catalytic carboxylate: generate at
`{0.0, 0.1, 0.2, 0.3}` or whatever the RFdiffusion2 RASA conditioning actually accepts —
**check the interface and record it**. Burial is a window to hit, not a value to maximise.

Generate equal numbers per cell. Do not concentrate on the cell that looks best.

---

## 6. Stage specifications

### 6.1 `atom_map.yaml` — do this first, everything depends on it

Superpose `6NT.sdf` onto the ligand in 7K4Z, then onto `ts.sdf` and `substrate.sdf`, and write
out the explicit atom correspondence for the three positions the metrics need:

- `abstraction_centre` — 6NT's acidic N–H nitrogen; the substrate/TS C3
- `oxyanion_centre` — the atom bearing developing negative charge
- `nitro_group` — the nitro N and both O

Derive this from the actual coordinates and record which atom names map to which. **Do not
hard-code atom names into `geometry.py`; read them from this config.** Getting this mapping
wrong silently invalidates every geometric score downstream.

### 6.2 Stage 1 — backbone (RFdiffusion2)

Input: theozyme cell + `ts.sdf` coordinates.
Output: backbone PDBs.

- Use **atomised, unindexed motif** conditioning: supply only the functional-group atoms
  (carboxylate OD1/OD2, donor heavy atom, aromatic ring atoms), not full residues with fixed
  sequence positions. RFdiffusion2's own ablation found inferring rotamers and indices during
  sampling beat supplying native values.
- Use the **RASA** conditioning on ligand atoms to sweep burial, and the **ORI** token to place
  the site relative to the scaffold core.
- Target length: sample uniformly in [150, 255] so stage 6's 260 limit is never binding.
- Volume: start with 100 backbones per cell (600 total) to shake out the pipeline; scale after.

Acceptance: ≥ 80% of outputs are parseable PDBs with the motif atoms present and the catalytic
heavy-atom RMSD to the input motif under 1.5 Å.

### 6.3 Stage 2 — sequence (LigandMPNN)

Input: backbone PDB + `6NT.sdf` as ligand context (**not** `ts.sdf` — LigandMPNN needs a real
molecule).
Output: sequences + packed side chains.

- **Fix the catalytic residues.** Pass the base, donor and wall positions in the fixed-residue
  list. Verify in the output that they were not redesigned — assert this, do not trust it.
- Use a soluble-weighted model variant if available.
- 8 sequences per backbone, sampling temperature from the repo default; record it.

Acceptance: catalytic residues unchanged in 100% of outputs; length ≤ 260.

### 6.4 Stage 3 — physics (tmol), runs on everything

Input: designed PDB + the four ligands.
Output: `scores["tmol.*"]`.

Implement:
- `tmol.relax_energy` — batched relax, total RosettaGenFF energy, normalised per residue.
- `tmol.dE_ts`, `tmol.dE_sub`, `tmol.dE_prod` — interaction energy for each ligand placed in the
  same pocket. The TS is held rigid at its designed pose; the substrate and product are repacked.
- **`tmol.ddE_ts_sub` = `dE_ts` − `dE_sub`** — the ranking metric. Large and negative is good.
  Near zero means we built a binder, not a catalyst.
- `tmol.product_release` = `dE_prod`, penalised if strongly negative (product inhibition).
- `tmol.base_sasa` — SASA of the catalytic carboxylate (CG, OD1, OD2), reported as relative SASA
  against an Asp in a Gly-Asp-Gly reference.

tmol parameterises ligands automatically in ~1 s; if that fails on `ts.sdf` because of the
non-integer bond orders, fall back to assigning integer bond orders to the TS template and record
the choice. Batch via `PoseStack`.

Acceptance: runs on the Milestone 1 set, and `ddE_ts_sub` separates HG3.17 from
`HG3.17_scramble`.

### 6.5 Stage 4 — active-site ensemble (PLACER) — **v1**

Input: designed PDB + sequence + `6NT.sdf`. PLACER takes the backbone as input rather than
predicting it, which is exactly what we have here.
Output: ensemble of ~30–50 models, and:

- `placer.geom_frac` — fraction of ensemble members satisfying **both** the base–abstraction-centre
  distance and the donor–oxyanion distance within the thresholds calibrated in §4.3.
- `placer.pose_spread` — RMSD spread of the ligand across the ensemble.
- `placer.donor_chi_spread` — circular standard deviation of the donor side-chain torsions. This
  is the Gln50-rigidity question measured directly; low spread is good.

Budget: 1–3 s per model, so ~1 min per variant at 30 models. **Check early whether PLACER batches
on the GPU** — this single fact sets the population size for stage 8.

### 6.6 Stage 5 — fold and interface (Boltz-2), survivors only

Input: sequence + `6NT.sdf` (and a `substrate.sdf` control arm).
Output:

- `boltz.ipsae_p2l` — ipSAE in the **protein→ligand** direction: align on protein residues, score
  all ligand heavy atoms, **no PAE cutoff**, fixed d0 = 2 Å. Implement from
  github.com/dunbracklab/IPSAE; do not use plain ipTM for ligands, where a single atom can
  dominate the score.
- `boltz.plddt_global`, `boltz.pae_interface`
- `boltz.seed_consistency` — ligand pose RMSD spread across 4 seeds.
- `boltz.ablation_gap` — the same design's composite score minus its base→Ala mutant's.
  **Rank on gaps, not absolutes**; absolute confidence is not calibrated on de novo folds.

Do **not** use the Boltz-2 affinity head for ranking. It is trained on inhibitor IC50 data and
tighter binding is not the objective.

### 6.7 Stage 6 — compliance gate

Three hard checks, run before any final ranking:

1. `len(sequence) <= 260`
2. MMseqs2 or BLASTp against UniProt + PDB seqres. Record **both** global identity and the
   identity over the aligned region plus coverage. A 95%-identical 100-residue stretch inside a
   250-residue design fails the spirit of the rule even if global identity passes — flag it.
3. Foldseek search against the PDB, then TM-align or US-align on the top 20 hits for exact
   TM-scores. **TM-score is asymmetric** — compute both normalisations and report the higher.

Output a one-line pass/fail per design with the numbers, into `scores["compliance.*"]`.

### 6.8 Stage 7 — apo ensemble (BioEmu-2) — **v2, finalists only**

Sample the equilibrium ensemble of the top ~10 designs in both apo and holo form. Metric:
fraction of the apo ensemble in the catalytically productive conformation, using the same
geometric criteria as stage 4.

This is the HG3-vs-HG3.17 question: unbound HG3 through HG3.14 had the pocket residues blocking
productive binding, while HG3.17's unbound state is productive in ~40% of crystal molecules.
Needs Linux + NVIDIA GPU, and 4+ GPUs for context-parallel sampling.

### 6.9 Stage 8 — maturation loop — **v2**

Wraps stages 2–5.

- Population 50 per generation, seeded from the top 20 v0 designs.
- Moves: LigandMPNN resampling at temperature (biased toward the parent sequence to limit
  mutational load), plus RFdiffusion2 partial diffusion for backbone moves.
- Tier order **by cost**: tmol on all 50 → PLACER on the top ~20 → Boltz-2 at 4 seeds on the top ~8.
  Roughly 35 min per generation on one GPU.
- Selection: **Boltzmann selection** with the temperature lowered each generation (exploration →
  exploitation), and **forward-thinking selection** — a variant's potential is the average score
  of itself *and its first-generation descendants*. Both taken from AI.zymes.
- Log every generation's fitness and each component separately. The fitness-vs-generation plot is
  the deliverable even if activities come back dead.

---

## 7. Running jobs, and timeboxing

- One job script per stage in `jobs/`, all sourcing `env/paths.sh`.
- **Batch, never one job per design.** Chunk to ~50 designs per task so model load time is
  amortised. With SLURM use job arrays; without a scheduler, a fixed worker pool over the same
  chunks (one worker per GPU).
- Send logs to `runs/<run_id>/logs/` and keep them — the failure modes here are all in the logs.
- Get partition, QOS and GPU-request syntax from `sinfo -s` or the organisers. Do not copy
  resource flags out of examples; wrong flags usually mean the job silently runs on CPU.
- Make each stage idempotent and restartable (§0.3). Assume jobs get killed.
- **Timebox every new-tool spike to 2 hours.** If a tool will not run, record the failure in
  TOOL_NOTES.md and switch to the documented fallback below.

Fallbacks, in order of preference:
| If this fails | Do this |
|---|---|
| tmol ligand parameterisation | geometry-only gate; drop ΔΔE |
| PLACER | Boltz-2 seed spread as the preorganisation proxy |
| RFdiffusion2 | PNovo (github.com/google-deepmind/alphaprotein-novo) — it ships a Kemp eliminase example campaign, `GDM_KE_1483` |
| BioEmu-2 | drop stage 7; it is the sacrificeable phase |

---

## 8. Milestones and acceptance

| M | Deliverable | Acceptance |
|---|---|---|
| 0 | Envs built, every tool's bundled example runs, TOOL_NOTES.md written | each tool either runs or is documented as unavailable with its fallback chosen |
| 1 | **Scoring stack validated on the HG series** (§4) | Spearman ρ ≥ 0.7 vs log10(kcat/KM); ablations separate |
| 2 | Stages 0–2 produce 600 designs across 6 cells | ≥ 80% parse, catalytic residues intact, all ≤ 260 aa |
| 3 | Stages 3, 5, 6 run end-to-end; ranked compliant list | ≥ 20 designs pass all gates |
| 4 | PLACER integrated (stage 4) | metrics reproduce the §4.2 ablation separation |
| 5 | Maturation loop, 3 generations | fitness trajectory plotted per component |
| 6 | BioEmu-2 on finalists | apo productive fraction for the top 10 |

**Submission set** — pick as a designed experiment, not the top five rows of one sorted table:
1. best composite
2. best by `tmol.ddE_ts_sub` (often a different design)
3. best from a different oxyanion-donor cell
4. best from the other end of the RASA sweep
5. best from the maturation loop, or a second fold family
6. *(if six are allowed)* the base→Ala ablation of design 1 — a carboxylate in any hydrophobic
   pocket accelerates this reaction, so without the ablation a hit is uninterpretable

---

## 9. Open questions to resolve empirically — record answers in TOOL_NOTES.md

1. Does RFdiffusion2's RASA conditioning take per-atom values, and in what range?
2. Does PLACER batch on the GPU, or is it one model per forward pass? (Sets population size.)
3. Can tmol parameterise the TS template, or must bond orders be integerised?
4. Does Boltz-2 accept an arbitrary SDF ligand, or does it require a CCD code?
5. What is the actual wall-clock of one Boltz-2 prediction on this cluster's GPU? The ~20 s figure
   is from the paper's hardware and sets the whole stage-8 budget.
6. Does 6NT appear in the Boltz-2 / AF3 training cut-off window? If yes, the substrate control arm
   in stage 5 is load-bearing, not optional.

---

## 10. Reference numbers to calibrate against

| Enzyme | kcat/KM (M⁻¹s⁻¹) | Note |
|---|---|---|
| designed, pre-evolution (2008 series) | 6–160 | 8 actives of 59 designs |
| HG3 | 146 | no oxyanion donor |
| Des27 (fully computational) | 130 | kcat 0.07 s⁻¹, KM 0.5 mM, Tm 79 °C |
| HG3.3b | 2,200 | |
| Des27.7 | 12,700 | best fully computational, 141 mutations from any natural protein |
| HG3.7 | 27,000 | Gln50 introduced |
| HG3.14 | 52,000 | |
| HG4 | 103,000 | ensemble-based design, 8 mutations, no screening |
| HG3.17 | 126,000 | 17 rounds of directed evolution |

Realistic target for five first-shot designs with no wet-lab iteration: **any measurable turnover
above background**, stretch 10²–10⁴ M⁻¹s⁻¹. Published first-round hit rates are 4–14%.
