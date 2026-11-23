# 03 — Reproducing the pipelines

Prerequisites: a Slurm GPU cluster like ours, the environments in [02](02_ENVIRONMENTS_AND_CLUSTER.md), and the organisers' reference
paths. Everything below was run this way; seeds, thresholds and tool flags are the ones we used. **Check `squeue` before every GPU
submission** (team quota = one 8-vCPU node or two 4-vCPU nodes).

## 0. Workspace

```bash
git clone <this repo> && cd kemp-eliminase-pipelines
scripts/make_workspace.sh  $WORK            # symlink tree with the layout the scripts expect ($WORK = $KEMP_ROOT)
scripts/fetch_benchmark.sh $WORK            # literature crystal structures from RCSB (inbound only)
# create $WORK/env/kemp-orch as described in docs/02, then
cd $WORK && source env/paths.sh && export KEMP_ROOT=$WORK
PYTHONPATH=src env/kemp-orch/bin/python -m pytest tests/ -q          # expect 30 passed
```

The tests encode the real observed values (donor distance separates active from inactive; base distance does not). Run them
before trusting any change to the scoring code.

## 1. Foundation checks (CPU)

```bash
env/kemp-orch/bin/python -I src/derive_atom_map.py benchmark/pdb   # 6NT atom roles from 11 ligand copies (prints the table behind configs/atom_map.yaml)
env/kemp-orch/bin/python -I src/milestone1.py .               # geometry stack on the HG series (docs/reference/MILESTONE1.md)
env/kemp-orch/bin/python -I src/check_motif.py . inputs/HG3_17_kemp_motif.pdb
env/kemp-orch/bin/python -I src/motif_to_cif.py <motif.pdb> <out.cif>   # RFD4 needs .cif inputs
```

## 2. Backbones

### 2a. Conditioning-only (our early batches and the control)
The RFD4 specs live in `inputs/` (`kemp_rfd4.yaml`, `kemp_rfd4_full_buried.yaml`, `kemp_rfd4_tips_buried.yaml`, `kemp_rfd4_his.yaml`);
`full_buried` = whole motif + BURIED flags, `tips_buried` = tip atoms only + BURIED. The conditioning recipe (unindexed motif, fixed
absolute coordinates, `BURIED`, hydrophobic-contact and nitro-O acceptor flags) is documented in
[reference/TOOL_NOTES.md](reference/TOOL_NOTES.md) — read it before editing a spec; two of the three first attempts failed silently.

```bash
sbatch --chdir=$WORK --output=$WORK/logs/stage1-%j.log \
  --export=ALL,KEMP_SPEC=$WORK/inputs/kemp_rfd4_full_buried.yaml,KEMP_N=32,KEMP_BS=2,KEMP_SEED=201 jobs/stage1_rfd4.sbatch
env/kemp-orch/bin/python -I src/finalize_run.py . <raw out_dir> runs/<name>        # run folder with attempts/ + success/ (docs/reference/RUN_FOLDER_CONVENTION.md)
```
Stage-1 gate (ligand buried, base/donor/wall at calibrated distances): `configs/thresholds.yaml`.

### 2b. One-job unsteered control (the 250-backbone run)
Mirrors the steered run's settings (his with-Trp spec, `batch_size=4`, `compile_lengths=null`, `~metrics`):
```bash
sbatch --chdir=$WORK --export=ALL,KEMP_N=250,KEMP_K=100,KEMP_MODE=topk,KEMP_SEED=7,KEMP_BOLTZ_SEED=0 overnight/run_overnight.sbatch
```
Stages GEN → SELECT (top-K by number of the 7 computable core criteria) → MPNN (4 sequences + native; K R H D E N Q S T Y C forbidden
within 6 Å of the base carboxylate O; motif fixed; seed 111) → BOLTZ (no MSA, 6NT, 3 recycles, 1 sample, seed 0) → SCORE.
Idempotent: resubmit with `KEMP_RUN=<run dir>` to resume. See `pipelines/04_overnight_unsteered_control/README.md`.

### 2c. Reward-steered beam search (cluster-side campaign)
```bash
python code/make_r4.py ; python code/make_r5.py                  # experiment YAMLs for the arms
sbatch -J <name> -o <log> --export=ALL,N_REP=3 code/run_r5_job.sbatch <campaign> <seed> <arm> [<arm>...]
python code/eval_r5.py                                           # design gate + per-arm summary
```
Arms: `unsteered | steer_base | steer_base_tmol | steer_base_channel` (round 4/5) and `steer_tmol | steer_tmol_ch18 | unsteered` (overnight campaign).
Needs the RFD4 clone patched with `rfd4_patch/` (custom metrics into `rfd4/src/rfproteina/metrics/custom_metrics/`, `pyproject.diff`).
Details, caveats and the reward list: [04_GENERATION_METHODS.md](04_GENERATION_METHODS.md) and `pipelines/02_rfd4_steering/README_AGENT_HANDOFF.md`.

### 2d. Partial diffusion
```bash
sbatch --dependency=afterany:<job> --exclude=gpu-dy-l4big-1,gpu-dy-l4small-1,gpu-dy-l4small-2 \
  --export=ALL,KEMP_JOB61=<multi-seed run dir>,KEMP_JOB60=<overnight run dir>,KEMP_BOLTZ_STOP=<epoch>,KEMP_MS_STOP=<epoch> \
  partial/run_partial.sbatch
```
Parents must be RFD4 `*_generated.cif` (they carry the `C_CRD` fixed-atom column). Parameters used: 12 parents × 8 replicates × noise levels
0.3 and 0.5 (`partial_trajectory_fraction_steps`). Stage logic and filters: `pipelines/05_partial_diffusion/README.md`.

## 3. Sequences and refolds

```bash
sbatch --export=ALL,KEMP_RUN=<finalized stage-1 run>,KEMP_NSEQ=4 jobs/stage2_ligandmpnn.sbatch        # LigandMPNN, motif fixed, verified
sbatch --export=ALL,KEMP_RUN=<run>,KEMP_TAG=<tag> jobs/refold_boltz.sbatch                            # Boltz-2 refold + re-gate
```
Boltz-2: `--model boltz2 --diffusion_samples 1 --recycling_steps 3 --no_kernels`, YAML with `msa: empty` and the ligand as `ccd: 6NT`.
Strict gate used in the campaign: pLDDT ≥ 80, Cα RMSD ≤ 2 Å, ligand ipTM ≥ 0.6, ligand RMSD ≤ 2 Å (after Cα alignment), catalytic tips ≤ 1 Å.

## 4. Host-style scoring (multi-seed)

```bash
sbatch --export=ALL,KEMP_CANDS=<candidates.csv>,KEMP_OUT=<out>,KEMP_SEEDS="0 1 2 3 4" scoring/boltz_seeds.sbatch   # 5 independent seeds, separate runs
env/kemp-orch/bin/python -I scoring/host_score.py <candidates.csv> <seeds root> <out dir>
env/kemp-orch/bin/python -I scoring/rank_v2.py <out csv>          # geometry-multiplicative stand-in used for the final ranking
```
Candidate lists: `scoring/pool_v1/candidates.csv` (106) and `scoring/pool_overnight_top24/candidates.csv`.
Formulas: [05_SCORING_AND_RANKING.md](05_SCORING_AND_RANKING.md).

## 5. Evaluation suite (finalists only)

```bash
# 1) substrate (H5J) and product (6VP) Boltz refolds for the finalists (same settings as the 6NT refolds)
python pipelines/06_evaluation_suite/our_additions/stage_finalists.py     # edit FINALISTS first
sbatch pipelines/06_evaluation_suite/our_additions/arms_then_placer.sbatch   # Boltz arms -> prepare + ipSAE2 -> PLACER (GPU, ~40 min for 12 designs)
# 2) CPU parts on the head node (tmol, PyRosetta, aggregate): see suite/jobs/run_full_suite.sh (dry-run by default, --submit to run)
bash pipelines/06_evaluation_suite/suite/jobs/run_full_suite.sh --root <boltz outputs with 6NT/H5J/6VP per protein> --submit
```
Apply `our_additions/geometry_wall_ring_fallback.patch` to the suite if PLACER crops lack wall-ring atoms (otherwise the analysis crashes).
Suite docs: `pipelines/06_evaluation_suite/PIPELINE.md`.

## 6. Eligibility

```bash
python pipelines/08_eligibility/make_candidates_v6.py          # builds candidates.csv with ABSOLUTE structure paths (relative ones crashed a run)
sbatch pipelines/08_eligibility/uniref_search.sbatch           # UniRef90 identity, GPU, --mem=30G, node l4-1 (expect up to ~1 h)
python -I pipelines/08_eligibility/check_eligibility.py candidates.csv --out results --threads 6 [--reuse] --uniref-tsv results/uniref90_hits.tsv
```
Thresholds: identity ≥ 0.70 fails, TM > 0.70 fails, warn within 10 %. The published HG3.17 sequence/structure is always included as a
control and **must** FAIL (it does: identity 1.00, TM 0.983).

## 7. BioEmu and literature correlation (optional, not used for selection)

```bash
sbatch pipelines/09_bioemu/run_bioemu.slurm <system.yaml> <n_samples>     # bioemu2 sample-equilibrium
python pipelines/09_bioemu/analyze_ensembles.py ...                         # RMSF, ligand pose, approach geometry per sample
python pipelines/09_bioemu/literature/ensemble_metrics.py ...; python pipelines/09_bioemu/literature/correlate.py
```

## 8. Final selection and upload

```bash
python -I pipelines/11_selection_and_submission/merged80_pool.py            # merged top-80 pool (Maxi top 30 + our scored candidates, dedup by sequence)
python -I pipelines/11_selection_and_submission/rerank_merged80.py          # eligibility filter + judging x literature factor -> rerank_top20.csv, designs_rerank_top20.fasta
python3 pipelines/11_selection_and_submission/validate-submission.py designs.fasta
./pipelines/11_selection_and_submission/submit-designs.sh designs.fasta     # only with the team's go-ahead; re-uploads replace the file
python -I pipelines/11_selection_and_submission/export_final_structures.py  # structures of the final 20
```
Inputs these need (not all in the repo — see [08](08_PROVENANCE_AND_DATA_LOCATIONS.md)): `ranking_v2.csv`, `maxi_top30_scored.csv`, eligibility CSVs,
the suite's literature table `kemp_eval_dataset_full.csv` and the campaign's `evaluated_sequences_full_suite.csv`.

## What is deterministic and what is not

| Step | Seeded? |
|---|---|
| RFD4 generation | yes (`seed=`); steered runs use 7 (r4), 11, 12 (r5); our control 7 |
| LigandMPNN | seed 111 |
| Boltz-2 | `--seed`; consistency uses seeds 0–4 as separate runs |
| tmol / PyRosetta relax | stochastic, ≈ 5 REU noise between identical runs |
| PLACER | sampling; 20 models per ligand state |
| Eligibility | deterministic given the databases |
