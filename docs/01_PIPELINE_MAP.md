# 01 — Pipeline map

Each stage: what goes in, what comes out, which script, which environment, what it costs. Directory names refer to this repo.

```
              ┌──────────────────────────── theozyme + ligand (core/, pipelines/01/inputs) ───────────────────────────┐
              ▼                                                                                                       │
 S1  RFD4 backbones (GPU)  ──►  design gate (CPU)  ──►  S2 LigandMPNN (GPU/CPU)  ──►  S3 Boltz-2 refold, 6NT, 1x (GPU)
     conditioning / steering /       burial, motif RMSD,                motif fixed, KRHDENQSTYC                 CA ≤ 2 Å, ligand RMSD ≤ 2 Å,
     partial diffusion               base environment, clash            forbidden near base                      motif ≤ 1 Å, ipTM / pLDDT
                                                                                                                  │
                                         ┌────────────────────────────────────────────────────────────────────────┘
                                         ▼
 S4  deeper evaluation of survivors                                                    S5  selection
     • 5-seed Boltz refolds → host-style stand-in score (07)                             • merged top-80 pool (11)
     • substrate (H5J) + product (6VP) refolds → ipSAE2 (06)                             • eligibility filter (08)
     • PLACER holo/apo, tmol, PyRosetta, extra features (06)                             • judging × literature factor (05, 11)
     • BioEmu ensembles, literature correlation (09)                                     • top 20 → FASTA → validator → upload
```

| Stage | Purpose | Code (this repo) | Env / resources | Notes |
|---|---|---|---|---|
| **Theozyme & geometry foundation** | atom map, calibrated thresholds, geometry/burial scoring, tests | `core/src/kemp/`, `core/src/check_motif.py`, `check_burial.py`, `derive_atom_map.py`, `milestone1.py`, `core/configs/*.yaml`, `core/tests/` | `kemp-orch` venv, CPU | the scoring core; validated on the HG3 series ([reference/MILESTONE1.md](reference/MILESTONE1.md)) |
| **S1a RFD4 conditioning-only** | backbones around the fixed motif + ligand | `pipelines/01_backbones_sequence_refold/jobs/stage1_rfd4.sbatch`, `…/inputs/kemp_rfd4*.yaml`; unsteered control: `pipelines/04_overnight_unsteered_control/` | RFD4 `.venv-cu12`, 1 L4 GPU, ≈ 87 s/design incl. model load | recipe: [reference/TOOL_NOTES.md](reference/TOOL_NOTES.md) |
| **S1b RFD4 steering** | reward-guided beam search (4 branches × beam 3) | `pipelines/02_rfd4_steering/code/{make_r4,make_r5,run_r5}.py`, `rfd4_patch/custom_metrics/*.py`, `rfd4_patch/configs/` | RFD4 venv + tmol; cluster-side (not our GPU quota) | see [04](04_GENERATION_METHODS.md) |
| **S1c Partial diffusion** | re-noise good designs, regenerate the tail | `pipelines/05_partial_diffusion/` | 1 GPU (one chained job, `STAGES.log` records stage times) | none of its children reached the final 20 |
| **Design gate** | cheap CPU reject before any GPU refold | steered: `pipelines/02_rfd4_steering/code/eval_r5.py`; unsteered: `04…/select_designs.py` | CPU | pass_geom / pass_seq / pass_core |
| **S2 LigandMPNN** | sequences for each backbone | `core/src/run_ligandmpnn.py`, `run_stage2_ligandmpnn.py`, `04…/mpnn_batch.py`, jobs `stage2_ligandmpnn.sbatch` | `env-ligandmpnn`; T 0.1, 4 sequences, seed 111 | native RFD4 sequence + 4 redesigns (+ 1 shell-fixed for steered) |
| **S3 Boltz-2 refold** | does the sequence refold to the design with the ligand in place? | `core/src/prepare_refold.py`, `finalize_refold.py`, `04…/boltz_batch.py`, jobs `refold_*.sbatch` | `env-boltz2`, no MSA, 3 recycles, 1 sample | the dominant loss stage |
| **S4a Host-style scoring** | stand-in for the hosts' ranking | `pipelines/07_host_style_scoring/` | CPU + multi-seed Boltz (GPU) | geometry² · iPDE · burial · consistency |
| **S4b Evaluation suite** | ipSAE2, PLACER, tmol, PyRosetta, extra features | `pipelines/06_evaluation_suite/suite/`, `our_additions/` | 4 envs; PLACER on GPU, tmol/PyRosetta on head-node CPU | finalists only (≈ 40 min GPU for 12 designs × 3 ligands) |
| **S4c Eligibility** | identity < 0.70, TM ≤ 0.70 | `pipelines/08_eligibility/` | MMseqs2 (CPU + GPU UniRef90), Foldseek, US-align | local only; control (published HG3.17) must FAIL |
| **S4d BioEmu & literature** | ensemble behaviour; literature correlation | `pipelines/09_bioemu/` | BioEmu-2 on GPU; analysis CPU | null result after selection-aware correction |
| **S5 Selection** | pool → filter → rank → FASTA | `pipelines/11_selection_and_submission/` | CPU | `validate-submission.py`, `submit-designs.sh` |

## Hand-offs between stages (file contracts)

| From → to | Artifact |
|---|---|
| S1 → gate | RFD4 `*_generated.cif` (protein + ligand; conditioning columns incl. `C_CRD` kept — needed as partial-diffusion parents) |
| gate → S2 | `selection.csv` + backbone PDBs |
| S2 → S3 | `sequences.csv` (`seq_id, sequence, base_pos, donor_pos`) → Boltz YAML inputs (`msa: empty`, ligand `ccd: 6NT`) |
| S3 → S4 | `boltz_out/**/predictions/<id>/<id>_model_0.cif` + `confidence_*.json` + `pae_*.npz` |
| S4 → S5 | per-candidate tables: `ranking_v2.csv` (stand-in score), `eligibility_*.csv`, evaluation-suite tables |
| S5 → upload | `designs.fasta` with `base=`/`oxy=` headers, checked by `validate-submission.py` |

## Resource model (what constrained us)

- One Slurm partition (`gpu`), dynamic AWS nodes. Team budget = **8 vCPU of GPU workers at a time** (one `l4-1` or both `l4small`).
- `l4big` is outside the allocation. MMseqs2 against UniRef90 needs ≈ 28–30 GB RAM → only the 32 GB node fits.
- So GPU stages were chained, never parallel; CPU-only stages (tmol, PyRosetta, scoring, eligibility structure search) ran on the head node.
- Details: [reference/CLUSTER.md](reference/CLUSTER.md).
