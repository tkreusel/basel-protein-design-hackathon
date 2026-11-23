# 01 Pipeline map

```
theozyme + ligand (core/, pipelines/02/inputs)
  -> S1  RFD4 backbones           conditioning only | reward-steered beam search | partial diffusion
  -> S1' design gate              burial, pocket, clashes, motif RMSD, base environment
  -> S2  LigandMPNN               motif fixed, polar and charged residues excluded near the base
  -> S3  Boltz-2 refold with 6NT  no MSA; gate on fold, ligand pose, motif tips, geometry
  -> S4  evaluation               substrate and product refolds, ipSAE2, PLACER, tmol, PyRosetta, host-style score, BioEmu
  -> S5  selection                merged pool, eligibility, literature-calibrated ranking, FASTA, validator
```

| Stage | Purpose | Code | Environment | Notes |
|---|---|---|---|---|
| Foundation | atom map, calibrated thresholds, geometry and burial scoring, tests | `core/src/kemp/`, `core/src/check_motif.py`, `check_burial.py`, `derive_atom_map.py`, `milestone1.py`, `core/configs/`, `core/tests/` | CPU | validated on the HG series ([reference/MILESTONE1.md](reference/MILESTONE1.md)) |
| S1a conditioning only | backbones around the fixed motif and ligand | `pipelines/02_rfd4_steering/src/kemp/gen.py` (campaign); `pipelines/01.../jobs/stage1_rfd4.sbatch`, `pipelines/04...` (early rounds and the one-job control) | RFD4 venv, one GPU | recipe in [reference/TOOL_NOTES.md](reference/TOOL_NOTES.md) |
| S1b steering | reward-guided beam search, 4 branches and beam 3 | `pipelines/02_rfd4_steering` (`config/arms.yaml`, `rfd4_patch/`) | RFD4 venv and tmol | [04](04_GENERATION_METHODS.md) |
| S1c partial diffusion | re-noise good designs and regenerate the tail | `pipelines/05_partial_diffusion` | one GPU | no child reached the final 20 |
| S1' design gate | cheap CPU rejection before any refold | `pipelines/02.../src/kemp/filter_pre.py`; `pipelines/04.../select_designs.py` | CPU | `pass_geom`, `pass_seq`, `pass_core`, fallback tiers |
| S2 LigandMPNN | sequences per backbone | `pipelines/02.../src/kemp/mpnn.py`; `core/src/run_ligandmpnn.py` | LigandMPNN env | T 0.1; native + 4 redesigns (+ 1 shell-fixed for steered arms) |
| S3 Boltz-2 refold | does the sequence refold to the design with the ligand in place | `boltz_prep.py`, `boltz_gate.py` | Boltz-2 env, no MSA, 3 recycles | the dominant loss |
| S4a host-style score | stand-in for the hosts' ranking | `pipelines/07_host_style_scoring` | CPU + multi-seed Boltz | geometry squared x iPDE x burial x consistency |
| S4b evaluation suite | ipSAE2, PLACER, tmol, PyRosetta | `pipelines/06_evaluation_suite`, driven by `jobs/stage_B.sh` | four environments | PLACER needs a GPU |
| S4c eligibility | identity below 0.70, TM at most 0.70 | `pipelines/08_eligibility` | MMseqs2, Foldseek, US-align | local only; the published HG3.17 control must fail |
| S4d BioEmu, literature | ensemble behaviour, literature correlation | `pipelines/09_bioemu` | BioEmu-2 | no signal after correction |
| S5 selection | pool, filter, rank, FASTA | `pipelines/11_selection_and_submission`; campaign: `kemp.final_scoring`, `kemp.final_selection` | CPU | `validate-submission.py` |

## File contracts between stages

| From to | Artifact |
|---|---|
| S1 to gate | RFD4 `*_generated.cif` (protein and ligand; the conditioning columns including `C_CRD` are kept because partial diffusion needs them) |
| gate to S2 | `gate.csv` and `sel.csv` (campaign) or `selection.csv` (early rounds), plus the design CIFs |
| S2 to S3 | `sequences.csv` (`seq_id, sequence, fixed, ...`) and Boltz YAML inputs (`msa: empty`, ligand `ccd: 6NT`) |
| S3 to S4 | `boltz_out/**/predictions/<id>/<id>_model_N.cif`, `confidence_*.json`, `pae_*.npz`; curated into a forest by `select_models.py` |
| S4 to S5 | per-candidate tables (`ranking_v2.csv`, `eligibility_*.csv`, `candidates.csv`, evaluation-suite tables) |
| S5 to upload | `designs.fasta` with `base=` and `oxy=` headers, checked by `validate-submission.py` |

## Resources

On the hackathon cluster one Slurm partition with dynamic AWS nodes was available and the team budget was 8 vCPU of GPU workers at a time, so GPU stages were chained and
CPU-only stages ran on the head node ([reference/CLUSTER.md](reference/CLUSTER.md)). The campaign ran on a larger cluster with whole 8-GPU nodes
([02](02_ENVIRONMENTS_AND_CLUSTER.md)).
