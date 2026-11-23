# Kemp eliminase de novo design — pipelines, configs and docs

Everything we built for the Kemp-eliminase design hackathon (Nov 2026): backbone generation with
**RFD4-Proteina** (conditioning-only and reward-steered), sequence design with **LigandMPNN**, refolding with **Boltz-2**,
a physics/ensemble evaluation suite (**ipSAE2, PLACER, tmol, PyRosetta**), a stand-in for the hosts' scoring, an eligibility
gate (sequence identity / TM-score against published structures), BioEmu ensembles, a literature calibration, and the
final selection that produced the 20 submitted designs (which won).

> **Status:** this is a snapshot of working scripts, copied from the project folders as they were used.
> Paths inside scripts are the cluster paths we ran on (`/mnt/work/projects/...`, `/mnt/work/repos/...`, `/mnt/ref/...`);
> see [docs/02_ENVIRONMENTS_AND_CLUSTER.md](docs/02_ENVIRONMENTS_AND_CLUSTER.md) for the path map and what to change to re-run elsewhere.
> Large data (model weights, raw RFD4 output, Boltz outputs, PDB mirrors, search databases) is **not** in this repo;
> [docs/08_PROVENANCE_AND_DATA_LOCATIONS.md](docs/08_PROVENANCE_AND_DATA_LOCATIONS.md) says where each lives and how to regenerate it.

## The task in one paragraph

Design ≤ 20 protein sequences (50–260 aa) that, co-folded with the transition-state analogue **6NT** by Boltz-2,
place a catalytic carboxylate (Asp/Glu) next to the abstracted atom (N3) and an oxyanion donor (Gln/Asn/Ser/…) next to the
developing oxyanion (N1). Designs must be novel: < 70 % sequence identity to published enzymes and TM-score ≤ 0.7 to
published structures. The hosts rank by catalytic contact distances and approach geometry (terms multiply and dominate),
interface confidence (iPDE), active-site burial and pose consistency across repeated predictions
([docs/reference/SUBMISSION_hosts_spec.md](docs/reference/SUBMISSION_hosts_spec.md)).

## The pipeline at a glance

```
 theozyme (HG3.17: Asp127 + Gln50 [+ Trp44] + 6NT)                                  core/, pipelines/01
        │
        ▼
 RFD4-Proteina backbones ──┬─ conditioning only (unsteered)                          pipelines/01, 04
                           ├─ reward-steered beam search (steer_tmol, +channel band)  pipelines/02
                           └─ partial diffusion around good designs                  pipelines/05
        │ design gate (burial, motif RMSD, base environment, clash)
        ▼
 LigandMPNN (motif fixed; K R H D E N Q S T Y C forbidden near the base)             pipelines/01, 04
        │
        ▼
 Boltz-2 refold with 6NT (no MSA), 1x gate → 5 seeds → substrate/product arms        pipelines/01, 07
        │
        ├─► host-style stand-in score  (geometry² · iPDE · burial · consistency)     pipelines/07
        ├─► evaluation suite (ipSAE2, PLACER, tmol, PyRosetta, features)             pipelines/06
        ├─► eligibility gate (PDB/UniRef90 identity, Foldseek + US-align TM)         pipelines/08
        └─► BioEmu ensembles + literature calibration                                pipelines/09
        │
        ▼
 merged top-80 → eligibility filter → judging × literature factor → top 20           pipelines/11 → results/
```

## Repository map

| Path | What is in it |
|---|---|
| [`docs/`](docs/) | numbered documentation (start with `00`), plus `docs/reference/` = original project notes kept verbatim |
| [`core/`](core/) | shared library and foundation: geometry (`src/kemp`), motif/burial/atom-map scripts, calibrated thresholds, tests, `env/paths.sh`, the 6NT ligand |
| [`pipelines/01_backbones_sequence_refold/`](pipelines/01_backbones_sequence_refold/) | stage 1–5 Slurm jobs and RFD4 input specs (theozyme conditioning) |
| [`pipelines/02_rfd4_steering/`](pipelines/02_rfd4_steering/) | reward-steered RFD4 (beam search), custom metrics, round drivers, and the campaign result tables from the large cluster run |
| [`pipelines/04_overnight_unsteered_control/`](pipelines/04_overnight_unsteered_control/) | one-job unsteered control: generate → select → LigandMPNN → Boltz → score |
| [`pipelines/05_partial_diffusion/`](pipelines/05_partial_diffusion/) | partial diffusion around best designs with a CPU rejection filter |
| [`pipelines/06_evaluation_suite/`](pipelines/06_evaluation_suite/) | the ipSAE2 / PLACER / tmol / PyRosetta suite and our additions and patches |
| [`pipelines/07_host_style_scoring/`](pipelines/07_host_style_scoring/) | stand-in for the hosts' score, multi-seed Boltz consistency, calibration against crystals |
| [`pipelines/08_eligibility/`](pipelines/08_eligibility/) | identity and TM-score novelty gate with local tools only |
| [`pipelines/09_bioemu/`](pipelines/09_bioemu/) | BioEmu-2 ensemble analysis and the literature-correlation study |
| [`pipelines/11_selection_and_submission/`](pipelines/11_selection_and_submission/) | merged pool, literature-calibrated re-ranking, FASTA/validator/submit scripts, tables |
| [`results/`](results/) | the 20 submitted designs' structures (Boltz refold and RFD4 backbone) |

(Numbers `03` and `10` are intentionally unused: sequence design and refolding live inside `01`/`04`; the literature calibration
lives with BioEmu in `09` and `11`.)

## Result

The submission was 20 sequences on **10 distinct backbones** (US-align clustering of the RFD4 designs, identical at TM 0.6/0.8/0.9):
10 designs from the big-cluster overnight campaign (9 steered, 1 unsteered), 6 from its earlier steered round (round 4/5, all steered), 4 from our own unsteered runs. The list, with
generation method per design, is in [docs/06_RESULTS_AND_FUNNELS.md](docs/06_RESULTS_AND_FUNNELS.md).

## Where to start

| If you want to … | read |
|---|---|
| understand what we did and why | [docs/00_PROJECT_OVERVIEW.md](docs/00_PROJECT_OVERVIEW.md), [docs/01_PIPELINE_MAP.md](docs/01_PIPELINE_MAP.md) |
| re-run something | [docs/03_REPRODUCING.md](docs/03_REPRODUCING.md) |
| understand conditioning vs steering vs partial diffusion | [docs/04_GENERATION_METHODS.md](docs/04_GENERATION_METHODS.md) |
| understand the scores and the final ranking weights | [docs/05_SCORING_AND_RANKING.md](docs/05_SCORING_AND_RANKING.md) |
| see the funnel (how many designs survived each stage) | [docs/06_RESULTS_AND_FUNNELS.md](docs/06_RESULTS_AND_FUNNELS.md) |
| avoid our mistakes | [docs/07_LESSONS_AND_PITFALLS.md](docs/07_LESSONS_AND_PITFALLS.md) |

## Third-party tools and licences

RFD4-Proteina (NVIDIA), LigandMPNN, Boltz-2, ipSAE, PLACER, tmol, BioEmu, MMseqs2, Foldseek, US-align and PyRosetta are
used as installed on the cluster and are **not** redistributed here. PyRosetta is free for non-commercial academic use under
the RosettaCommons licence — check your use qualifies. Check each tool's licence before reusing outputs commercially.
