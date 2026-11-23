# Team Substrate|Ventures: de novo Kemp eliminase design

Benedict Wolf, Friederike Reinhardt, Maximilian Vassen, Paul Brunner, Tobias Kreusel

Code, configuration and documentation of our entry to the Protein Design Hackathon at D-BSSE (ETH Zurich, November 2026). The task was to design up to 20 protein sequences
(50 to 260 residues) that Boltz-2 co-folds with the transition-state analogue 6NT such that a catalytic carboxylate sits next to the abstracted atom N3 and an oxyanion donor next to N1.
Designs had to be novel (below 70 % sequence identity to published enzymes, TM-score at most 0.7 to published structures). The hosts ranked by catalytic contact distances and approach
geometry, interface confidence (iPDE), active-site burial and pose consistency over repeated predictions. Our 20 sequences, on 10 distinct backbones, won the event. No design has been tested experimentally.

## Approach

Backbones come from RFD4-Proteina conditioned on the HG3.17 active site (Asp127, Gln50, optionally Trp44) around 6NT, either with conditioning alone or with an inference-time beam search on
geometry rewards and the TMol interface energy. Sequences are designed with LigandMPNN with the motif fixed and polar and charged residues excluded around the catalytic base, and refolded with
Boltz-2 (no MSA). Survivors are evaluated with substrate and product refolds, ipSAE2, PLACER ensembles, tmol and PyRosetta; candidates are filtered for novelty against the PDB and UniRef90 and ranked
by a stand-in for the hosts' score multiplied by a literature-derived factor.

```
theozyme (Asp127, Gln50, [Trp44], 6NT)
  -> RFD4 backbones: conditioning only | steered beam search | partial diffusion
  -> design gate (burial, pocket, clashes, motif RMSD, base environment)
  -> LigandMPNN (motif fixed; K R H D E N Q S T Y C excluded near the base)
  -> Boltz-2 refold with 6NT -> refold gate -> 10x TSA, 10x substrate, 10x product
  -> ipSAE2, PLACER, tmol, PyRosetta; host-style score; BioEmu
  -> eligibility (identity, TM-score) -> literature-calibrated ranking -> top 20
```

The campaign that produced most of the final designs generated 10,500 backbones in six motif/arm combinations. Steering raised the fraction of backbones passing the base-environment gate from 1 to 3 %
to 13 to 25 %, but the pass criteria are also the rewards, and the dominant loss is downstream: the median ligand RMSD between Boltz-2 refold and design is 3 to 4.5 A. Of 11,385 refolded sequences, 153 pass
the strict refold gate. Details and numbers are in [pipelines/02_rfd4_steering](pipelines/02_rfd4_steering) and [docs/06_RESULTS_AND_FUNNELS.md](docs/06_RESULTS_AND_FUNNELS.md).

## Repository

| Path | Content |
|---|---|
| [`pipelines/02_rfd4_steering/`](pipelines/02_rfd4_steering) | the campaign pipeline: configuration, stage modules, job scripts, tests, figures, result tables. Start here to re-run the generation |
| [`pipelines/06_evaluation_suite/`](pipelines/06_evaluation_suite) | ipSAE2, PLACER, tmol and PyRosetta evaluation |
| [`pipelines/01_backbones_sequence_refold/`](pipelines/01_backbones_sequence_refold), [`04_overnight_unsteered_control/`](pipelines/04_overnight_unsteered_control), [`05_partial_diffusion/`](pipelines/05_partial_diffusion) | early rounds, the one-job unsteered control, partial diffusion |
| [`pipelines/07_host_style_scoring/`](pipelines/07_host_style_scoring) | stand-in for the hosts' score, multi-seed consistency, calibration on crystals |
| [`pipelines/08_eligibility/`](pipelines/08_eligibility) | identity and TM-score novelty gate with local databases |
| [`pipelines/09_bioemu/`](pipelines/09_bioemu) | BioEmu-2 ensembles and the literature correlation |
| [`pipelines/11_selection_and_submission/`](pipelines/11_selection_and_submission) | merged pool, literature-calibrated ranking, FASTA validation and upload |
| [`core/`](core) | geometry library, motif and burial checks, calibrated thresholds, tests, the 6NT ligand |
| [`results/`](results) | structures of the 20 submitted designs |
| [`docs/`](docs) | method, scoring, results and lessons; `docs/reference/` holds working notes from the event |
| [`archive/`](archive) | superseded round 4/5 steering code |

## Using it

The campaign pipeline and the evaluation suite are configured through `pipelines/02_rfd4_steering/config/site.env` and contain no machine-specific paths.
The scripts under `core/` and the other pipelines were copied as run and still refer to the hackathon cluster (`/mnt/work/...`); [docs/02](docs/02_ENVIRONMENTS_AND_CLUSTER.md) lists what to change.
Model weights, raw outputs and search databases are not included; [docs/08](docs/08_PROVENANCE_AND_DATA_LOCATIONS.md) says where each lives and how to regenerate it.

```bash
cd pipelines/02_rfd4_steering
cp config/site.env.example config/site.env && $EDITOR config/site.env
bash jobs/install_rfd4.sh
PYTHONPATH=src python -m pytest tests -q        # checks the gates, seeds and rankings against the campaign's own tables
```

| To | Read |
|---|---|
| understand what we did and why | [docs/00](docs/00_PROJECT_OVERVIEW.md), [docs/01](docs/01_PIPELINE_MAP.md) |
| re-run a stage | [docs/03](docs/03_REPRODUCING.md) |
| compare conditioning, steering and partial diffusion | [docs/04](docs/04_GENERATION_METHODS.md) |
| see how candidates were scored and ranked | [docs/05](docs/05_SCORING_AND_RANKING.md) |
| see the funnels and the submitted list | [docs/06](docs/06_RESULTS_AND_FUNNELS.md) |
| avoid our mistakes | [docs/07](docs/07_LESSONS_AND_PITFALLS.md) |

## Limitations

The ranking rests on an in-silico stand-in for the hosts' score, whose formulas are not public, and on a literature factor from 19 enzymes that did not survive a selection-aware permutation test.
Boltz-2 ipTM is saturated for these designs, PLACER fractions are uncalibrated, and the steering rewards overlap with the pass criteria.

## Third-party tools

RFD4-Proteina, LigandMPNN, Boltz-2, ipSAE, PLACER, tmol, BioEmu, MMseqs2, Foldseek, US-align and PyRosetta are used as installed and are not redistributed. PyRosetta is free for non-commercial
academic use under the RosettaCommons licence. Check each tool's licence before reusing its outputs commercially.
