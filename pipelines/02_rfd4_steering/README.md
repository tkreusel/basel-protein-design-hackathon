# 02 RFD4 campaign: steered and unsteered backbone generation, sequence design, refolding, evaluation

This directory holds the pipeline that produced the large design campaign (`20261122_kemp_c1`): 10,500 RFD4-Proteina backbones,
13,000 LigandMPNN sequences refolded with Boltz-2, and a full physics and ensemble evaluation of 528 sequences.
It replaces the earlier hand-over snapshot of the round 4/5 code, which is kept in [`../../archive/early_rounds_r4_r5`](../../archive/early_rounds_r4_r5).

## What it does

```
RFD4-Proteina, 3 arms x 2 motifs        gen.py            conditioning only, or beam search on a reward set
  -> design gate                        filter_pre.py     burial, pocket, clashes, motif RMSD, base environment (fallback to 30 designs/run)
  -> LigandMPNN                         mpnn.py           native + 4 full redesigns (+ 1 shell-fixed for steered arms)
  -> Boltz-2, TSA (6NT), 1 sample       boltz_prep.py
  -> refold gate                        boltz_gate.py     pLDDT, CA RMSD, ligand ipTM, ligand RMSD, motif tips, refold geometry (fallback to 10/run)
  -> survivors: Boltz-2 9x TSA, 10x substrate (H5J), 10x product (6VP)
  -> forest + report                    select_models.py  samples ranked by pocket RMSD to the design
  -> evaluation suite                   ../06_evaluation_suite   ipSAE2, PLACER (holo and apo), tmol, PyRosetta
  -> collection and selection           collect.py, report.py, flags.py, final_scoring.py, final_selection.py, export_results.py
```

Motifs: `hgnt` (HG3.17 Asp127 and Gln50) and `hgwt` (plus Trp44), both positioned around the transition-state analogue 6NT.
Arms:

| arm | generation |
|---|---|
| `unsteered` | conditioning only (burial, hydrophobic contact, nitro acceptor, fixed unindexed tips), 300 replicates |
| `steer_tmol` | beam search (4 branches, beam 3, checkpoints at denoising steps 30, 42, 52, 60) on twelve geometry rewards plus the TMol interface energy; 100 replicates x 3 beam members |
| `steer_tmol_ch18` | as `steer_tmol` plus a bottleneck-radius band of 1.7 to 1.9 A and a weaker ligand-SASA term |
| `<arm>_sheet` | the same arms with beta-sheet conditioning on the expandable segment, as a small diversity add-on |

All thresholds, reward weights and seeds are in [`config/`](config/): `arms.yaml`, `gates.yaml`, `motifs.yaml`, `sequence_design.yaml`.
The three beam members of a steered replicate are one backbone (CA RMSD about 0.03 A); the unit of replication is the lineage, not the design.

## Layout

| path | content |
|---|---|
| `config/` | arms and rewards, gate thresholds, motif definition, LigandMPNN settings, `site.env.example` |
| `src/kemp/` | one module per stage (see above); `config.py` and `structure.py` are shared helpers |
| `jobs/` | shell drivers: `stage_A.sh`, `stage_B.sh`, `stage_P.sh`, `stage_X.sh`, `worker.sh`, `submit.sh`, `status.sh`, `archive_runs.sh`, `repair.sh`, `install_rfd4.sh` |
| `inputs/` | RFD4 input specs (four per motif/sheet combination), the HG3.17 motif, the 6NT TMol parameters, the other theozymes explored |
| `rfd4_patch/` | Kemp reward metrics, the metric registry, and the `pyproject` patch applied to the RFD4 clone |
| `tests/` | 29 pytest checks and one shell check; fixtures are excerpts of the campaign |
| `figures/` | the campaign figures (one script each) and the shared style |
| `results_summary/` | per-run funnel, per-arm summary, the final top 30, notes on how to read the evidence |

## Running it

```bash
cd pipelines/02_rfd4_steering
cp config/site.env.example config/site.env      # set scratch, results, interpreters, tool locations
bash jobs/install_rfd4.sh                       # clone RFD4 at the pinned commit, patch, install metrics, write the experiment configs
export KEMP_PIPELINE=$PWD; source jobs/lib.sh

sbatch jobs/smoke.sbatch                        # full chain at tiny scale; check it before a large launch
bash jobs/submit.sh --motifs "hgnt hgwt" --batches "1 2 3" --sheet 1 --nodes 2          # dry run: prints the task list
bash jobs/submit.sh --motifs "hgnt hgwt" --batches "1 2 3" --sheet 1 --nodes 2 --submit # whole-node workers, one per GPU
bash jobs/status.sh
python -m kemp.export_results                   # (PYTHONPATH=src) tables, sequences and structure folders in $KEMP_RESULTS
```

A task is `(motif, arm, batch)`, 300 designs. Every stage skips outputs that exist, so a killed job resumes where it stopped.
Workers choose work in the order PLACER, stage B, stage A, and (optionally) the extended evaluation; locks are directories created with `mkdir`,
and a lock whose Slurm job is gone is taken over. We used whole-node jobs with eight workers because single-GPU jobs waited many hours in the queue;
two or three worker processes per GPU gave roughly twice the throughput because the GPUs were CPU-bound.
PLACER ran only on H200 GPUs (torch 2.3 with cu121 gave "no kernel image" on the others); set `PLACER_GPU_PATTERN` for your hardware.

Measured on H200 GPUs: steered generation 70 to 100 s per replicate (300 designs in 2 to 3 h), unsteered about 7 s per design, stage B about 1 to 2 min per
design for ipSAE2, tmol and PyRosetta, PLACER about 30 s per 20-model ensemble.

## Tests

```bash
PYTHONPATH=src python -m pytest tests -q        # needs pandas, PyYAML, biotite, pytest
bash tests/test_jobs.sh
```

The checks reproduce what the campaign produced from excerpts of its own tables: the design gate and the refold gate (every flag and tier of three
runs, including the fallback tiers), the seeds of all generated chunks, the experiment YAML files byte for byte, the motif, omit and shell residue sets of
the sequence design, the stage-5 report, the extended set, the shortlist flags, the composite ranking (139 of 139 ranks), the panel rank (528 of 528) and the
geometry of the final candidates. They do not run RFD4, Boltz-2, LigandMPNN or the evaluation suite; those need the GPU environments.

## Results

Counts per arm over the non-sheet runs (35 runs; designs are the generated backbones):

| motif | arm | designs | lineages | pass geometry (%) | pass core (%) | sequences refolded | strict Boltz pass | designs with a strict pass |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| hgnt | unsteered | 600 | 600 | 31 | 1.2 | 925 | 10 | 6 |
| hgnt | steer_tmol | 3000 | 1000 | 50 | 22.9 | 4116 | 102 | 67 |
| hgnt | steer_tmol_ch18 | 3000 | 1000 | 37 | 13.4 | 2442 | 34 | 25 |
| hgwt | unsteered | 600 | 600 | 22 | 2.8 | 650 | 0 | 0 |
| hgwt | steer_tmol | 3000 | 1000 | 40 | 25.3 | 2478 | 4 | 4 |
| hgwt | steer_tmol_ch18 | 3000 | 1000 | 26 | 14.8 | 2670 | 26 | 16 |

Steering raises the share of backbones that satisfy the base-environment rules by roughly an order of magnitude, but those rules are also the reward terms,
so this is the expected outcome of the search and not evidence of catalytic competence. The loss that matters comes later: Boltz-2 without an MSA often places
the ligand away from the designed pose (median ligand RMSD to the design 3 to 4.5 A), and 153 of 11,385 sequences pass the strict refold gate (219 more enter
through the fallback). Of the 528 sequences that received the full evaluation, 143 have ipSAE2(TSA) above 0.8, 55 of them also strict; the tier-1 set that
satisfies every gate has 49 members, 41 of which bind the oxyanion donor to N1 in the Boltz TSA models. The top 30 of these are in
[`results_summary/final_top30`](results_summary/final_top30). How far the evidence reaches is discussed in
[`results_summary/INTERPRETATION.md`](results_summary/INTERPRETATION.md).

## Limitations

- No design was tested experimentally. Boltz-2 ipTM is saturated (0.85 to 0.99), PLACER fractions are not calibrated (the crystal of HG3.17 scores about 0.05 holo), and tmol is both a reward and an evaluation term.
- Unsteered designs are gated on geometry only, steered designs on the full core filter; both flag sets are stored for every design so the arms can be compared under one filter afterwards.
- The fallback tiers exist so that every run yields sequences; they are marked (`gate_tier`, `boltz_tier`) and the final selection uses tier 1 only.
- The composite score is a ranking heuristic over seven metrics. Differences between neighbouring ranks are within noise.
- PyRosetta Sc is unreliable and is not used; tmol relaxes are stochastic (about 5 REU).
