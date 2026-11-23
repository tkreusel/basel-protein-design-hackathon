# 02 — RFD4 reward-steered generation (large-cluster campaign)

Contents are the **handoff from the steering-campaign owner** (code as of 2026-11-21 ~21:10, round 4/5) plus the **campaign's final summary tables**.

| Path | What |
|---|---|
| `README_AGENT_HANDOFF.md` | full description: goal, motifs, steering semantics, arms, rewards, cluster facts, results and caveats — **read first** |
| `code/make_r4.py`, `make_r5.py` | write the per-arm RFD4 experiment YAMLs |
| `code/run_r4.py`, `run_r5.py` (+ `*_job.sbatch`) | generate → flatten designs → per-design metrics CSV, one process per job |
| `code/eval_r4.py`, `eval_r5.py`, `filter_designs.py` | aggregate to arm summaries; the design gate (core criteria) |
| `code/downstream/` | LigandMPNN prep (`prep_mpnn.py`), Boltz chunks (`prep_boltz.py`, `run_boltz_chunk.sbatch`), `eval_boltz.py`, `select_passers.py` |
| `code/rfd3/` | the earlier RFD3 Feynman–Kac attempt (abandoned: pocket proxy over-tightened, ESS collapsed to 1–2/12) |
| `rfd4_patch/custom_metrics/` | the reward metrics (`kemp_base`, `kemp_pocket`, `kemp_chem`, `kemp_channel`, `kemp_tmol`) — drop into `rfd4/src/rfproteina/metrics/custom_metrics/` |
| `rfd4_patch/configs/inference_experiment/` | beam-search experiment configs per arm; `pyproject.diff` removes `bionemo-ir` (needs GLIBC 2.34) |
| `inputs/` | the three motif specs (`kemp4_hg317-tips_{with,no}-trp.yaml`, `kemp4_gdm1483-glu-ser.yaml`), theozyme structures, `6NT.tmol` |
| `early_rounds_metrics/` | metrics of the earlier rounds |
| `campaign_results/` | **the overnight campaign** (`steer_tmol`, `steer_tmol_ch18`, `unsteered` × `hgnt`/`hgwt`, 10,500 backbones): `funnel.csv` (per run), `arm_summary.csv`, `CAMPAIGN_README.md`, `PIPELINE_OVERVIEW.md`, `INTERPRETATION_AND_RECOMMENDATION.md` |

**Gap to know about:** the code that ran the *final overnight campaign* (300 designs/run on the big cluster, arms `steer_tmol` and `steer_tmol_ch18`, beam search with the
channel band 1.7–1.9 Å, shell-fixed LigandMPNN variant, strict Boltz gate, ipSAE2/PLACER/tmol/PyRosetta evaluation, `make_results_folder.py`) lives in the owner's own
repository on his cluster and was not available to us beyond the summary tables and the round-4/5 handoff. The *method* is documented in
`campaign_results/PIPELINE_OVERVIEW.md` and `docs/04_GENERATION_METHODS.md`; the steering mechanism, rewards and gates above are, to our knowledge, the same code lineage (not verified against his repository).
