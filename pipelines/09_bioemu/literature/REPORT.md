# BioEmu + Boltz/ipSAE metrics vs literature kcat/KM (2026-11-22)

Data: BioEmu-2 ensembles (1000 samples each, protein + 6NT) of 12 HG3-family enzymes (HG3 ... HG3.17, HG3.R1-R5,
HG3.17wR5, HG3.5w17) and Boltz-2 / ipSAE2 tables (10 models each) from the zips in `bioemu/new_zip/`; labels
`/mnt/work/manually_created_dataset.csv` (HG3.3/7/14 back-calculated from a figure). 1A53-core was not analysed with BioEmu
(different scaffold, Asp129 is Val there). Crystal numbering + 2 in the 312-residue sequences (Asp127 -> 129, Gln50 -> 52),
confirmed by Maxi's session. Scripts: `ensemble_metrics.py`, `correlate.py`; table `correlations.csv`; per-enzyme `metrics/*.json`.

## Result: nothing survives a correction for picking the best metric

74 tests (52 metrics, two enzyme sets; BioEmu metrics only on the 12 HG3 variants). Spearman rho against log10 kcat/KM, n = 12 (HG3 family):

| metric | rho | raw p | selection-aware p* |
| --- | --- | --- | --- |
| BioEmu median CA RMSF (global flexibility) | +0.70 | 0.011 | 0.18 |
| BioEmu base burial (CA within 12 A of Asp129 CG) | -0.56 | 0.061 | 0.55 |
| BioEmu base O...N3 median, ligand-in-pocket frames | -0.50 | 0.094 | 0.69 |
| Boltz ligand ipTM / ipSAE2 (mean of 10 models) | +0.50 | 0.094 | 0.65 |
| Boltz ipDE | -0.45 | 0.14 | 0.79 |
| BioEmu pocket RMSF | +0.45 | 0.14 | 0.82 |
| BioEmu near-attack fraction (base O-N3 <= 3.5 A, buried ligand) | +0.12 | 0.72 | 1.0 |

*permutation test: labels shuffled 20,000 times, p = share of shuffles whose BEST |rho| over the whole metric family
(BioEmu 24 metrics / Boltz 25) is at least as large. BH q over all tests >= 0.78.

## What the ensembles show

- Ligand occupancy of the designed pocket varies a lot between variants (27-82 % of samples) but not with activity; the most
  active enzyme HG3.17 has the LOWEST occupancy (27 %, the ligand leaves in most samples) while HG3.R5 (best) has 78 %.
- Near-attack frames (base O within 3.5 A of N3, ligand buried): 15-32 % in 11 of 12 enzymes (7 % in HG3.17), no trend (rho +0.12).
- Base-geometry score and pocket occupancy: rho -0.25 and -0.21, not significant.
- Positive RMSF correlation (more flexible ensemble = more active) is the opposite of a preorganisation story and is
  confounded by the evolution path (later variants carry more mutations); it is not a usable design criterion.

## Our designs against the literature range (BioEmu, same metrics)

| | pocket occupancy | near-attack fraction | base O-N3, pocket median (A) |
| --- | --- | --- | --- |
| literature (12 HG3 variants) | 0.27-0.82 | 0.07-0.32 | 3.8-4.9 |
| our three Maxi designs | 0.15-0.31 | 0.01-0.04 | 5.4-6.4 |

BioEmu rates the designs below every literature enzyme on near-attack frames. Caveats: designs are 127 aa (burial counts are not
comparable), BioEmu's ligand handling is untested for de novo folds, and the literature set shows no calibration to activity.

## Caveats

- n = 12, one scaffold, evolutionarily ordered variants: any monotone metric of "number of mutations" will correlate.
- Labels from two sources disagree between families (Maxi: HG3 146 vs 4100 here); ranks inside the HG3 family agree.
- Maxi's session found the same pattern with PLACER, tmol and Rosetta (n = 19): only arm differences (e.g. tmol TSA - product
  pocket energy, Boltz substrate - product PAE) carry signal; single-arm absolute values do not.
