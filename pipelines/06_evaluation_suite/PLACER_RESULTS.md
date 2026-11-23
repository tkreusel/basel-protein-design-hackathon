# PLACER results (20 samples per ensemble; 17 replicates; 12 entries)

apo pre-formed = fraction of the apo ensemble (side chains predicted without ligand) in which the catalytic carboxylate sits within the
crystal-derived window of the reference TSA pose: Asp O...N3 2.18-2.77 A, O...H 1.25-1.76 A, O...H-N angle 147-180 deg, C-O...H angle 116-140 deg
(all four must hold). 20 models -> binomial uncertainty about +-0.1.

| entry | kcat/KM | apo pre-formed | apo Asp-N3 d (A) | apo base chi SD (deg) | holo in-geometry | holo prmsd (A) | holo lig RMSD to input (A) | holo pose spread (A) | holo donor in range |
|---|---|---|---|---|---|---|---|---|---|
| 7K4Q HG3 | 146 | 0.02 | 2.06 +-0.48 | 21 | 0.00 | 3.30 | 2.22 | 2.60 | - |
| 5RGA HG3 | 146 | 0.02 | 2.04 +-0.14 | 2.8 | 0.00 | 3.25 | 2.42 | 2.72 | - |
| 7K4X HG3.7 | 27,000 | 0.05 | 2.10 +-0.08 | 5.9 | 0.02 | 3.37 | 2.43 | 2.70 | 0.25 |
| 7K4Z HG3.17 | 126,000 | 0.20 | 2.07 +-0.13 | 4.2 | 0.05 | 3.33 | 2.32 | 2.59 | 0.35 |
| 8RD5 HG3.R5 | n/a | 0.55 | 2.29 +-0.15 | 12 | 0.00 | 3.45 | 2.27 | 2.40 | 0.12 |
| 4BS0 HG3.17 E47N/N300D | n/a | 0.05 | 2.15 +-0.12 | 8.5 | 0.10 | 3.28 | 2.32 | 2.66 | 0.42 |
| 29SB KABLE2.5 (NMR) | n/a | 0.00 | 2.93 +-0.10 | 3.0 | 0.00 | 3.71 | 3.54 | 3.12 | - |
| 7K4Q Boltz | 146 | 0.15 | 2.12 +-0.08 | 4.3 | 0.00 | 3.36 | 2.37 | 2.52 | - |
| 5RGA Boltz | 146 | 0.25 | 2.17 +-0.10 | 3.0 | 0.00 | 3.36 | 2.72 | 2.83 | - |
| 4BS0 Boltz | n/a | 0.45 | 2.46 +-0.11 | 5.1 | 0.00 | 3.28 | 2.08 | 2.41 | 0.35 |
| 8RD5 Boltz | n/a | 0.50 | 2.44 +-0.09 | 3.6 | 0.00 | 3.08 | 1.89 | 2.26 | 0.75 |
| 29SB Boltz | n/a | 0.10 | 2.67 +-0.39 | 28 | 0.15 | 3.45 | 2.50 | 2.78 | 0.55 |

Correction (made after the first analysis): `holo donor in range` was 1.00 everywhere because donor distances were taken from the input
structure instead of each model; it is now computed per model (values above).
