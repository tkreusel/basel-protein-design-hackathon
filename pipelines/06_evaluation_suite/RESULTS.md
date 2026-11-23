# Results: 17 replicates, 12 entries (crystals/NMR + Boltz refolds)

kcat/KM (M-1 s-1) from the project plan (Broom et al. 2020): HG3 146, HG3.7 27,000, HG3.17 126,000. **No literature value entered for HG3.R5, HG3.17 E47N/N300D or KABLE2.5.**

## Table 1: crystal / NMR structures

| entry | variant | kcat/KM | donors<3.5A | dBase-N3 | PLACER apo pre-formed | PLACER holo in-geometry | tmol ddE(TSA-S) | Rosetta ddG(TSA-S) | tmol dE(P-S) | H-bond E (tmol) | packstat |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 7K4Q | HG3 | 146 | 0.0 | 2.61 | 0.03 | 0.00 | -3.92 | -2.03 | 1.61 | -1.86 | 0.63 |
| 5RGA | HG3 (277 K) | 146 | 0.0 | 2.59 | 0.03 | 0.00 | -2.68 | -2.53 | 2.56 | -1.83 | 0.75 |
| 7K4X | HG3.7 | 27,000 | 2.0 | 2.58 | 0.05 | 0.03 | -5.09 | -2.05 | 2.91 | -4.07 | 0.79 |
| 7K4Z | HG3.17 | 126,000 | 1.0 | 2.64 | 0.20 | 0.05 | -3.35 | -1.35 | 1.42 | -3.73 | 0.85 |
| 8RD5 | HG3.R5 | - | 1.5 | 2.66 | 0.55 | 0.00 | -3.18 | -1.36 | 1.67 | -3.81 | 0.78 |
| 4BS0 | HG3.17 E47N/N300D | - | 1.0 | 2.59 | 0.05 | 0.10 | -5.25 | -1.93 | 2.07 | -4.15 | 0.78 |
| 29SB | KABLE2.5 (de novo, NMR) | - | 0.0 | 2.66 | 0.00 | 0.00 | -4.67 | -0.36 | 0.65 | -1.92 | 0.65 |

## Table 2: Boltz-2 refolds (single chain, 1 sample each)

| entry | kcat/KM | ipSAE2 p->l | ligand ipTM | pLDDT catalytic | PAE cat->lig (A) | dBase-N3 | donors | PLACER apo | tmol ddE | Rosetta ddG | tmol dE(P-S) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 7K4Q_boltz | 146 | 0.90 | 0.98 | 0.99 | 0.74 | 2.57 | 0.0 | 0.15 | -4.08 | -1.30 | -0.48 |
| 5RGA_boltz | 146 | 0.87 | 0.98 | 0.99 | 0.86 | 2.57 | 0.0 | 0.25 | -5.01 | -1.74 | -1.78 |
| 8RD5_boltz | - | 0.89 | 0.98 | 0.98 | 0.79 | 2.60 | 1.0 | 0.50 | -3.96 | -0.57 | -0.21 |
| 4BS0_boltz | - | 0.85 | 0.97 | 0.99 | 1.00 | 2.60 | 1.0 | 0.45 | -2.34 | 0.04 | -0.53 |
| 29SB_boltz | - | 0.90 | 0.96 | 0.98 | 0.79 | 2.63 | 1.0 | 0.10 | -3.93 | -2.32 | -0.89 |

Full per-metric tables: `work/report/validation.md`, `validation_ladder.csv`, `validation_classes.csv`, `summary_replicates.csv`.

