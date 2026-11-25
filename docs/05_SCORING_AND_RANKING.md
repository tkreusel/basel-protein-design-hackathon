# 05 — Scoring and the final ranking

The hosts' formulas are not public. Everything here is a stand-in built from their written criteria and calibrated on crystal
structures. Component values are always kept so the ranking can be re-weighted.

## A. Hosts' criteria (from `reference/SUBMISSION_hosts_spec.md`)

1. catalytic contact distances (base → N3, donor → N1) — "carries the most weight";
2. approach geometry — angle and in-plane term, the worse of the two decides;
3. iPDE (interface predicted distance error; lower is better);
4. active-site packing (base enclosed in a pocket);
5. pose consistency across repeated predictions.
"Geometry outweighs everything else combined; the terms multiply." Each contact has an optimum (too short is penalised harder than too long);
designs are folded with MSAs.

## B. Judging stand-in (`pipelines/07_host_style_scoring/`)

total = geometry² × s_iPDE × s_burial × s_consistency (mean over Boltz seeds for the geometry terms)

| Term | Definition |
|---|---|
| contact | Gaussian on base O…N3 distance (reference 2.55 Å) × Gaussian on donor…N1 distance (reference 3.05 Å); σ floor 0.25 Å; the short side uses σ × 0.6 (too short penalised harder) |
| approach | for base and donor: min(angle score, out-of-plane score); product of the two |
| geometry | contact × approach |
| s_iPDE | clip((4 − iPDE) / 3.5, 0, 1) (iPDE from Boltz `complex_ipde`; our 425 refolds: 0.42–5.19 Å, median 1.20) |
| s_burial | min(1, CB-count within 12 Å of the base / 45.9) (45.9 = crystal reference) |
| s_consistency | seed-to-seed reproducibility: median pairwise ligand RMSD (weight 2) plus spreads of base distance, donor distance, burial, iPDE |

`configs/host_reference.yaml` (core) holds the crystal-derived ideal geometry (Asp127/Gln50 contacts from 8 crystal copies; crystals score 0.89–0.95 in the
earlier mean-of-Gaussians version). The first stand-in `host_score.py` used an additive mean (`composite_abs`); `rank_v2.py` is the multiplicative version used for the final ranking because the hosts say terms multiply.

Judgement calls (not validated against the hosts): the 2.55 / 3.05 Å optima, squaring geometry, the 0.6 short-side factor, the iPDE and burial scalings.
Bug history: `geom_scores` once overwrote raw distances with 0–1 scores, inflating the consistency term; fixed 2026-11-22 04:30 UTC (see `pipelines/07_host_style_scoring/README.md`).

## C. Literature factor L

Source: univariate regressions of the evaluation-suite metrics against kcat/KM over 19 literature enzymes (ranked by leave-one-out Q²).
For a design and metric: z = (value − mean over the 19 enzymes) / sd, clipped to ±3, oriented by the sign of the regression correlation, weighted by LOO Q².
L = Φ( Σ wᵢzᵢ / Σ wᵢ ) over the metrics available for that design (0–1; 0.5 if none).

Metrics in the final ranking (the 8 available for every design from Boltz arms + ipSAE2):

| Metric | Sign | Weight (Q²) |
|---|---|---|
| ipSAE protein→ligand, substrate arm | + | 0.467 |
| ligand ipTM, substrate − product | + | 0.411 |
| ligand ipTM, TSA − substrate | − | 0.404 |
| catalytic-residue→ligand PAE, substrate arm | − | 0.398 |
| PAE catalytic→ligand, TSA − substrate | + | 0.342 |
| ipSAE, TSA − substrate | − | 0.325 |
| ipSAE, substrate − product | + | 0.271 |
| ligand ipTM, substrate arm | + | 0.259 |

An earlier 12-metric version also used tmol E_protein (sub − prod, −, 0.281), lk_ball (TSA − prod, +, 0.258), the PyRosetta field `rlxcst.field.CH.env_nobase` (−, 0.325)
and the Boltz catalytic→pocket PAE (sub − prod, −, 0.549); those need tmol/PyRosetta/features runs that only existed for ~42 designs, so the final ranking used the common 8.
Not used: PLACER donor median (its positive correlation, +0.775, contradicts the criteria and is probably an outlier effect) and the ipSAE iPAE12 rows (not in the campaign tables).

Honest status of L: a selection-aware permutation test (labels shuffled 20,000×, maximum |ρ| over the metric family) found no BioEmu or
literature metric that survived correction (best raw p = 0.011 → corrected 0.18). L is a tie-breaker built on 19 enzymes, not a validated predictor.

## D. The final ranking — weights

final = judging × L — no additive weights; the two factors multiply (as the hosts say terms do).

| Factor | Where its weights live | Effective influence |
|---|---|---|
| Judging (geometry², iPDE, burial, consistency) | inside the formula (geometry squared ≈ twice the combined weight of the rest) | spans ~10⁵× across the 80-design pool (sd of log = 4.1) → ≈ 85–90 % of the variation in log(final) |
| Literature factor L | LOO-Q² weights above | spans ~10× (sd of log = 1.3–1.8) → ≈ 10–15 % |

Sensitivity: the final ranking correlates with the judging-only ranking at Spearman 0.988; top 10 share 9 designs and top 20 share 17. Adding the
twelfth metric (catalytic→pocket PAE) changed the earlier ranking very little (Spearman 0.998). Beyond about rank 12 the scores are weakly separated
(judging scores ≈ 0.001–0.003 vs 0.03–0.12 for the top five), so order within ranks 13–20 is close to arbitrary.

Eligibility is applied before ranking (TM ≤ 0.70, identity < 0.70 against PDB and UniRef90), identical sequences are collapsed to the best-ranked copy, and
designs with no substrate/product arm data get L = 0.5.

## E. Consistency, seeds and samples

- Host-style scores use 5 independent seeds (separate runs, not diffusion samples of one trunk pass — those would share a trunk and overstate consistency).
- Campaign designs: 10 Boltz TSA samples.
- Substrate/product arms for the literature factor: single seed, single sample (noisier than the geometry terms).
- No design was checked with several diffusion samples within a run.

## F. Reproduce
`pipelines/07_host_style_scoring/rank_v2.py` → `ranking_v2.csv`; `pipelines/11_selection_and_submission/{merged80_pool,rerank_merged80}.py` → final table.
Script header docstrings state their exact inputs.
