# How to read the campaign evidence, and what to test next

## What "passing" means (tiers, strongest first)
1. `final_strict_all_gates_ipsae0.8/` : passed the arm's DESIGN gate (tier 1, no fallback) + STRICT Boltz-1x gate (pLDDT>=80, CA<=2 A, ligand ipTM>=0.6, ligand RMSD<=2 A, catalytic tips<=1 A, refold geometry) + ipSAE2(TSA)>0.8; ranked by the 7-metric composite. This is the folder to inspect first.
2. `final_top25_ipsae0.8_composite/` : ipSAE2>0.8 ranked by composite, but includes sequences that were never gated (extended set) or came through the fallback; their violations are listed per sequence (motif RMSD>1 A or refold geometry).
3. `shortlist_all_criteria/`, `passing_strict_boltz_gate/`, `tsa_preferred_over_substrate/`, `top50_by_panel_rank/` : earlier cuts (see README).
`recommended_for_testing.csv` (in the results folder) is a diversified pick from tier 1 (<=2 sequences per backbone, catalytic tips <=1.0 A, ligand RMSD <=1.5 A).

## What the evidence is worth
Reasonably informative (structural self-consistency): design geometry vs theozyme (snapshots: design tips 0.1-0.4 A, Boltz tips 0.4-0.9 A), Boltz refold agreement (CA RMSD, pocket RMSD, ligand RMSD, catalytic-tip RMSD), ligand ipTM/ipSAE2/ipDE consistency across 10 samples, PyRosetta dG_sep (independent of the TMol reward).
Weak / exploratory: tmol TSA-minus-substrate/product energy deltas (single relax, noise ~5 REU; TMol is also the steering reward), PLACER fractions (not calibrated: the best known enzyme HG3.17 scores ~0.05 holo; within-backbone variance ~50 % of total), Boltz TSA-vs-substrate ipTM differences (saturated ipTM 0.85-0.99; product 6VP is almost indistinguishable from the TSA 6NT).
Not measured at all: catalysis, water network / protonation, loop dynamics, expression, stability, selectivity, actual TS-vs-ground-state stabilisation (6NT is a neutral TSA analogue, 6VP is the neutral phenol).
The composite is a ranking heuristic, not a kcat predictor. Differences between neighbouring ranks are within noise; treat the top ~10-15 as an unordered shortlist.

## Inspecting a candidate
1. Overlay design_RFD4.cif, boltz_TSA_6NT_best.cif and the theozyme (inputs/HG3_17_kemp_motif.pdb, align on 6NT): Asp base carboxylate O to N3 ~2.5-2.7 A and Gln50 NE2 toward the nitro/N1 side as in HG3.17; Boltz must reproduce it (tips <=1 A).
2. Pocket: ligand buried (SASA <=0.10), no unsatisfied buried polars on the nitro, channel bottleneck 1.7-1.9 A in the channel arm; look for ordered-water access.
3. Base environment: no polar side chains within 4.5 A of the base carboxylate (the designed intent); shell-fixed variants may keep native polar residues at 4.5-6 A.
4. Consistency across ligands: substrate (H5J) and product (6VP) models should sit in the same pocket; a large ligand RMSD between TSA and substrate refolds = no pocket specificity.
5. Red flags: only the extended/fallback tier; ligand ipTM high but ligand RMSD>1.5 A; catalytic tips>1 A; PyRosetta dG_sep > -14; Cys, Met-rich or very hydrophobic surface; all candidates from one backbone (redundant); hgwt with Trp44 clashing in the Boltz model.

## Orthogonal checks before wet-lab work
- A second co-folding oracle with ligand (for example Chai-1, or ESMFold2 with a CCD ligand): ligand pose and catalytic tips on the top 12.
- AF2 initial guess for fold sanity of the protein only (no ligand support).
- Sequence checks: free Cys, N-glyc motifs, aggregation/hydrophobic patches, ESM pseudo-log-likelihood, nearest natural homologue (foldseek/BLAST) to confirm novelty vs HG3.
- Re-run tmol/PyRosetta with 3-5 relaxes for the finalists to put error bars on dG_sep and the composite terms.
- Controls to order with the designs: HG3.17 (positive), catalytic Asp->Asn/Ala (negative), a scrambled sequence of one design, and one unsteered design (to test whether steering matters).

## Recommended set (recommended_for_testing.csv)
12 sequences from 7 backbones (<=2 per backbone) in the order of the composite: 10 hgnt + 2 hgwt; 7 channel-arm + 5 steered-only; mix of full and shell-fixed LigandMPNN variants so the shell-fixed idea gets tested. Backups = next-ranked sequences of the same backbones (different sequences = lower redundancy risk).
