# 06 — Results and funnels

Numbers come from our logs and from the large-cluster campaign's tables (`pipelines/02_rfd4_steering/campaign_results/`). Counts for stages that
were assessed with thresholds we re-applied ourselves (section 2) are labelled as such.


## 1. Maxi's cluster campaign (35 finished runs, 6 motif/arm combinations; source tables/funnel.csv)
| Stage | Count | Lost to |
|---|---|---|
| RFD4 backbones generated (unsteered / steered / steered+channel) | 10,500 | |
| design gate "pass_geom" (ligand+base SASA<=0.10, pocket ratio, no clash, motif RMSD<=0.3 A) | 3,754 (36 %) | burial / motif drift |
| design gate "pass_core" (+ no polar/Lys/Arg near base, nitro H-bond far from base, TMol dG<0) | 1,654 (16 %) | sequence-independent base-environment rules (unsteered arm: 7 of 600 hgnt) |
| LigandMPNN sequences (native + 4 redesigns + shell-fixed) refolded with Boltz-2, 6NT, 1x | 11,385 | |
| strict Boltz gate (pLDDT>=80, CA<=2 A, ligand ipTM>=0.6, ligand RMSD<=2 A, motif tips<=1 A, refold geometry) | 153 (1.3 %) (+219 relaxed "fallback") | median ligand RMSD vs design 3-4.5 A: Boltz puts the ligand elsewhere |
| full suite (10x TSA, 10x substrate, 10x product; ipSAE2, PLACER, tmol, PyRosetta) | 528 sequences (368 main + 160 extended) | |
| ipSAE2(TSA) > 0.8 | 143 (55 also strict) | |
| tier-1 design gate + strict gate + ipSAE2>0.8 + all 7 composite metrics | 49 | |
| oxyanion donor bound (Gln NE2..N1 <= 3.5 A, best sample and >=50 % of samples) | 41 | donor |
| ranked by composite z -> top 30 | 30 | |

## 2. Our overnight unsteered control (job 60, Maxi's with-Trp spec, evaluated under his definitive gates)
| Stage | Count | Lost to |
|---|---|---|
| RFD4 backbones | 250 | |
| ligand SASA<=0.10 | 249 | |
| + base SASA<=0.10 | 158 | base exposed (-91) |
| + no clash | 133 | clash (-25) |
| + motif RMSD<=0.3 A | 55 | motif drift (-78) |
| top 100 designs -> 500 sequences, Boltz x1: CA<=2 A | 437 | fold (-63) |
| ligand RMSD<=2 A | 68 | ligand pose != design (-369) |
| motif tips<=1 A | 8 | (-60) |
| ligand ipTM / pLDDT | 7 | (-1) |
| ipSAE2(TSA) > 0.8 | 4 (128__rfd4, 230__mpnn4, 226__rfd4, 126__mpnn1) | |
Caveat: only 27 of the 100 refolded designs were among the 55 that pass the geometry gate; 28 passing designs were never refolded.

## 3. Our earlier rounds (2026-11-21, 32-design batches)
spec without burial: 32 -> 12 pass (20 ligand on surface); + BURIED flag, full motif: 32 -> 31; tip-atom + BURIED: 32 -> 16 (6 donor too far, 4 wall too far); His-motif spec: 32 -> 3 (not used).
LigandMPNN x4 per backbone, Boltz refold funnel (jobs 27/28), then 5-seed host-style scoring of a 106-candidate pool (30 GDM-theozyme designs excluded).

## 4. Partial diffusion (job 62)
12 parents x 8 reps at 2 noise levels -> 192 children -> 112 pass the backbone filter -> 64 selected -> 256 sequences -> 16 beat their parent on the composite -> 14 re-scored with 3 seeds. None reached the final 20 (best child ranked 22 in the merged pool).

## 5. Final selection (merged pool → top 20)
| Stage | Count | Lost to |
|---|---|---|
| scored candidates (our pools 144 + Maxi overnight top 30) | 174 rows / 170 unique sequences | |
| dominant loss term (host-style score) | contact distances 143, approach angle/plane 31 | |
| base O..N3 within ~2.2-3.0 A | 69 | |
| donor within ~2.7-3.3 A | 35 | |
| both | 24 | |
| judging stand-in >= 0.01 / >= 0.001 | 11 / 21 | |
| merged top 80 | 80 | |
| eligibility (TM<=0.70, identity<0.70 vs PDB and UniRef90) | 8 removed (7 over TM 0.70, 6 of them by less than 0.01 (0.702-0.708; one at 0.721); 1 with 96 % UniRef90 identity) | |
| judging score x literature factor -> top 20 | 20 | |
Final 20: 10 Maxi overnight, 6 Maxi round 4/5, 4 ours (3 overnight control, 1 earlier shortlist); 10 distinct backbones (structural clustering of the RFD4 designs, US-align TM).
Caveat: the UniRef90 check is missing for the 45 newest sequences (job cancelled at the deadline); TM and PDB identity are done for all 80.

## 6. The submitted 20 (slot order = final rank)

`SubstrateVentures01…20` in `results/final_top20_structures/` and `pipelines/11_selection_and_submission/results/`. "Judging" = the stand-in score (05 B);
"Lit. factor" = L over the common 8 metrics (05 C; "n/a (0.5)" = no substrate/product refolds, neutral). Geometry columns are 5-seed means for our designs and the campaign's 10 TSA samples for his.

| Slot | Steering | Seq. origin | Source design id | Length | base / oxy | base O..N3 (A) | donor (A) | Judging | Lit. factor | TM max |
|---|---|---|---|---|---|---|---|---|---|---|
| 01 | steered | shell1 | `hgnt__steer_tmol_ch18__b02_c0025_no_trp_13_beam3_6__shell1` | 131 | 79 / 49 | 2.66 | 3.1 | 0.0867 | 0.813 | 0.597 |
| 02 | steered | mpnn2 | `hgnt__steer_tmol__b07_c0050_no_trp_3_beam3_0__mpnn2` | 173 | 148 / 31 | 2.52 | 2.95 | 0.1212 | 0.565 | 0.64 |
| 03 | steered | mpnn1 | `hgnt__steer_tmol__b07_c0025_no_trp_3_beam3_5__mpnn1` | 185 | 141 / 69 | 3.04 | 3.07 | 0.0693 | 0.875 | 0.438 |
| 04 | unsteered | rfd4 | `kemp4_with_trp_128__rfd4` | 161 | 31 / 89 | 2.68 | 2.86 | 0.0394 | 0.619 | 0.533 |
| 05 | unsteered | mpnn3 | `kemp4_with_trp_126__mpnn3` | 215 | 35 / 72 | 2.66 | 3.19 | 0.0351 | n/a (0.5) | 0.433 |
| 06 | steered | shell1 | `hgnt__steer_tmol_ch18__b04_c0025_no_trp_23_beam3_6__shell1` | 165 | 68 / 96 | 2.74 | 2.98 | 0.0312 | 0.534 | 0.524 |
| 07 | steered | mpnn2 | `hgnt__steer_tmol__b07_c0050_no_trp_3_beam3_9__mpnn2` | 173 | 148 / 31 | 2.96 | 3.98 | 0.0220 | 0.446 | 0.639 |
| 08 | steered | mpnn3 | `hgnt__steer_base_tmol__no_trp_0_beam3_6__mpnn3` | 127 | 77 / 50 | 3.41 | 3.01 | 0.0167 | n/a (0.5) | 0.632 |
| 09 | steered | mpnn4 | `hgnt__steer_tmol__b07_c0050_no_trp_3_beam3_0__mpnn4` | 173 | 148 / 31 | 2.63 | 3.42 | 0.0297 | 0.281 | 0.64 |
| 10 | steered | mpnn3 | `hgnt__steer_base_tmol__no_trp_0_beam3_5__mpnn3` | 127 | 77 / 50 | 2.68 | 3.32 | 0.0742 | 0.073 | 0.629 |
| 11 | steered | rfd4 | `hgnt__steer_tmol_ch18__b04_c0025_no_trp_23_beam3_6__rfd4` | 165 | 68 / 96 | 2.85 | 3.0 | 0.0064 | 0.661 | 0.565 |
| 12 | unsteered | mpnn4 | `full-af14b3c83c24_s4` | 197 | 91 / 49 | 2.6 | 5.75 | 0.0053 | 0.767 | 0.627 |
| 13 | unsteered | mpnn1 | `hgnt__unsteered__b01_c0150_no_trp_19__mpnn1` | 171 | 133 / 114 | 2.63 | 2.94 | 0.0152 | 0.204 | 0.515 |
| 14 | steered | shell1 | `hgnt__steer_tmol_ch18__b04_c0025_no_trp_23_beam3_11__shell1` | 165 | 68 / 96 | 2.74 | 3.01 | 0.0088 | 0.205 | 0.469 |
| 15 | steered | mpnn4 | `hgnt__steer_base_tmol__no_trp_0_beam3_5__mpnn4` | 127 | 77 / 50 | 3.92 | 4.54 | 0.0033 | n/a (0.5) | 0.63 |
| 16 | steered | mpnn3 | `hgnt__steer_tmol__b07_c0025_no_trp_3_beam3_4__mpnn3` | 185 | 141 / 69 | 2.81 | 2.81 | 0.0027 | 0.605 | 0.444 |
| 17 | unsteered | mpnn1 | `kemp4_with_trp_107__mpnn1` | 168 | 62 / 41 | 2.69 | 4.65 | 0.0018 | n/a (0.5) | 0.6 |
| 18 | steered | mpnn1 | `hgnt__steer_base_tmol__no_trp_0_beam3_5__mpnn1` | 127 | 77 / 50 | 2.77 | 3.32 | 0.0017 | n/a (0.5) | 0.623 |
| 19 | steered | mpnn1 | `hgnt__steer_base_tmol__no_trp_0_beam3_6__mpnn1` | 127 | 77 / 50 | 2.66 | 3.4 | 0.0017 | n/a (0.5) | 0.623 |
| 20 | steered | rfd4 | `hgnt__steer_base_tmol__no_trp_0_beam3_5__rfd4` | 127 | 77 / 50 | 3.18 | 2.97 | 0.0009 | 0.497 | 0.633 |

Generation method per design (exact description incl. seeds and motif) is generated by `pipelines/11_selection_and_submission/build_final_table_lit.py` / `rerank_merged80.py` (column `method` of `rerank_top20.csv`).

**Diversity:** the 20 sequences sit on **10 distinct backbones** (US-align on the RFD4 design backbones, single linkage; same result at TM 0.6/0.8/0.9):
backbone A = slots 08, 10, 15, 18, 19, 20 (round 4/5 beam members); B = 06, 11, 14; C = 02, 07, 09; D = 03, 16; singletons 01, 04, 05, 12, 13, 17.
Origins: 10 overnight-campaign designs (9 steered, 1 unsteered), 6 round 4/5 (steered), 4 ours (3 from the unsteered overnight control, 1 earlier shortlist).

## 7. Which ligand takes the pocket? (Boltz-2, TSA + substrate + product cofolded together)

Each submitted design with 6NT (chain B), H5J substrate (C) and 6VP product (D) in one complex; no MSA, seed 0, 5 diffusion samples
(`runs/…top20_all-three-ligands-together_boltz2`, summary reproduced here). "In pocket" = ligand centroid within 4.0 Å of where the TSA sat in that
design's single-ligand refold (after Cα superposition) and an atom within 4.5 Å of the base carboxylate O.

| Slot | TSA | Substrate | Product | Preferred | mean ipTM (TSA / SUB / PROD) |
|---|---|---|---|---|---|
| 01_steered_shell1 | 0.2 | 0.0 | 1.0 | PROD | 0.418 / 0.451 / 0.565 |
| 02_steered_mpnn2 | 0.0 | 0.0 | 1.0 | PROD | 0.356 / 0.363 / 0.706 |
| 03_steered_mpnn1 | 1.0 | 0.0 | 0.0 | TSA | 0.523 / 0.373 / 0.436 |
| 04_unsteered_rfd4 | 0.6 | 0.2 | 0.2 | TSA | 0.429 / 0.382 / 0.367 |
| 05_unsteered_mpnn3 | 0.0 | 0.0 | 1.0 | PROD | 0.465 / 0.424 / 0.506 |
| 06_steered_shell1 | 1.0 | 0.0 | 0.0 | TSA | 0.535 / 0.394 / 0.363 |
| 07_steered_mpnn2 | 0.0 | 0.0 | 1.0 | PROD | 0.304 / 0.38 / 0.697 |
| 08_steered_mpnn3 | 0.0 | 0.0 | 1.0 | PROD | 0.437 / 0.444 / 0.606 |
| 09_steered_mpnn4 | 1.0 | 0.0 | 0.0 | TSA | 0.5 / 0.419 / 0.419 |
| 10_steered_mpnn3 | 0.0 | 0.0 | 1.0 | PROD | 0.454 / 0.481 / 0.645 |
| 11_steered_rfd4 | 0.2 | 0.4 | 0.4 | SUB | 0.445 / 0.447 / 0.456 |
| 12_unsteered_mpnn4 | 1.0 | 0.0 | 0.0 | TSA | 0.631 / 0.42 / 0.412 |
| 13_unsteered_mpnn1 | 1.0 | 0.0 | 0.0 | TSA | 0.59 / 0.52 / 0.406 |
| 14_steered_shell1 | 1.0 | 0.0 | 0.0 | TSA | 0.548 / 0.388 / 0.367 |
| 15_steered_mpnn4 | 0.2 | 0.0 | 0.8 | PROD | 0.436 / 0.395 / 0.513 |
| 16_steered_mpnn3 | 1.0 | 0.0 | 0.0 | TSA | 0.575 / 0.459 / 0.364 |
| 17_unsteered_mpnn1 | 0.8 | 0.0 | 0.0 | TSA | 0.525 / 0.431 / 0.517 |
| 18_steered_mpnn1 | 1.0 | 0.0 | 0.0 | TSA | 0.537 / 0.411 / 0.441 |
| 19_steered_mpnn1 | 1.0 | 0.0 | 0.0 | TSA | 0.532 / 0.411 / 0.454 |
| 20_steered_rfd4 | 0.0 | 0.0 | 1.0 | PROD | 0.395 / 0.388 / 0.673 |

TSA is preferred by 11 of 20 designs (9 with 5/5 samples), the product by 8, the substrate by 1 (never cleanly). One seed, 5 samples; Boltz tends to fill the pocket with one
ligand and park the others, so this shows preference, not binding energy; the product's higher ipTM partly reflects its smaller size.

## 8. Other results worth citing

- **Literature correlation / BioEmu** (`pipelines/09_bioemu/`): no metric survives selection-aware correction (best BioEmu global Cα RMSF ρ = +0.70, raw p = 0.011, corrected p = 0.18). BioEmu ensembles of three campaign designs: protein rigid, ligand bimodal (15–23 % bound), base-end swings 2.5–2.8 Å.
- **Steering vs unsteered** (campaign): `pass_core` 16–30 % (steered) vs 1–3 % (unsteered); the unsteered arm lost most designs to base-environment rules.
- **Our control vs the campaign:** under the campaign's definitive gates our 250 unsteered backbones gave 4 designs through ipSAE2 > 0.8 (`runs/…overnight-gate-survivors`); one of them (`kemp4_with_trp_128__rfd4`) is slot 04.
- **Partial diffusion:** 16 of 256 child sequences beat their parent on the composite; none reached the final 20.
- **Eligibility:** of the merged top-80, 8 were removed (7 over TM 0.70 — 6 of them by less than 0.01 (0.702–0.708), one at 0.721 — and one with 96 % UniRef90 identity). Our own pool had no PDB-seqres identity above the limit. UniRef90 was not completed for the 45 newest sequences before the deadline.
