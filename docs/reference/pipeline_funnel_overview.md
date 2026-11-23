# Pipeline funnel overview (for the pitch) - numbers from our logs and Maxi's campaign tables, 2026-11-22

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

## 5. Final selection
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

## Distinct backbones in the final 20 (US-align on the RFD4 design backbones, single linkage; identical result at TM 0.6, 0.8 and 0.9)
10 backbones: {08,10,15,18,19,20} (6 designs, Maxi round 4/5 beam members), {06,11,14} (3), {02,07,09} (3), {03,16} (2), and singletons 01, 04, 05, 12, 13, 17.
