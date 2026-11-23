# Kemp eliminase structures in the PDB

Found by full-text search of RCSB for "Kemp eliminase" (97 entries), then
metadata-checked entry by entry. **6NT presence was verified by grepping the
downloaded coordinate file, not inferred from the title** — several entries say
"complexed with transition state analog" but use a different analogue.

Reproduce: `src/find_kemp_structures.sh` (search) + `src/derive_atom_map.py`.

## Downloaded and scored — `benchmark/pdb/` (6NT-bound)

All share the HG3 xylanase scaffold numbering (base Asp127, donor site 50,
wall Trp44), so `geometry.py` reads them without modification.

| PDB | Variant | Res (Å) | base→N3 | donor→N1 | composite | kcat/KM |
|---|---|---|---|---|---|---|
| 7K4Q | HG3 | 1.28 | 2.25 | 6.05 | 0.150 | 146 |
| 7K4X | HG3.7 | 1.60 | 2.49 | 2.96 | 1.000 | 27,000 |
| 7K4Z | HG3.17 | 1.08 | 2.51 | 3.10 | 0.937 | 126,000 |
| 8RD5 | HG3.R5 | 1.50 | 2.57 | 3.03 | 1.000 | — |
| **4BS0** | HG3.17 Q50 E47N/N300D **(parent)** | **1.09** | 2.53 | 3.10 | 0.939 | — |
| **6SRW** | HG3.17 **Q50F** E47N/N300D | 1.45 | 2.37 | — none | 0.150 | — |
| **6TU6** | HG3.17 **Q50M** E47N/N300D | **1.05** | 2.55 | — none | 0.360 | — |
| 6SS1 | HG3.17 **Q50A** E47N/N300D | — | 2.61 | — none | 0.400 | — |
| 6SS3 | HG3.17 **Q50K** E47N/N300D | — | 2.52 | 7.37 | 0.400 | — |
| **7K4U** | **HG3 K50Q** | 1.30 | 2.54 | 3.20 | 1.000 | — |
| 6NW4 | evolved computational design | 3.00 | — | — | — | — |
| 9TZD | KempTIM1 (de novo TIM barrel) | 1.20 | n/a | n/a | n/a | — |

Bold = added after the first Milestone 1 pass.

**6NW4 fails loud** (`no atoms found for ligand 6NT atom N3`) — at 3.0 Å the
triazole is incompletely modelled. It is excluded rather than scored with a
substituted default. This is `geometry.py`'s fail-loud policy working as
intended.

**9TZD is a different scaffold** (de novo TIM barrel), so residues 44/50/127 are
meaningless there. Its catalytic set must be located from coordinates before it
can be scored — see "worth doing next" below.

### Why 4BS0 + 6SRW/6TU6/6SS1/6SS3 matters

This is a **clean single-variable series**: five structures, identical
HG3.17 E47N/N300D background, differing only at position 50
(Q → F / M / A / K). It removes the caveat on the earlier Milestone 1 result,
where 6SS1/6SS3 carried extra mutations relative to 7K4Z and so were not
controlled comparisons. They *are* controlled against 4BS0.

Result: parent composite **0.939** vs every position-50 mutant **≤ 0.400**, and
the donor metric is absent entirely for F/M/A. Criterion 2 is now properly
established rather than approximated.

**7K4Q vs 7K4U** is a second clean pair: a single K50Q substitution on the HG3
background moves the donor metric from **6.05 → 3.20 Å**. That is precisely the
mutation that created HG3.7's oxyanion donor, and the metric sees it.

## Apo / conformational-state structures — `benchmark/apo/`

These bear directly on why Milestone 1 criterion 1 fails. The HG3.7→HG3.17 gain
is a preorganisation effect invisible in a single holo structure (plan §6.8),
and these are the structures that show it:

| PDB | What it is | Res (Å) |
|---|---|---|
| 7K4T | HG3.17, apo | **0.999** |
| 7K4W | **HG3.17 in the inactive state** | 1.90 |
| 9U01 | KempTIM1, apo | 1.25 |

Also available, not yet downloaded: 7K4P (HG3 apo, 1.08), 7K4R (HG3 K50Q apo,
1.57), 7K4S (HG3.7 apo, 2.00), 7K4V (HG3.17, 1.30), 7K4Y (HG3.17 at 343 K, 1.80).

**7K4W is the interesting one.** An explicitly deposited "inactive state" of
HG3.17 gives us the unproductive conformation directly from crystallography.
Measuring our geometric metrics on 7K4T/7K4V (productive) vs 7K4W (inactive)
calibrates what "unproductive" looks like numerically — which is the threshold
PLACER's `placer.geom_frac` and BioEmu's apo productive fraction would need
anyway. **We may be able to get the preorganisation discriminator from crystals
instead of from BioEmu-2 (plan stage 7, v2 scope).**

## The other literature designs from plan §10

| Plan §10 entry | Structure(s) | Notes |
|---|---|---|
| HG3 (146) | **7K4Q** (+6NT), 7K4P (apo) | have it |
| HG3.7 (27,000) | **7K4X** (+6NT), 7K4S (apo) | have it |
| HG3.17 (126,000) | **7K4Z**, **4BS0**, 7K4T/V/W/Y | have it, richly |
| HG3.3b (2,200) | **none found** | no deposited structure located |
| HG3.14 (52,000) | **none found** | no deposited structure located |
| HG4 (103,000) | **not identified** | 6NW4 "Evolution of a computationally designed Kemp eliminase" (3.0 Å) is the candidate but the title does not name the variant |
| Des27 (130) / Des27.7 (12,700) | **9HVB, 9HVH** — "High-efficiency Kemp eliminases by complete computational design" (2.0/2.1 Å) | almost certainly this series; **neither contains 6NT** |
| 2008 series, 6–160 | KE07: 4Z08, 5D2V–5D38, 6C7H/M/T/V, 6C8B, 6CT3, 6DKV, 6DNJ (directed-evolution series). KE70: 3NPU/V/W/X, 3NQ2/8, 3NQV, 3NR0. KE59: 3UXA, 3UXD, 3UYC, 3UY7/8, 3UZ5, 3UZJ | large evolution series available |
| — | HG-2: **3NYD** (1.23 Å) | TSA-bound but **not** 6NT |
| — | 1A53-core: 8FMC (apo), 8FOQ (+TSA, not 6NT) | another scaffold family |
| — | HG649: 8USE and the 8US* series | |
| — | KABLE2.5: 29SB, 9SMT (NMR solution) | de novo designed |
| — | **KempTIM1/4: 9TZD (+6NT, 1.20), 9U01 (apo), 9QKX** | **de novo TIM barrels — closest published analogue to our own task** |
| — | HG3.17 + 5-cyanobenzotriazole: 8KH3 | alternative TSA |
| — | P450 BM3-TMK + 5-nitrobenzisoxazole: 7Y9J/7Y9K | different scaffold, real substrate bound |

Two full-text false positives: **7EPN/7EPO** are ketosteroid isomerase (a
natural enzyme), not Kemp eliminases.

## Worth doing next with these

1. **HG3.3b and HG3.14 have no structures**, so Milestone 1 criterion 1 stays at
   n = 3 from crystals alone. Raising n needs either published kcat/KM for
   variants we *do* have structures for (4BS0, 7K4U, 6SRW, 6TU6, 6SS1, 6SS3,
   8RD5 — all currently `None`), or Boltz-2 predicted structures for HG3.3b and
   HG3.14. **Sourcing the activity values is far cheaper than predicting
   structures** and would take criterion 1 to n ≈ 8 immediately.
2. **Score 7K4T/7K4V vs 7K4W** to calibrate productive-vs-inactive geometry from
   crystals. Potentially substitutes for stage 7.
3. **Locate KempTIM1's catalytic residues** in 9TZD from coordinates. It is a de
   novo TIM barrel with 6NT bound at 1.20 Å — an independent check that our
   metrics transfer off the HG3 scaffold, which is exactly what they must do to
   score our own designs.
4. The KE07/KE70/KE59 evolution series offers many more geometry points if their
   analogue ligands can be mapped onto 6NT's roles via `atom_map.yaml`.

> Activity values are deliberately left as `None` in `src/milestone1.py` where we
> have not verified them against the literature. Filling one in from memory
> would silently corrupt the only quantitative gate we have.
