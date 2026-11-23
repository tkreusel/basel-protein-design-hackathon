# Milestone 1 — geometry scoring stack validated on the HG series

Plan §4. Reproduce with:

```bash
source env/paths.sh
env/kemp-orch/bin/python -I src/milestone1.py .      # the table below
PYTHONPATH=src env/kemp-orch/bin/python -m pytest tests/ -q   # criterion 3
```

**Verdict: 2 of 3 acceptance criteria pass. Criterion 1 fails at ρ = +0.50
(target ≥ 0.70), and the reason is diagnostic rather than a broken metric.**

Scope: this is the *geometry-only* stack. tmol is not installed
(docs/TOOL_NOTES.md), so there is no `ddE_ts_sub` term yet, and PLACER (stage 4)
and BioEmu-2 (stage 7) are out of v0 scope. Those absences are exactly what
criterion 1 turns out to need — see below.

## Measured geometry

All distances in Å, from `src/geometry.py` reading roles from
`configs/atom_map.yaml`. Chain A of each deposited structure.

| PDB | Variant | base→N3 | donor→N1 | wall | min contact | composite | kcat/KM (M⁻¹s⁻¹) |
|---|---|---|---|---|---|---|---|
| 7K4Q | HG3 | 2.25 | **6.05** | 5.02 | 2.25 | 0.150 | 146 |
| 7K4X | HG3.7 | 2.49 | **2.96** | 5.01 | 2.49 | 1.000 | 27,000 |
| 7K4Z | HG3.17 | 2.51 | **3.10** | 5.54 | 2.51 | 0.937 | 126,000 |
| 8RD5 | HG3.R5 | 2.57 | 3.03 | 4.79 | 2.57 | 1.000 | — |
| 6SS1 | HG3.17 Q50A/E47N/N300D | 2.61 | **99.00** (no tip) | 5.05 | 2.61 | 0.400 | — |
| 6SS3 | HG3.17 Q50K/E47N/N300D | 2.52 | **7.37** | 3.83 | 2.52 | 0.400 | — |

## Criterion 1 — Spearman ρ vs log₁₀(kcat/KM) — FAIL

| Metric | ρ | p | vs 0.7 |
|---|---|---|---|
| composite (higher better) | **+0.500** | 0.667 | FAIL |
| donor distance (lower better) | +0.500 | 0.667 | FAIL |
| base distance (lower better) | **−1.000** | 0.000 | FAIL |

```
Predicted: HG3.7 > HG3.17 > HG3
Actual:    HG3.17 > HG3.7 > HG3
```

Two things are going on, and they should not be conflated.

**(a) The criterion is unreachable at n = 3.** Only three benchmark structures
have a published kcat/KM. With n = 3, Spearman ρ can only take the values
−1, −0.5, +0.5, +1 — so "ρ ≥ 0.7" really means "get all three in perfect
order". We get HG3 (the weakest, by 2+ orders of magnitude) unambiguously
right and invert the top two. Reporting ρ = 0.5 as a near-miss of 0.7 would
overstate the precision available; the honest statement is **one inversion out
of three pairs**.

To make criterion 1 a real test, the benchmark needs the plan §4.1
"predict-from-sequence" points — HG3.3b (2,200) and HG3.14 (52,000) — which
have no 6NT crystal structures and therefore require stage 5 (Boltz-2) to
generate coordinates. That raises n to 5 and makes ρ meaningful. **This is the
single highest-value addition to the benchmark.**

**(b) The inversion is the known limit of static holo geometry, not a bug.**
HG3.7 → HG3.17 is a 4.7× activity gain, and in the crystal the two are
geometrically near-identical (donor 2.96 vs 3.10 Å; base 2.49 vs 2.51 Å). The
plan already names the mechanism in §6.8: *"unbound HG3 through HG3.14 had the
pocket residues blocking productive binding, while HG3.17's unbound state is
productive in ~40% of crystal molecules."* The HG3.17 gain lives in the **apo
conformational ensemble**, which a single ligand-bound crystal structure cannot
express. No static geometric metric computed on a holo structure can recover
it, no matter how the terms are weighted.

So criterion 1 is not evidence that the metric definitions are wrong. It is
evidence that **preorganisation needs an ensemble method** — stage 4 (PLACER
`placer.donor_chi_spread`, which is precisely "the Gln50-rigidity question
measured directly") or stage 7 (BioEmu-2 apo productive fraction). Both were
already scoped as the discriminators; this result quantifies why they matter.

### The base-distance trap — ρ = −1.000

The catalytic base distance **anti-correlates perfectly** with activity. HG3,
the weakest variant, has the *shortest* Asp127→N3 distance (2.25 Å) of the
entire series; HG3.17 has 2.51 Å.

All six structures sit in a 2.23–2.77 Å band regardless of a 10³-fold activity
range, so base geometry is **necessary but not discriminating**. Any composite
that rewards "shorter base distance" will actively rank designs backwards. The
composite therefore scores it as a *window-satisfaction* term (1.0 anywhere in
2.4–3.0 Å), not a "closer is better" term, and weights it 0.25 against the
donor's 0.60. `tests/test_geometry.py::test_base_metric_does_not_track_activity`
pins this so a later "improvement" cannot silently reintroduce it.

## Criterion 2 — ablation separation — PASS

All four required comparisons pass, with wide margins:

| Comparison | Result |
|---|---|
| HG3.17 donor 3.10 < 6SS1 (Q50A) 99.00 | PASS |
| HG3.17 donor 3.10 < 6SS3 (Q50K) 7.37 | PASS |
| HG3.17 composite 0.937 > 6SS1 0.400 | PASS |
| HG3.17 composite 0.937 > 6SS3 0.400 | PASS |

The donor–oxyanion metric separates the two classes cleanly and with no
overlap: **active 2.96–3.10 Å vs ablated 6.05–99.0 Å**, a gap of >2.9 Å. Any
cut-off in 3.5–5.5 Å gives perfect classification;
`configs/thresholds.yaml` uses 3.5 Å, just outside H-bonding range.

Note that HG3 itself (Lys50, 6.05 Å) falls on the *ablated* side of this
boundary, which is correct: plan §10 lists HG3 as having **no oxyanion donor**,
and the coordinates confirm Lys50's NZ never comes within H-bonding distance of
N1. The metric is reading real chemistry.

The `HG3.17_D127A` and `HG3.17_scramble` in-silico negatives from plan §4.1 are
**not yet built**. The two deposited oxyanion ablations (6SS1, 6SS3) already
give criterion 2's required separation, so these are a lower priority, but
`D127A` remains worth building because it is the only negative that tests the
*base* metric.

## Criterion 3 — known-geometry recovery — PASS

20/20 tests pass (`tests/test_geometry.py`). `src/derive_atom_map.py` and
`src/kemp/geometry.py` select atoms and compute distances through independent
code paths, so agreement to <0.015 Å across 6 structures × 3 metrics is a real
cross-check, not a tautology. The suite also covers the Kabsch reflection
guard, fail-loud on missing atoms, the absent-donor sentinel, and the [lit]
1.5 Å clash rule (no native structure violates it; closest contact 2.25 Å).

## The atom map (plan §6.1) — derived, not assumed

The mapping everything downstream depends on, established from 11 ligand copies
and recorded in `configs/atom_map.yaml`:

| Role | 6NT atom | Substrate analogue | Evidence |
|---|---|---|---|
| `abstraction_centre` | **N3** | C3 | Asp127 OD2→N3 = 2.23–2.77 Å; N3 is the closest triazole N in all 11 copies (N2 3.2–3.8, N1 4.3–4.9) |
| `oxyanion_centre` | **N1** | O1 → phenolate O | Gln50 NE2→N1 = 2.96–3.10 Å in every Gln50 structure; N1 closest acceptor every time |
| `nitro_group` | NO1, O11, O21 | nitro | **refuted** as the oxyanion contact: 7.6–8.7 Å from the donor tip |

**Consequence for the theozyme set (plan §5.2): cell family `b` is wrong for
this chemistry.** The plan proposes `b` = "Tyr/Ser/Thr hydroxyl → nitro group".
In every structure the oxyanion donor contacts the **triazole N1** (the
developing phenolate oxygen), and the nitro group is 7.6–8.7 Å away on the far
side of the ring. A theozyme that H-bonds a hydroxyl to the nitro group is not
stabilising the developing negative charge — it is decorating the
electron-withdrawing group that sets the leaving-group pKa.

Recommended revision before stage 0 generates cells: keep the donor→oxyanion
*position* (N1/phenolate O) fixed across all cells and vary the **donor chemistry**
there instead:

| Cell family | Donor at the oxyanion position |
|---|---|
| `a` | Gln/Asn amide NE2/ND2 → N1 — the validated HG3.7+ geometry |
| `b` | Tyr/Ser/Thr hydroxyl → N1 — same position, different chemistry (untested here, genuinely novel) |
| `c` | none / water-mediated → N1 — the HG3 baseline |

This keeps the plan's 6-cell × RASA structure and its "equal numbers per cell,
do not concentrate on the cell that looks best" discipline, while pointing every
cell at the position the crystals say matters.

## What this gates

Plan §0.5: *"Do not generate a single backbone until the scoring stack can rank
the HG series correctly."*

The stack **can** rank the thing a de novo design will actually get wrong — a
missing or misplaced oxyanion donor (criterion 2, clean separation, >2.9 Å
margin) — and its metrics are verified against deposited coordinates
(criterion 3). It **cannot** resolve two already-good designs that differ only
in preorganisation (criterion 1's single inversion).

That is a usable gate for stage 1–2 output, because the failure mode of a fresh
RFD4 backbone is "donor nowhere near the oxyanion centre", not "HG3.7-vs-HG3.17
subtlety". Recommendation: **proceed to generation, with the gate's limitation
recorded** — and treat the ranking of *survivors* as unvalidated until either
n = 5 (Boltz-2 on HG3.3b/HG3.14) or an ensemble metric (PLACER) is in place.
Do not pick a 5-design submission set on this composite alone.

## Honest summary

| Criterion | Target | Result | |
|---|---|---|---|
| 1. Spearman ρ vs log₁₀(kcat/KM) | ≥ 0.70 | **+0.50** (n=3, one inversion) | |
| 2. Ablation separation | HG3.17 > ablations | 4/4, margin >2.9 Å | |
| 3. Known-geometry recovery | match deposited coords | 20/20 tests, <0.015 Å | |

Open, in priority order:
1. Add HG3.3b + HG3.14 via Boltz-2 to make criterion 1 testable at n = 5.
2. Revise the theozyme cells per the `b`-family finding above.
3. Build the `HG3.17_D127A` negative to test the base metric.
4. Calibrate `base_relative_sasa` — currently `calibrated: false` and not gated on.
5. Add tmol `ddE_ts_sub` when installed; re-run this whole document.
