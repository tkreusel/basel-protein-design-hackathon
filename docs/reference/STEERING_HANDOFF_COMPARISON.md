# Steering handoff vs our pipeline (2026-11-21)

Source: `20261121_rfd4_steering_handoff.tar.gz` (Maxi, bwForCluster Helix),
unpacked to `imports/steering_handoff_20261121/` (treated as data; none of its
code was executed). Its README is `README_AGENT_HANDOFF.md`.

## What they did

RFD4 with **inference-time beam search** (4 branches x beam 3; checkpoints at
denoising steps 30/42/52/60, rewards computed on the predicted clean
structure). Rewards are chemistry descriptors of the active site. Same model
(RFD4-Proteina), same motif (HG3.17 Asp127 / Gln50 / Trp44 + 6NT).

Conditioning in every arm (steered or not): tip atoms fixed, `C_RAS BURIED`
(motif + ligand), `C_HPHOB` (base CB/CG, ligand ring C), `C_HB ACCEPTOR` on the
nitro oxygens.

Reward terms: base polar side-chain atoms within 4.5 A (weight -4), Lys/Arg
within 6 A of the base (-3), base SASA (-2), base apolar fraction, pocket
volume band, ligand SASA, nitro H-bond far from base (+) / near base (-),
clashes, motif RMSD (base, partner), solvent-channel band 1.4-1.9 A (channel
arm), TMol dG (tmol arm).

Downstream (built, **not yet run on real data**): LigandMPNN with
K R H D E N Q S T Y C omitted within 6 A of the base carboxylate, Boltz-2 refold.

## Their results (279 designs, 9 lineages per steered arm)

Conditioning-only arm: 0 of 18 designs pass their core criteria per motif;
base polar atoms 1.8-2.4, cations 0.9-1.4. Steered arms: 3-12 of 27 designs
pass, cations ~0, polar ~1.0-1.1 (the fixed Trp NE1). Only 3 designs pass core +
channel >= 1.4 A. Their own caveat: pass criteria equal the reward terms
(circular); no refold results yet.

## Our designs under their criteria

Their base-environment metrics were reimplemented from the definitions
(`src/steering_criteria.py`) and **validated on their own files: 198/198
HG3.17-motif designs reproduce their published values exactly** (base SASA,
polar atoms, cations, nitro H-bonds).

Stage-1 designs (partial core: ligand SASA, base SASA, polar <= 1, no cation,
nitro H-bond, no clash; volume ratio, motif RMSD and TMol not applied):

| our run | n | our gate | polar<=1 | no cation | nitro ok | all |
| --- | --- | --- | --- | --- | --- | --- |
| full motif, no burial | 32 | 12 | 9 | 23 | 8 | 0 |
| full motif + BURIED | 32 | 31 | **0** | 17 | 15 | **0** |
| tip atoms + BURIED | 32 | 16 | 13 | 22 | 9 | 3 |

Mean base polar atoms: 2.28 / **3.19** / 2.00. Burial conditioning alone
raised polar atoms around the base (2.28 -> 3.19).

Refolded shortlist (31): ligand SASA 29, base SASA 25, polar<=1 **8**, no cation
18, nitro 15, no clash 24, **all six: 0**; polar<=1 and no cation: 6.

## Interpretation

1. Our gate measured burial and catalytic distances. Their criteria measure the
   chemistry around the base. They are independent: 31/32 of our job-21
   "successes" have 0 designs with a clean base environment.
2. Steering demonstrably moves those descriptors (their data), so the
   cation/polar neighbours we found around the base (42 of 71 promising
   designs had a Lys/Arg within 6 A) are fixable at generation time, and
   cheaply at sequence time with their omit-list.
3. We have what they lack: a refold-based, non-circular validation funnel
   (fold, pocket, constellation, donor), scRMSD, ipSAE, and the finding that
   tip-atom conditioning plus LigandMPNN gives 16/64 strict passes vs 4/124 for
   full motif. Their downstream is untested.
4. Convergent results: both found tip-atom conditioning is the better recipe,
   both measure a 1.4 A channel bottleneck, both find real crystals have 0 cations
   near the base.

## Caveats

- Their pass criteria are partly calibrated on the crystals and equal the rewards:
  passing is not evidence of activity.
- n is small and lineages (3 designs each) are the real sample size.
- Our check is on RFD4 designs (stage 1) and Boltz refolds (shortlist); side
  chains in refolds are exactly the thing in doubt.
