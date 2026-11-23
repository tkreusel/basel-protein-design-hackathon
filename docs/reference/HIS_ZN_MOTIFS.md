# Alternative base motifs: histidine, and zinc (2026-11-21)

## Histidine as the base

Two theozymes, specs in `inputs/kemp_rfd4_his.yaml`. Conditioning is held
identical to the best Asp recipe (tip atoms, absolute C_CRD, unindexed,
`C_RAS BURIED` on the ligand), so differences in yield come from the motif.

| key | geometry | ND1 -> N3 |
| --- | --- | --- |
| `kemp_his_hg317` | HG3.17 motif, Asp127 -> His (imidazole rigid-fitted on the carboxylate, RMSD 0.09 A), Gln50 donor + Trp44 wall kept | **2.57 A** (Asp OD2 was 2.53) |
| `kemp_his_gdm` | handoff's GDM_KE_1483 Glu -> His (ND1 variant), Ser partner | 3.00 A; Ser OG -> N1 3.08 A |

The handoff's second file (`ND1atOE1`) was **not used**: its nearest ring N is
4.12 A from N3, too far to abstract the proton.

Built with `src/make_his_motif.py`. Selectors were verified to resolve to the
intended atoms (6 imidazole atoms, 3 each Gln/Trp, 12 ligand atoms).

### Chemistry that changes how designs must be judged

- **A neutral imidazole is a weaker base than a buried carboxylate.** Expect
  lower intrinsic activity unless the environment compensates.
- **Only the attacking N accepts.** With ND1 as acceptor, the other ring N (NE2)
  must carry the proton, so it needs an H-bond acceptor partner (tautomer lock).
  Without one the tautomer is not defined.
- **Our base-environment criteria flip.** For Asp we penalise polar atoms near
  the base. For His we must instead (a) require an acceptor near NE2 and (b)
  penalise H-bond donors on ND1, which would occupy the lone pair. The existing
  `base_environment.py` / `steering_criteria.py` assume carboxylate oxygens.
- **Known tight contact:** in the HG3.17 His fit, CE1 sits 2.95 A from Trp44 NE1
  (a C...N contact below vdW distance). The motif is conditioned so the model
  cannot move it; watch for strain or poor yield in `kemp_his_hg317`.

### Not yet generalised

`extract_designs.py`, `geometry.py` and the access/base-environment scripts
hard-code Asp/Gln/Trp identities and carboxylate atom names. They must be
parameterised per motif before the His runs can be gated. Not done yet.

## Zinc as the base: open design decision

Not built. Feasibility facts checked:

- Boltz-2 knows `ZN`, `OH` and `HOH` as components.
- RFD4 benchmark inputs already carry a conditioned `ZN` ion (5 of the 37 MCSA
  enzymes, e.g. M0630_1j79).
- Nothing yet checked on literature precedent for Zn-based Kemp catalysis; do
  not assume it.

A zinc base is almost certainly **Zn-bound hydroxide** acting as the general
base (as in carbonic anhydrase), not the bare ion. That needs:

1. a tetrahedral Zn site (e.g. 3 His, or 2 His + Glu/Asp) with ~2.0-2.2 A
   Zn-N distances, taken from a real structure rather than constructed;
2. an explicit hydroxide/water position as the 4th ligand, oriented at N3;
3. new success criteria: coordination geometry and distances, a free 4th site
   pointing at the ligand, and Zn not exposed.

Boltz-2 will not enforce coordination geometry and cannot model hydroxide, so
the refold test would check only that the three protein ligands bind Zn and the
open site faces the ligand.
