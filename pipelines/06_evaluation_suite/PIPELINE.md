# Kemp eliminase in-silico evaluation pipeline (Maxi)

Preorganisation-based scoring of Kemp eliminases bound to a transition-state analogue (6NT) using
**PLACER** (ensembles), **ipSAE2** (Boltz-2 interface confidence) and **tmol** (Rosetta beta2016 in PyTorch).
Location: `/mnt/work/projects/Maxi/pipeline/`. **No Slurm job has been submitted.**

## 1. Status

| Piece | State | How it was checked |
|---|---|---|
| Input standardisation (crystal / NMR / Boltz → one format) | done | all 5 Maxi files, 18 ligand copies/models |
| Ligands: **6NT (TSA), H5J (substrate), 6VP (product)** as CCD codes with CCD atom names, + tmol parameters | done | role atoms checked against the CCD; isomer H5J/6VP never confused by the matcher; tmol builds all three with no chemistry warnings |
| Geometry (base, syn/anti, oxyanion donors, wall, burial, water-site proxy) | done | reproduces known crystal values; independent heavy-atom check of syn/anti |
| Theozyme reference constraints (from the crystals) | done | `config/theozyme_reference.yaml` |
| tmol protocol (relax, constraints, ΔΔE, repulsion, field, BUNS/void, apo, ablations) | done | CPU integration test, see §6 |
| tmol constraint term | done | energy equals numpy deviation to 4 decimals, before and after relax |
| PLACER holo + apo, ensemble analysis | done | CPU smoke test, **2 samples only** (format + code path) |
| ipSAE2 + site-focused Boltz confidence | done | run on an existing Boltz refold; matches the standalone ipsae value |
| **PyRosetta 2026.47** (own env `/mnt/work/envs/pyrosetta`, CPU only) | installed, imports, filters present | free for non-commercial academic use under the RosettaCommons licence: please make sure your use qualifies |
| Rosetta params for 6NT / H5J / 6VP | done | CCD atom names + the same MMFF charges as tmol |
| PyRosetta stage (relax, interface dG, Sc, packstat, BUNS, independent ΔΔG, ablations) | done, tested on one crystal copy (CPU, 6 min / 18 min with ablations) | relaxed score −3.0/res; ΔΔG signs agree with tmol. **Sc is unreliable** (0.78 unrelaxed vs 0.38 relaxed for a ligand that moved ~0.3 Å); treat as unvalidated |
| Boltz single-chain refolds (5) registered + ipSAE2 + geometry quick look | done (CPU) | `PROPOSAL.md` |
| Aggregation + report | written, merges all four tools | exercised only on single-replicate outputs |
| Rosetta `.cst` reader | written, **unverified against your file** | round-trip on a synthetic `.cst` only |
| Native substrate/product refold glue (`run_tmol_entry`) | written, **untested** | no Boltz substrate/product refolds exist yet |
| Any GPU run, any run with ≥ 3 PLACER samples | **not done** | head node has no GPU; you asked not to submit |

## 2. What I need from you (blocking questions are marked ⛔)

1. ⛔ **Your theozyme CST file** (or the QM TS geometry). I found none on disk. Until then the constraint set is
   *derived from the crystal TSA complexes* (defensible, but it is not your theozyme). The `.cst` reader follows the
   documented enzdes block format; it needs the real file to confirm atom names/types.
2. **Boltz refolds for the three ligands** (being generated): every protein is run with **6NT, H5J and 6VP** (CCD codes). The arm of each
   model is read from the ligand actually present in the CIF, so file names only need to carry the protein: `<entry>_<6NT|H5J|6VP>`,
   `<entry>__<tsa|sub|prod>` or just `<entry>` all work (case-insensitive match to `config/entries.yaml`, e.g. `7k4q_H5J` -> `7K4Q`).
   Refold entries get the suffix `_boltz` so crystal and refold results stay separate. Set `BOLTZ_ROOT` to the output folder.
   A name that disagrees with its ligand is reported and the ligand wins.
3. ~~PyRosetta?~~ Installed (§4). Still open: **your `.cst` file** is needed to use enzdes scoring (§2.1); REMARK 666 mapping of CST blocks to residues
   is currently "block 1 = catalytic base" and must be confirmed against the real file.
4. **Activity labels.** Only HG3 = 146 M⁻¹s⁻¹ is known for the Maxi set. A rank test needs ≥ 4 distinct values.
   `config/entries.yaml` has an opt-in switch for 7K4X (27,000) and 7K4Z (126,000) from the sibling project; say if you
   want them on, and send values for 8RD5 / 4BS0 / 29SB.
5. NMR entry 29SB (10 models): run 3 of them, or all?  (Default `--max-replicates 3`.)

## 3. The pipeline

```
Boltz-2 refolds (yours) + Maxi crystals
        │  00_prepare_and_ipsae.sh  (CPU, minutes, head node)
        ▼
 work/inputs/manifest.json  ── standardised PDBs (protein chain A, ligand chain B, canonical atom names, no H)
        ├─► ipSAE2 (CPU)       protein→ligand ipSAE, ipTM, site PAE/pLDDT per Boltz model; sd over models
        ├─► 10_placer.sbatch   (GPU) PLACER holo + apo ensembles, 50 models each
        ├─► 20_tmol.sbatch     (GPU) tmol protocol per replicate
        ├─► 30_placer_analysis.sbatch (GPU) ensemble statistics + tmol electric field per model
        └─► 40_pyrosetta.sbatch (CPU only, no GPU reserved) Rosetta relax / interface / Sc / packstat / BUNS / independent ddG
                 ▼
          aggregate (CPU) → work/report/{summary_entries.csv, summary_replicates.csv, report.md}
```
Four environments, one per tool (they cannot be merged): RFD4 `.venv-cu12` (tmol + analysis), `placer_env`, `ipsae2`, `/mnt/work/envs/pyrosetta`.
`jobs/submit_all.sh` prints the chain (dry run); GPU jobs are chained, not parallel, because of the 8-CPU team budget.

### tmol protocol per TSA replicate
| tag | what | answers |
|---|---|---|
| `raw` | as given, tmol protonation | baseline |
| `rlx` | full Cartesian FastRelax, CA restrained (0.5 Å), ligand free | does the site hold by itself? |
| `rlxcst` | same, theozyme constraints held (weight ×10) | `strain.E_noncst` = energy cost of imposing the theozyme |
| `tsa_pk`, `sub`, `prod` | substrate / product superposed on the TSA in the **same** relaxed pocket, 8 Å pocket relax | `delta.*`: TS-stabilisation ΔΔE = ΔE(TSA) − ΔE(substrate); product binding |
| `sub_ts`, `prod_ts` | ligand held at the TSA pose, only side chains repack | charge/shape at the TS geometry |
| `native_sub/prod` | your independent Boltz substrate/product refolds, same full relax | second, independent ΔΔE (`delta_native.*`) |
| `apo` (opt.) | protein-only relax | reorganisation energy |
| `abl.*` (opt.) | base / donor / wall → Ala | stability change and loss of ligand interaction |

### Three-ligand inputs (6NT / H5J / 6VP)
| Stage | What runs per arm |
|---|---|
| prepare | arm from the ligand in the CIF; catalytic site detected per arm |
| ipSAE2 | every arm separately (`ipsae.tsa`, `ipsae.sub`, `ipsae.prod`); `ipsae.delta_tsa_sub`; site PAE/pLDDT |
| PLACER | **holo ensemble for each ligand**; the apo ensemble once per entry (it does not depend on the ligand) |
| tmol | full protocol on the 6NT refold; the H5J and 6VP refolds get the same full relax + energies + geometry (`native_sub.*`, `native_prod.*`), giving **independent** `delta_native.*` next to the same-pocket-swap `delta.*` |
| PyRosetta | same: `ros.native_sub.*`, `ros.native_prod.*`, `ros.delta_native.*` |
| geometry | substrate arm gets the proton-transfer geometry (base to C3-H2); the product has no abstractable proton so it gets donor/wall/burial metrics only |

## 4. Your tmol list: what is native, what is a proxy, what does not make sense

| Request | Status | Notes |
|---|---|---|
| Per-constraint deviations (O···H, O···H–C3 linearity, syn/anti, dihedrals) | **native-style, implemented** | numpy on the relaxed pose with tmol hydrogens, plus tmol's own constraint term (verified equal). Reported as dev, z = dev/tol, pass/fail. Your `.cst`: reader written, unverified. |
| H-bond donors to the leaving-group atom | **implemented** | every donor within 3.6 Å, best side-chain donor, real D–H···A angle, tmol `hbond` energy per residue to the ligand. For 6NT the "phenolate O" is N1 (neutral, so weaker than the real TS). |
| ΔΔG(TS) − ΔΔG(substrate) | **implemented as constrained interaction-energy difference** | same pocket, same protocol, plus the independent-refold variant. Ranking heuristic only (your caveat is right). |
| Product binding / product inhibition | **implemented** | `prod` free and `prod_ts`; flag when product binds tighter than substrate. Product is typed by tmol as an alcohol-like O (`Oal`), charge −0.83: approximate. |
| Repulsion terms | **implemented** | `fa_ljrep` split into intra-residue (block diagonal) and inter-residue, whole protein, pocket and ligand–protein. tmol's intra term is *not identical* to Rosetta `fa_intra_rep`. Cheap; do not over-trust (16-active benchmark). |
| Total score per residue | **native** | plus per-residue energies and strained-pocket-residue counts |
| Buried unsatisfied polars | **proxy** | geometric BUNS from tmol hydrogens + SASA; not the Rosetta filter |
| Holes / packstat | **not in tmol → proxy** | enclosed empty volume next to the ligand (grid probe). Real PackStat needs PyRosetta. |
| Shape complementarity to the TS | **not in tmol, not implemented** | reported instead: ligand rel-SASA, enclosure, contacts (same definitions as the tobi_v1 project). Real Lawrence–Colman Sc needs Rosetta/CCP4. |
| `cartesian_ddg` point mutants | **approximation** | ablations X→Ala (pack + relax of the pocket); *not* the Rosetta protocol. ~3 relaxes per replicate. |
| Electric-field projection on C–H and C–O | **implemented, unscreened** | tmol partial charges, vacuum Coulomb (+ ε = 4r check), environment without the base, MV/cm, positive = favourable. PBSA/APBS not installed; use for ranking and ensemble variance only. |
| Field variance across the ensemble | **implemented** | per PLACER model (top-N by prmsd) |
| Water: H-bonds to the carboxylate, trapped waters | **proxy / not possible** | tmol has no explicit water. Implemented: lone-pair water-site accessibility (agrees with crystal waters in 7 of 8 copies, a small, partly non-independent check) and tmol's implicit `lk_bridge`/`lk_ball` terms. 3D-RISM/GIST/Dowser/SPaDES would be separate software; high-energy trapped waters cannot be assessed here. |

### With PyRosetta installed, the not-in-tmol items become:
| Item | PyRosetta implementation | Status |
|---|---|---|
| Shape complementarity to the ligand | Lawrence–Colman Sc via `ShapeComplementarityCalculator` | **real Sc**, but Rosetta ignores non-polymer residues, so ligand atoms are added by hand with element radii (found: filter returns −1/0 otherwise). 4BS0: 0.78 |
| Packstat | `InterfaceAnalyzer` interface packstat + `PackStat` filter | real; whole-pocket-local packstat is not available |
| Holes | **not possible**: needs the external `DAlphaBall` binary, absent from the wheel (and its absence kills the process) | proxy `void.volume_A3` stays |
| Buried unsatisfied H-bonds | `BuriedUnsatHbonds` filter (pocket, ligand) + interface delta-unsat | real Rosetta filter; compare with the tmol-based approximation |
| cartesian ddG | own cartesian pocket relax of X→Ala, mean of `nstruct` | still an approximation of the `cartesian_ddg` app (not used) |
| Binding energy | `InterfaceAnalyzer` dG_separated (repack-separated) for TSA / substrate / product in the same pocket | **independent cross-check of the tmol ΔΔE** |
| enzdes CST energies | `AddOrRemoveMatchCsts` (needs REMARK 666 block→residue mapping) | written; **untested** (needs your CST) |

## 5. Things found along the way that change how to read results

* **HG3-lineage crystals bind the TSA through the anti lone pair of Asp127**, not syn: O′–C–O···N3 ≈ 150–170°
  (4BS0, 7K4Q, 8RD5; confirmed with a heavy-atom-only calculation, independent of hydrogen placement). Syn/anti is
  reported, but I did **not** gate on "syn is better". If your CST demands syn, the HG3 crystals would fail it.
  (A neutral N–H azole TSA is not the true TS, so this is evidence, not proof.)
* The O···H–X torsions (about the H-bond axis) vary by ±40–55° between crystal copies, so they carry little
  discrimination and have wide tolerances.
* KABLE2.5 (29SB, NMR): Asp49 is the base, there is no donor to N1 within 3.6 Å, and the O···H–N angle is ~110°
  (not linear). It is kept as a separate reference group because its scaffold differs.
* HG3 (7K4Q, 5RGA) has no donor to N1 except the weak Cys/Met84 backbone NH; HG3.17-type (4BS0, 8RD5) has Gln50 NE2
  at 3.0–3.1 Å, as expected.
* **Ligand naming:** an early version used stand-in codes (`SUB`/`PRD`, then `ZKB`/`ZCN`) with my own atom names; for real CCD codes
  tmol/atomworks merged the CCD definition with mine and silently capped/re-protonated atoms. The pipeline now uses the genuine CCD codes
  **6NT, H5J, 6VP with their own CCD atom names**, so every tool sees one definition. **6VP in the CCD is the neutral phenol**
  (5-nitro-2-hydroxybenzonitrile), not the phenolate anion: product-binding energies therefore under-represent the charged product.
  H5J and 6VP are isomers (both C7H4N2O3); the atom matcher distinguishes them by connectivity and is unit-tested for that.
* tmol 0.2.0 `fast_relax` leaves the constraint weight at 0 afterwards (energy reports then show 0 constraint energy). Fixed in
  `relax()`, covered by a test.
* A tmol context only knows ligand types present when it was built; the session now builds from the protein plus all three ligands.
* PyRosetta install: the official installer calls whichever `pip` is first on `PATH` (it used the system one and failed with
  "not a supported wheel"); it works with the env's `bin` first. Needs `-corrections::beta_nov16 true` for `beta_nov16`.
  InterfaceAnalyzer reads cached energies (cross-term = 0.0 on freshly modified poses unless the pose is scored first; fixed).
* **Bug found and fixed in my PyRosetta stage:** `FastRelax(constrain_relax_to_start_coords)` leaves a coordinate-restraint weight on the score
  function it is given, adding ~5,200 REU to every later score (relaxed structure looked worse than raw: +5.7/res vs the true −3.0). Relax now runs on a clone.
  Energies from the *first* PyRosetta test (before the fix) were contaminated and were discarded.
* `InterfaceAnalyzer` cross-term energy stayed 0.0 even after scoring first: dropped from the output; use `dG_sep`.
* **tmol relax is stochastic:** two identical runs differed by several REU (ΔΔE TSA−substrate −5.2 vs −3.6; "strain" +3.5 vs −12). Treat
  differences below ~5 REU as noise; the second crystal copy of each entry gives a first noise estimate.
* Boltz refolds of the five native sequences recover the site (base–N3 2.5–2.8 Å, Gln at 2.8–3.0 Å in the two Gln50-type enzymes, no donor in HG3);
  ipSAE2 p→l 0.85–0.90, ligand ipTM 0.96–0.98, catalytic-residue pLDDT 0.98–0.99. These are native sequences of structures Boltz may have seen, so
  this is a pipeline test, not evidence. For KABLE2.5 (29SB) Boltz places a Gln75 donor (3.04 Å) and a *syn* approach (−10°) that the NMR ensemble does not show.
* PLACER: only 1–3 s/model on a GPU; on the head node ~80–100 s/model, so run it on a GPU node. PLACER warns about nitro
  charge-state differences between the PDB and its CCD ligand; harmless here. With 2 samples the holo ligand landed 3–4 Å
  off the input pose: meaningless at n = 2, but the 50-model ensemble will show how much PLACER re-docks the TSA.

## 6. Caveats (please read before trusting any number)

* 6NT is a **neutral** TSA, not the delocalised-anion TS; Rosetta energies are not free energies. Deltas rank, they do not predict kcat/KM.
* Relaxation is a local minimisation: "substrate free" means refined near the TS pose, not a Michaelis-complex search.
  The independent substrate refold is the better Michaelis proxy.
* Reference constraints come from **8 HG3-lineage crystal copies** (tolerance = max(2 sd, floor)); no activity calibration exists
  for any composite score. Do not combine into one number until the validation ladder has been run.
* Boltz and PLACER outputs are not experimental structures; the pipeline measures geometry *of the predicted pose*.
* Head-node timings are CPU; GPU timings are unmeasured.

## 7. Run it (nothing here has been submitted)

```bash
cd /mnt/work/projects/Maxi/pipeline
BOLTZ_ROOT=/path/to/boltz/outputs bash jobs/00_prepare_and_ipsae.sh     # CPU, minutes
bash jobs/submit_all.sh                                                   # dry run: prints the sbatch chain
bash jobs/submit_all.sh --submit                                          # only when you say so
$PY_TMOL -m maxi_eval.aggregate                                           # after the jobs
$PY_TMOL tests/test_unit.py                                               # unit tests
```
Options: `ENTRIES=4BS0,7K4Q`, `MAXREP=3`, `NSAMPLES=50`, `CST=/path/theozyme.cst`, `EXTRA="--apo --ablate"`.

## 8. Full-dataset run (many proteins x 6NT / H5J / 6VP) and discriminative power

```bash
bash jobs/run_full_suite.sh --root /mnt/work/projects/Maxi/boltz_outputs            # DRY RUN: readiness, queue, plan, cost
bash jobs/run_full_suite.sh --root ... --submit                                    # after sign-off: 1 PLACER GPU job + head-node CPU shards + finisher
python -m maxi_eval.check_boltz_ready --root ...                                    # which proteins lack an arm / PAE files
```
* PLACER runs on enzyme-TSA, enzyme-substrate **and** enzyme-product complexes (holo) and once per protein in apo mode; ipSAE2 scores all three arms.
* tmol and PyRosetta run the TSA protocol plus independent H5J/6VP evaluation (`delta_native.*`). Ablations are off in the full run (`--ablate` to add them).
* A finisher process waits for the GPU job and CPU shards, then runs PLACER analysis, aggregate, validate and `validate_dataset`; it writes `work/report/FULL_SUITE_DONE`.
* **Labels:** `config/labels.csv` (`entry,kcat_KM[,class,note]`, entry without `_boltz`). No labels -> no discriminative-power claim.
* **Discriminative power** (`python -m maxi_eval.validate_dataset --group boltz`): per protein (replicates averaged), Spearman vs log10 kcat/KM with
  bootstrap CI and Benjamini-Hochberg q over the ~40 pre-declared metrics (`--panel all` = exploratory), AUROC for active vs inactive, partial rho controlling
  for oxyanion-donor presence, and a **leave-one-out composite** with metric selection inside the loop (the honest estimate). Refuses to rank below 12 labelled proteins.
  Unit-tested on simulated data (planted signal recovered; selection from 120 noise metrics does not look predictive when held out).
