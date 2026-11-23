# Figure guide

Campaign 20261122_kemp_c1. One script per figure plus `kstyle.py`; data are the tables in `$KEMP_RESULTS/tables`; regenerate with `bash make_all.sh`.
Colours (fixed): grey hollow/dashed = unsteered baseline, blue = steered (base/pocket/nitro/motif rewards + TMol), orange = steered + channel (bottleneck 1.7-1.9 A). Circles = hgnt, squares = hgwt.
All figures exclude the beta-sheet diversity runs unless stated. "Evaluated sequences" = LigandMPNN sequences that passed the Boltz-1x gate (strict or fallback) or were added by the extended evaluation, i.e. they have the full suite.

Glossary: ipSAE2 interface score from Boltz PAE (protein->ligand, 0-1, high good). ipDE Boltz interface distance error (A, low good). ligand ipTM Boltz confidence of the ligand placement. dG_sep PyRosetta separated-state binding energy (REU, low good). tmol dE tmol (Rosetta-like) interaction energy difference between ligands in the same relaxed pocket. PLACER apo/holo fraction of a 20-model PLACER ensemble whose catalytic Asp/ligand geometry lies in the crystal-derived window (apo: no ligand, holo: ligand predicted). TSA / substrate / product = 6NT / H5J / 6VP. Shortlist = strict Boltz gate + PLACER (apo>=0.3 or holo>=0.1) + ipSAE2(TSA)>=0.7 + TSA preferred over substrate (z>2 over 10 samples) + dG_sep<=-15 (a proposed definition, uncalibrated).

| fig | file | plotted | read it as | caveat |
|---|---|---|---|---|
| f01 | funnel | counts per 300 generated designs at each stage, per arm and motif (log y) | steepness = attrition; arms above grey = steering helps at that stage | unsteered gated on geometry only; sequence-level stages after the dashed line |
| f02 | boltz_confidence_hist | ipSAE2, ipDE, ligand ipTM, pocket RMSD of the TSA refold | right (ipSAE2/ipTM) or left (ipDE/RMSD) = better | n per arm 10-60, motifs separate |
| f03 | tsa_vs_substrate_product_hist | TSA minus substrate (row 1) / product (row 2) for ipTM, ipSAE2, ipDE, tmol dE, PyRosetta dG_sep | ipTM/ipSAE2 > 0 and ipDE/energies < 0 = TSA preferred | Boltz ipTM is saturated; PyRosetta (independent refolds) shows no TSA preference |
| f04 | energy_placer_hist | dG_sep, tmol strain, PLACER apo and holo fractions | lower energies / higher PLACER fractions better | PLACER ~0.1 binomial uncertainty; windows from HG3 crystals |
| f05 | ipsae_vs_dGsep_scatter | raw ipSAE2 vs dG_sep per sequence; rings = shortlist | bottom-right = best | TMol is a steering reward, PyRosetta is independent |
| f06 | placer_holo_vs_apo_scatter | raw PLACER holo vs apo fractions; rings = strict Boltz gate | upper-right = catalytic geometry in both | fractions are multiples of 0.05 (jittered 0.008) |
| f07 | channel_vs_tmol_designs | every generated design: bottleneck radius vs TMol dG, per arm; dark = pass core | does the channel reward move designs into the 1.7-1.9 A band, and at what energy cost | bottleneck is a discrete probe ladder (jittered) |
| f08 | boltz_gate_criteria | pass rate per Boltz-1x gate criterion and arm | the tall bars are not limiting; ligand RMSD (24 %) and motif RMSD (10 %) are | all folded sequences, both motifs pooled |
| f09 | scaffold_population | length vs helix fraction of core-passing designs | which scaffold families the arms produce | passing designs only |