# Final top-20 design structures (eligibility-filtered, current ranking)
- `boltz_6NT_refold/<slot>.cif` Boltz-2 refold with 6NT (what the hosts score): protein chain A, ligand chain B. Maxi overnight designs: his best of 10 TSA samples; others: our seed-0 refold.
- `rfd4_design/<slot>.cif` the RFD4 backbone + ligand the sequence was designed on (missing where not kept; see index.csv).
- Slot names = `final_designs.fasta` (SubstrateVentures01..20); `index.csv` maps slot -> design id, base/oxy residue (1-based), distances.
Catalytic residues: base Asp/Glu = base_res, oxyanion donor = oxy_res. Superpose on chain A CA atoms.
Regenerate after ranking changes: python -I tobi_v1/submission/export_final_structures.py
