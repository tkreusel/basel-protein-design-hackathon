# Kemp motif scaffolding (2026-11-21)
Motif: HG3.17 A44 (Trp ring), A50 (Gln side chain), A127 (Asp carboxylate) + ligand 6NT (TS analogue), all unindexed.
Pipeline: RFD3/RFD4 -> LigandMPNN (fixed motif) -> Boltz-2 (+ AF3-class refold) -> PLACER (preorganisation) -> filters.
Files: rfd3/inputs/kemp_rfd3.json (4 arms), rfd3/run_rfd3.sbatch (devel test), rfd3/check_outputs.py.
Env: fkguid_boltz WITHOUT PYTHONNOUSERSITE (torch comes from ~/.local); foundry env is broken.
Queue: partition devel (30 min max, A40) starts in ~1.5-2 min; gpu-single estimated start 2026-12-12.
