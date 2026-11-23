# 01 — Backbones, sequence design and refolding (stages 1–5)

`inputs/` RFD4 specs for the theozyme: `kemp_rfd4.yaml` (basic), `kemp_rfd4_full_buried.yaml` (whole motif + BURIED), `kemp_rfd4_tips_buried.yaml` (tip atoms + BURIED),
`kemp_rfd4_his.yaml` (His base), plus motif structures (`HG3_17_kemp_motif.pdb`, `kemp_motif.cif`, His/GDM variants) and `maxi/` (the steering campaign's with-Trp spec for the control run).

`jobs/`

| Job | Does | Needs |
|---|---|---|
| `stage1_rfd4.sbatch` | RFD4 backbone + sequence co-design for one spec | `KEMP_SPEC`, `KEMP_N`, `KEMP_BS`, `KEMP_SEED` |
| `stage2_ligandmpnn.sbatch` | LigandMPNN re-design of a finalized stage-1 run, catalytic residues fixed and verified | `KEMP_RUN`, `KEMP_NSEQ` |
| `refold_boltz.sbatch` / `refold_shard.sbatch` | Boltz-2 refold with 6NT of the stage-1 successes (single shard / sharded across two small GPUs) | `KEMP_RUN`, `KEMP_TAG` |
| `stage2_5_design_refold.sbatch` | both design arms (RFD4 sequences and LigandMPNN) refolded like-for-like | `KEMP_STAGE1` |
| `msa_ab_test.sbatch` | with-MSA vs no-MSA Boltz comparison (de novo designs fold identically) | |
| `rfd4_smoke.sbatch`, `probe_gpu.sbatch` | environment checks (GPU visible? RFD4 runs?) | |

Results of our batches: 12/32 (spec without burial), 31/32 (full motif + BURIED), 16/32 (tip atoms + BURIED), 3/32 (His). See `docs/06_RESULTS_AND_FUNNELS.md` §3.
