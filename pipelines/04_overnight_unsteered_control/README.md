# Overnight unsteered control (prepared 2026-11-21, NOT submitted)

Unsteered RFD4 control for Maxi's steered run, mirroring his handoff
(`imports/steering_handoff_20261121/`). One idempotent Slurm job:

```
GEN -> SELECT -> MPNN -> BOLTZ -> SCORE        sbatch overnight/run_overnight.sbatch
```

| Stage | What | Tested |
| --- | --- | --- |
| GEN | `rfd4-generate`, his with-Trp spec, conditioning only (`~metrics`, `batch_size=4`, `compile_lengths=null`) | spec resolves, overrides compose; **GPU path untested** (GPU was busy) |
| SELECT | his core criteria (7 of 10 computable), top-K | on his 18 unsteered designs: 0 pass all, as he published |
| MPNN | his settings; omit `KRHDENQSTYC` within 6 A of the base O; one process for all designs | 6 designs, CPU, 24 sequences, 0 violations |
| BOLTZ | no MSA, 6NT, RFD4 + 4 MPNN sequences per design, chunks of 125 | inputs generated; Boltz command identical to our earlier runs |
| SCORE | fold, site, contacts, base environment, access, ipSAE | reproduces earlier analysis exactly on 10/10 designs |

## Settings: mirrored vs assumed

Mirrored from his scripts: spec, `batch_size=4`, `compile_lengths=null`,
`~metrics`, LigandMPNN (v_32_010_25, T=0.1, 4 seqs, seed 111, atom context only),
omit list and 6 A zone, Boltz (no MSA, 3 recycles, 1 sample).

**Decided by Tobi (2026-11-21): use his settings for the rest, top 100.**

| Setting | Used |
| --- | --- |
| seed | 7 (his first seed; his unsteered arm used 7, 11, 12) |
| recycle strategy | repo default (`thinking_medium`), as his `run_r5.py` |
| Boltz seed | 0 (his `run_boltz_chunk.sbatch`) |
| designs going downstream | top 100 of 250, ranked by number of his criteria passed, NO pass/fail cut-off |
| criteria not computable | pocket volume ratio, pp clash, TMol dG |

The criteria do not steer anything: generation is conditioning only. All 250
designs are measured and saved (`selection/all_designs_metrics.csv`), so the full
unsteered distribution is reported; the criteria only rank which 100 are refolded.

## Not applied

His full `pass_core` needs custom metrics and TMol. We apply 7 of its 10
criteria, so pass counts are NOT comparable to his `pass_core` totals.
His evaluation-suite owner (mnt-ea) confirmed the suite does not patch the RFD4
clone and is meant for finalists only.

## Resource conflicts

Team quota is one GPU at a time. mnt-ea needs ~40 min of one L4 for a PLACER job;
that ordering is the user's call.

## Env vars

`KEMP_N=250 KEMP_K=100 KEMP_MODE=topk|all|core KEMP_SEED=201 KEMP_BS=4
KEMP_RECYCLE= KEMP_BOLTZ_SEED=0 KEMP_RUN=<dir to resume>`
