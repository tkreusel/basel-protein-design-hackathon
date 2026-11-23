# 11 — Selection and submission

| Script | Role |
|---|---|
| `merged80_pool.py` | merge Maxi's overnight top 30 (`maxi_top30_scored.csv`) with all our scored candidates (`ranking_v2.csv`), de-duplicate by sequence, top 80 by judging score → `merged80_pool.csv` |
| `rerank_merged80.py` | join eligibility by sequence (drop FAIL, TM > 0.70, UniRef90 identity ≥ 0.70), literature factor L (common 8 metrics), final = judging × L, de-duplicate sequences → `rerank_all.csv`, `rerank_top20.csv`, `designs_rerank_top20.fasta` |
| `lit_calibrated_rank.py`, `build_final_table_lit.py`, `build_final_table.py` | earlier 12-metric version and table builders (kept for provenance) |
| `export_final_structures.py` | copies the Boltz refold and RFD4 backbone of the final 20 |
| `validate-submission.py`, `submit-designs.sh` | the hosts' validator and uploader (from the organisers' bucket) |
| `results/` | `designs.fasta` (as submitted), `final_designs_name_mapping.csv` (slot → source design, steering, sequence origin), ranking tables |

Final names: `SubstrateVentures01…20_<steered|unsteered>_<mpnnN|rfd4|shell1>` (rank order); header `base=` / `oxy=` are 1-based residue numbers in the submitted sequence.
Weights and formulas: `docs/05_SCORING_AND_RANKING.md`.
