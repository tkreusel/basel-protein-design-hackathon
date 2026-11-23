# 08 — Eligibility gate (novelty)

Requirement: < 70 % sequence identity to published enzymes and TM-score ≤ 0.7 to published structures. All local; no sequence or structure leaves the cluster.

| Stage | Tool | Database | Output |
|---|---|---|---|
| 1 length/alphabet | — | — | 50–260 aa, 20 standard letters |
| 2 PDB identity | MMseqs2 `easy-search` (`-s 7.5 --max-seqs 2000 -e 1e-3`) | PDB seqres, protein chains | global identity = identity × aligned length / query length, plus aligned-region identity with coverage and ≥ 90 % stretches |
| 3 UniRef90 identity | MMseqs2 GPU search (`uniref_search.sbatch`, 30 GB RAM → node `l4-1`) | UniRef90-2022-05 | same columns; TSV is read by stage 4's script |
| 4 TM-score | Foldseek `easy-search` (TM-align hit order) then **US-align** on the top 20 distinct PDB entries, both normalisations (max reported) | PDB 30 %-cluster representatives (40,064 entries, built locally) | `tm_max`, `tm_hit` |
| verdict | | | FAIL if identity ≥ 0.70, aligned identity ≥ 0.70 at ≥ 50 % coverage, or TM > 0.70; WARN within 10 % of a limit; **control = published HG3.17 must FAIL** |

Files: `check_eligibility.py` (flags `--reuse --out --uniref-tsv --threads --skip-structure`), `make_candidates.py` / `make_candidates_v6.py` (build `candidates.csv`: `cand_id,sequence,structure`,
**absolute** structure paths), `uniref_search.sbatch`, `fetch_one.sh` (one PDB mmCIF from RCSB), `db_recipe/pdb30_entries.txt` (entry list for the structure DB),
`results/` (candidate lists v5/v6 and final `eligibility_v5/v6.csv`).

Join eligibility results to designs **by sequence**, not by name (50-character truncation collided ids once). See `docs/07_LESSONS_AND_PITFALLS.md`.
