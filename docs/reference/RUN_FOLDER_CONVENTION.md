# runs/ — how to read this folder

Written for a human opening this cold. You should not need to read any code to
find out what happened.

## Folder naming

```
<date>_<stage>_job<slurm id>_<size>_<extra>
2026-11-21_rfd4-backbones_job13_32designs_seed101
|          |                |     |         |
date       which stage      Slurm job        the seed, so it is reproducible
                                  how many designs were attempted
```

The Slurm job id is in the name so a run folder pairs with its log:
`runs/..._job13_.../` ↔ `logs/stage1-13.log`.

## Inside every run folder — always the same shape

| Path | What it is |
| --- | --- |
| `SUCCESS.md` | **Read this first.** Defines what success means for this stage, with the numeric thresholds and what they were calibrated against. Also states what success does *not* mean. |
| `summary.csv` | One row per attempt: every measurement, the verdict, and the failure reason. Open in a spreadsheet. |
| `attempts/` | **Every** attempt. The verdict is in the filename, so `ls` alone tells you what happened. |
| `success/` | Only the successes, copied and **ranked best-first**. `rank01_*` is the best. Plus `sequences.fasta`. |
| `manifest.jsonl` | Machine-readable, one DesignRecord per line. For scripts, not people. |
| `designs/<id>/record.json` | Per-design record. For scripts. |
| `backbones_pdb/` | PDB copies of every attempt, the hand-off to LigandMPNN. `manifest.jsonl` points here, so do not move it. For scripts. |
| `armA_rfd4_sequences.fasta` | Successes only, consumed by the next stage. Same content as `success/sequences.fasta` but without the ranking in the header. For scripts. |

**For people: `SUCCESS.md`, `summary.csv`, `attempts/`, `success/`.** The rest
is plumbing.

### `attempts/` filenames carry the verdict

```
attempt04_b8185dadc4d0_SUCCESS.cif
attempt07_1d45841068fb_FAIL-ligand-on-surface.cif
attempt11_af14b3c83c24_FAIL-ligand-on-surface.cif
```

Failure reasons used: `ligand-on-surface`, `base-too-far`, `donor-too-far`,
`wall-too-far`, `too-long`, `clash`, `gate-fail`, `ERROR` (crashed).

## The three top-level folders

- **`<date>_<stage>_...`** — real runs. These are results.
- **`_superseded/`** — runs kept only as evidence of a bug or a superseded
  method. **Not results.** Each carries a `WHY_SUPERSEDED.md`. The folder name
  states the problem, e.g.
  `..._SPEC-BROKEN-no-coordinate-conditioning`.
- **`_smoke-tests/`** — environment and installation checks. No designs worth
  keeping.

## Success definitions by stage

Each stage writes its own `SUCCESS.md`, but in short:

| Stage | Success means |
| --- | --- |
| `rfd4-backbones` (stage 1) | **Both**: (1) the active site survived — base, oxyanion donor and aromatic wall all at calibrated catalytic distance to the ligand; **and** (2) the ligand is in the right place — buried in a pocket, not on the surface. |
| `stage2-5` (sequence + refold) | not yet defined — will be: the active site survives an independent Boltz-2 refold from sequence alone |

Criterion (2) for stage 1 is not cosmetic. In the first real batch, **20 of 32
designs satisfied every catalytic distance while the ligand sat on the outside
of the protein**, up to 54% solvent-exposed. A distance-only check passed all
of them. Reference crystals have the ligand 99% buried.

## Generating a run folder

`src/finalize_run.py` creates this structure. It never judges anything itself —
it only lays out what `src/extract_designs.py::apply_motif_gate` already
decided, so `summary.csv` can never disagree with the manifest.

```bash
source env/paths.sh
env/kemp-orch/bin/python -I src/finalize_run.py . \
    <raw rfd4 out_dir> runs/<date>_<stage>_job<id>_<size>
```

## Where the raw model output goes

`finalize_run.py` reads the raw RFD4 output but does not move it; `summary.csv`
has a `cif` column pointing at the original file for provenance. Raw output for
the current run lives under the run folder's own tree where RFD4 wrote it.
