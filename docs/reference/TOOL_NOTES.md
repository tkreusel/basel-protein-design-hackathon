# Tool notes

Observed interfaces only: we did not guess CLI flags; we cloned the tool, read
`--help`, ran the bundled example and recorded what happened. A tool that
will not run gets its failure recorded and its fallback chosen, not a retry loop.


---

## RFD4-Proteina — installed, GPU blocked, fix identified

- **Location:** `/mnt/work/repos/RFD4-Proteina`, git `6306a21`, version `0.0.3`.
- **Entry points:** `rfd4-generate`, `rfd4-evaluate`, `rfd4-train`,
  `rfd4-train-ae` in `.venv/bin`.
- **Weights:** HF bundle `nvidia/RFD4-Proteina`, 1.7 GB, cached under
  `$HOME/.cache/huggingface` by the original clone. `configs/inference.yaml`
  selects the release checkpoint pair from the Hub, so a default run needs no
  checkpoint override.
- The repository ships its own contributor notes (`AGENTS.md`) and task recipes under `skills/`; they are more current than this file.

### Blocker: CUDA 13 torch vs CUDA 12.4 driver

`bash skills/rfd4-proteina-setup/scripts/preflight.sh` →
`FAIL: CUDA is unavailable in the selected environment.`

Confirmed on a real GPU node (job 2, `logs/probe-2.log`), so it is **not** a
head-node artifact:

```
RFD4 venv:  torch 2.14.1+cu132  cuda_available False
            "The NVIDIA driver on your system is too old (found version 12040)"
env-rfd3:   torch 2.7.1+cu126    cuda_available True   NVIDIA L4
```

Cluster driver is 550.127.08 (CUDA 12.4). A CUDA 13.x torch needs driver ≥ 580.

### Fix path (not yet applied)

The project supports CUDA 12 officially — `pyproject.toml` declares an
`all-cuda12` extra (line 111) alongside `all-cuda13`, differing mainly in
`cuequivariance-ops-torch-cu12` vs `-cu13`, and pins torch only as
`torch>=2.11,<3`.

`torch-2.14.1+cu126-cp312-cp312-manylinux_2_28_x86_64.whl` exists on
`download.pytorch.org/whl/cu126` — the **same torch version**, CUDA 12.6 build.
So the minimal change is to swap the CUDA variant and the `cuequivariance` ops
package, holding every other version constant.

Caveats to check when applying it:
- `setup.sh` exposes no CUDA selector (`--native` only picks venv vs cluster
  installer), so the extra must be selected through `uv` directly.
- Other compiled deps may also be cu13-built; `cuda-python>=13.3.1,<14` and
  `nvidia-cutlass-dsl` are cuda13-extra-only, but verify nothing else drags in a
  CUDA 13 runtime.
- **Rebuild into a separate venv**, not over `.venv` — the repo's `AGENTS.md`
  says not to reinstall a working environment, and keeping the original lets us
  fall back and keeps CPU-only paths usable.

### Interface facts worth keeping

- **Hydra `key=value`, never `--key value`.** `rfd4-generate inputs=<f>
  num_replicates=<N> out_dir=<d>`.
- **Enzyme design is a first-class task** with benchmark inputs at
  `data/benchmarks/mcsa/{mcsa_37, mcsa_37_fixed_ligand, mcsa_37_no_ligand}`
  (37 M-CSA active sites, `length` 180).
- **The mcsa input vocabulary matches the hackathon's `Kemp.json` exactly** —
  `ligand`, `unindex`, per-residue `select_fixed_atoms`. The organisers' Kemp
  motif spec should therefore port to RFD4 with little translation. This is the
  main reason RFD4 is worth unblocking.
- `recycle_strategy=medium/midpeak` is the documented schedule for MCSA /
  small-molecule tasks (120/150/180 passes). `AGENTS.md` is explicit that the
  `thinking_*` aliases are front-loaded, **not** task-aware defaults.
- Diversity/temperature knobs: `generation.model.bb_ca.simulation_step_params.tsr_k`
  (default 16) and `...local_latents.simulation_step_params.tsr_k` (default 4).
  Higher k = lower temperature.
- Partial diffusion (needed for plan §6.9 maturation):
  `generation.args.partial_trajectory_fraction_steps=<f>`, or
  `partial_t_backbone` / `partial_t_latent` in (0,1); near 1 stays close to input.
  `partial_t` is a corruption level, **not** Ångströms.
- **Conditioning rules that matter for the theozyme (§5.2, §6.2):**
  - Fixing any *sidechain* coordinate forces sequence conditioning — the atom37
    slot-occupancy pattern fingerprints the residue. Backbone-only (`N, CA, C, O`)
    leaks nothing and is how you fix a fold and redesign sequence.
  - Unindexed/floating motifs use `C_INTCRD` + `C_CRDGRP` (group `0` reserved for
    absolute `C_CRD`, `-1` for no known relative geometry). `C_CRD` and `C_INTCRD`
    are mutually exclusive per atom.
  - Ligands/atomised tokens are always sequence-conditioned and never generated.
  - **Aromatic bonds must be Kekulised** to `AROMATIC_SINGLE`/`AROMATIC_DOUBLE`;
    a bond left `AROMATIC` is silently cast to `SINGLE` and the ring's bond-order
    information is lost. Formal charges are read as written — "unknown" is not
    neutral. Directly relevant to `6NT.sdf` and `ts.sdf` (plan §5.1).
### The theozyme conditioning recipe — WORKING (got it wrong twice first)

**This is the recipe that works.** Verified: a de novo 185-residue scaffold
reproducing the HG3.17 site at base 2.57 Å / donor 3.23 Å / wall 5.43 Å
(crystal: 2.53 / 3.10 / 5.53).

```yaml
kemp_theozyme_unindexed:
  input: ./kemp_motif.cif          # .cif only; .pdb is rejected
  contig: "E1:150-250, A44, A50, A127, /0, L1305"
  conditions:
    C_SEQ:                          # fix catalytic identities
      true:
        - select: "sel('A/*/[44, 50, 127]')"
    C_CRD:                          # THE geometric target -- absolute
      true:
        - select: "sel('A/*/[44, 50, 127]') or sel('L')"
    C_IDX:                          # unindexed: model picks sequence positions
      false:
        - select: "sel('A/*/44')"
          group: wall
        - select: "sel('A/*/50')"
          group: donor
        - select: "sel('A/*/127')"
          group: base
```

**Use `C_CRD: true` (absolute), not `C_CRDGRP`.** The docs say it plainly and we
ignored it: `C_CRD` "works for a single rigid motif, but it can't express
'these atoms are rigid relative to each other, but I don't know or care where
that rigid body sits in space'". A theozyme **is** a single rigid motif, so
absolute conditioning is correct. `C_CRDGRP` exists for *several
independently placed* motifs. The global frame being arbitrary does not matter,
because the scaffold is generated around the motif.

What we learned the hard way:

| Attempt | Result |
| --- | --- |
| `C_SEQ` + `C_CRDGRP: 1`, no `C_CRD` | **validates, runs, silently wrong.** 0/10 designs kept the site; base→N3 8.1–23.0 Å against a 2.5 Å target. Group recorded, no geometric target existed. |
| `C_CRD: true` + `C_CRDGRP: 1` | rejected at prevalidation: *"Coordinate group 0 must match absolute coordinate conditioning; modify C_CRD first."* `C_CRD: true` reserves group **0**, conflicting with group 1. |
| `C_CRD: false` + `C_CRDGRP: 1` | validates, **still wrong** (base 13.1 Å). `C_CRD: false` really does mean *no* coordinate conditioning. |
| **`C_CRD: true`, no `C_CRDGRP`** | **works** — base 2.57, donor 3.23, wall 5.43 Å |

For reference, the canonical `motif_unindexed` benchmark input carries **no**
coordinate-condition column at all: just real coordinates,
`condition_sequence` set, `condition_index = False`, and **one shared**
`condition_index-group` for the whole motif. Its motif is a contiguous segment,
so a shared index group (which preserves relative residue offsets) suits it; a
discontiguous theozyme wants separate groups so each catalytic residue's
sequence position is independent.

Note our `C_IDX: false` with separate group labels still came back as
`condition_index_1_residue = True` in the output (while
`is_unindexed_insertion = True`). The geometry is right regardless, but the
unindexing may not be doing what we intend. Worth revisiting if fold diversity
looks poor.

**Lesson for the pipeline:** `prevalidate_inputs=true` catches *schema*
inconsistency but not a *semantically* empty motif. Only measuring the output
geometry catches that. `geometry.py` is therefore load-bearing as a generation
gate, not just a ranking metric — every batch must be motif-checked before the
designs are used.

Other interface facts confirmed while debugging:

- `input` must be `.cif`/`.cif.gz`/`.bcif`; a `.pdb` is rejected.
- The ligand needs its **own chain** in the input. With the ligand sharing the
  protein's chain id, a contig with a `/0` before it fails:
  *"Contig produces two separate output chains both named 'A'."* The MCSA
  benchmark inputs do the same (ligands on their own chains).
- Output provenance columns that actually matter for locating a moved motif:
  `condition_sequence_1_residue` (real resname vs `<M>`),
  `annotation_src-coordinate-conditioned-region_1_atom`,
  `condition_CoordinateGroup_1_atom`, `is_unindexed_insertion`.
- Our `C_IDX: false` came back as `condition_index_1_residue = True` in the
  output while `is_unindexed_insertion = True`. The two are not the same flag;
  unindexing did happen. Not yet fully understood — low priority.
- **All `*_refolding` metric tags need an AF3 apptainer/docker container.** No
  container runtime exists here (see CLUSTER.md), so `ame_refolding` and
  `motif_scaffold_refolding` are **unavailable**. Motif RMSD must be computed
  ourselves in `geometry.py`. Plain generation needs no container.
- Outputs are annotated CIFs under `inference_samples/` with provenance columns
  (`SA_SRC`, `SA_CRDSRC`, `SA_IDXGRPSRC`, `SA_UII`) tying output atoms to the
  input spec. `save_chain_of_interest=true` is required for later clustering.
- `CCD_MIRROR_PATH` and `PDB_MIRROR_PATH` are unset; atomworks warns. Set them if
  a CCD lookup is needed (relevant to fetching `6NT`).

---

## RFD3 — verified by the organisers, our fallback for stage 1

- Run via `$ENV_ROOT/env-rfd3/bin/rfd3`, Hydra-style `key=value`.
- Known-good Kemp invocation is in `/mnt/work/kemp_example/src/kemp-e2e.sbatch`.
  Engine args (`ckpt_path`, `out_dir`, `diffusion_batch_size`, `n_batches`,
  `seed`, `prevalidate_inputs`, `skip_existing`, `cleanup_guideposts`,
  `output_full_json`, `dump_prediction_metadata_json`) go **on the command line**;
  the input JSON's keys are all `DesignInputSpecification` fields.
- Input JSON format (`inputs/Kemp.json`): per-design-key `input` (motif PDB),
  `ligand` (CCD code), `contig` (`"260-260"` fixed, `"150-260"` sampled),
  `unindex` (comma-separated chain+resnum), `select_fixed_atoms` (per-residue atom
  name lists, or `ALL` for the ligand).
- Select design keys with `'json_keys_subset=[Kemp_HG3_17_base_donor]'`.
- **Open:** whether RFD3 exposes the RASA conditioning plan §5.2/§6.2 needs for
  the burial sweep. `logs/rfd3-help.txt` in the example workspace is the place to
  look.

## LigandMPNN — provided, known-good invocation

- `$ENV_ROOT/env-ligandmpnn/bin/python $LIGANDMPNN_ROOT/run.py`, conventional
  `--flag value`.
- Verified flags: `--model_type ligand_mpnn`, `--checkpoint_ligand_mpnn`,
  `--pdb_path`, `--out_folder`, `--fixed_residues` (space-separated
  `ChainResnum` list), `--ligand_mpnn_use_atom_context 1`,
  `--ligand_mpnn_use_side_chain_context 1`, `--batch_size`,
  `--number_of_batches`, `--pack_side_chains`, `--seed`, `--verbose`.
- Checkpoint `ligandmpnn_v_32_010_25.pt` = v32, 0.10 Å noise. Check whether a
  soluble-weighted variant is present (plan §6.3).
- Note for plan §6.3: `--pack_side_chains 1` is needed if stage 3/4 want packed
  rotamers; the example uses `0`.

## MMseqs2 — provided, GPU-accelerated

- `$MMSEQS` v18-8cc5c. Target DB `$MMSEQS_DB` (= `targetDB_gpu`), with
  `$MMSEQS_DB.idx` used as the search target.
- Verified pipeline: `createdb --dbtype 1 --shuffle 0` → `search --db-load-mode 2
  --gpu 1 --max-seqs 300 -e 0.001` → `result2msa --msa-format-mode 5` →
  `unpackdb --unpack-name-mode 0 --unpack-suffix .a3m`.
- Needs ~28 GB RAM → **must run on `gpu-dy-l4-1`** (the g6.2xlarge).
- Serves double duty: MSAs for Boltz-2 (stage 5) *and* the identity check in the
  compliance gate (stage 6). Confirm what `targetDB_gpu` actually contains —
  stage 6 needs UniProt + PDB seqres coverage to be meaningful.

## Boltz-2 — provided, known-good invocation

- `$ENV_ROOT/env-boltz2/bin/boltz predict <yaml> --model boltz2 --cache <dir>
  --out_dir <dir> --accelerator gpu --devices 1 --no_kernels --sampling_steps 200
  --recycling_steps 3 --diffusion_samples 1 --num_workers 0
  --preprocessing-threads 1 --output_format mmcif --seed 1`
- The cache dir must contain symlinks to `boltz2_conf.ckpt`, `boltz2_aff.ckpt`,
  `mols.tar` and `mols/` from `$BOLTZ_ASSETS`.
- `--no_kernels` is used by the verified example; keep it.
- Open (plan §9.4, §9.5): whether an arbitrary SDF is accepted or a CCD code is
  required, and the real wall-clock of one prediction on an L4.

## tmol — not attempted (next, timeboxed 2 h)

On PyPI (HTTP 200 from both node classes). Plan §6.4 ranking metric
`tmol.ddE_ts_sub` depends on it. Fallback if ligand parameterisation fails:
geometry-only gate, drop ΔΔE (plan §7).

## PLACER · BioEmu-2 · Foldseek · TM-align — not attempted

PLACER and BioEmu-2 are v1/v2 scope. Foldseek + TM-align/US-align are needed by
the stage 6 compliance gate and are **not installed** — they must be built before
Milestone 3.

---

## §9 open questions — current answers

| # | Question | Answer |
| --- | --- | --- |
| 1 | RFD4/RFD3 RASA conditioning: per-atom? what range? | open |
| 2 | Does PLACER batch on the GPU? | open (v1) |
| 3 | Can tmol parameterise the TS template? | open |
| 4 | Boltz-2: arbitrary SDF or CCD code only? | open |
| 5 | Boltz-2 wall-clock on this GPU? | open |
| 6 | Is `6NT` inside the Boltz-2/AF3 training window? | open — if yes the substrate control arm in stage 5 is load-bearing |
| — | Compute-node internet? | **yes**, pypi + RCSB both 200 |
| — | Container runtime for AF3 refolding metrics? | **none** — RFD4 refolding tags unusable |
