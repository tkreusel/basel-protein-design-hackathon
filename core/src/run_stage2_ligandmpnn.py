"""Stage 2 (arm B): re-design RFD4 backbones with LigandMPNN.

Arm A is RFD4's own co-designed sequence. Arm B is LigandMPNN's sequence on the
*same* backbone, so the comparison isolates the sequence-design step.

Plan 6.3 requirements enforced here:
  * 6NT supplied as ligand context (a real molecule, not the TS template).
  * Catalytic residues FIXED, and the output VERIFIED to confirm they were not
    redesigned -- "assert this, do not trust it".
  * Any design whose catalytic residues come back changed is recorded failed.

Orchestration shells out via subprocess; no tool imports another tool's
environment (plan 0.4).

Usage: python -I run_stage2_ligandmpnn.py <project root> <stage dir> [n_seq]
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp.record import DesignRecord, read_manifest  # noqa: E402

ENV_ROOT = os.environ.get("ENV_ROOT", "/mnt/ref/software/envs")
LIGANDMPNN_ROOT = os.environ.get("LIGANDMPNN_ROOT", "/mnt/ref/software/LigandMPNN")
CKPT = os.environ.get(
    "LIGANDMPNN_CHECKPOINT",
    "/mnt/ref/data_weights/ligandmpnn/ligandmpnn_v_32_010_25.pt",
)

ONE_TO_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
EXPECTED = {"base": "D", "donor": "Q", "wall": "W"}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def main(root: Path, stage_dir: Path, n_seq: int = 4) -> None:
    manifest = stage_dir / "manifest.jsonl"
    records = {r.design_id: r for r in read_manifest(manifest)}
    out_root = stage_dir / "ligandmpnn"
    out_root.mkdir(parents=True, exist_ok=True)

    arm_b_lines: list[str] = []
    n_done = n_skip = n_fail = 0

    for did, rec in records.items():
        # Consult the GATE, not just the crash flag.
        #
        # This previously read `if rec.status != "ok"`, which only skips
        # designs that errored. Designs that ran fine but FAILED the motif
        # gate (no active site, or ligand on the surface) were redesigned
        # anyway, while arm A correctly excluded them -- so the two arms held
        # different designs and any A-vs-B comparison was meaningless.
        # See runs/_superseded/..._job14_CANCELLED-ran-on-gate-failed-designs/
        if rec.status != "ok" or rec.gates.get("stage1_motif") != "pass":
            n_skip += 1
            continue
        pdb = rec.paths.get("backbone_pdb")
        if not pdb or not Path(pdb).exists():
            rec.fail("stage2", f"backbone pdb missing: {pdb}")
            n_fail += 1
            continue

        # Resumable (plan 0.3): skip designs that already have stage-2 output.
        dest = out_root / did
        done_marker = dest / "seqs"
        if done_marker.exists() and any(done_marker.glob("*.fa")):
            n_skip += 1
            print(f"  [skip] {did} already has stage-2 output")
            continue

        cat = rec.catalytic_residues
        wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat.get("wall")
        # LigandMPNN --fixed_residues takes space-separated "<chain><resnum>".
        # These are the RFD4 output positions, recovered from provenance.
        fixed = f"A{cat['base']} A{cat['donor']} A{wall}"

        cmd = [
            f"{ENV_ROOT}/env-ligandmpnn/bin/python",
            f"{LIGANDMPNN_ROOT}/run.py",
            "--model_type", "ligand_mpnn",
            "--checkpoint_ligand_mpnn", CKPT,
            "--pdb_path", str(pdb),
            "--out_folder", str(dest),
            "--fixed_residues", fixed,
            "--ligand_mpnn_use_atom_context", "1",
            "--ligand_mpnn_use_side_chain_context", "1",
            "--batch_size", "1",
            "--number_of_batches", str(n_seq),
            "--pack_side_chains", "0",
            "--seed", "1",
            "--verbose", "0",
        ]
        print(f"  [run]  {did} fixed='{fixed}'")
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            rec.fail("stage2", f"LigandMPNN rc={proc.returncode}: {proc.stderr[-400:]}")
            n_fail += 1
            print(f"  [FAIL] {did}: rc={proc.returncode}")
            rec.write(stage_dir / "designs" / did)
            continue

        fas = sorted((dest / "seqs").glob("*.fa")) if (dest / "seqs").exists() else []
        if not fas:
            rec.fail("stage2", "LigandMPNN produced no .fa output")
            n_fail += 1
            rec.write(stage_dir / "designs" / did)
            continue

        entries = read_fasta(fas[0])
        # entry 0 is the input/native row; designs follow.
        kept = 0
        for idx, (hdr, seq) in enumerate(entries[1:], start=1):
            # VERIFY the fixed residues survived (plan 6.3: assert, don't trust)
            bad = []
            for role, pos in (
                ("base", cat["base"]), ("donor", cat["donor"]), ("wall", wall)
            ):
                if pos - 1 >= len(seq) or seq[pos - 1] != EXPECTED[role]:
                    got = seq[pos - 1] if pos - 1 < len(seq) else "?"
                    bad.append(f"{role}@{pos}={got} want {EXPECTED[role]}")
            if bad:
                rec.errors.append(f"stage2 seq{idx}: catalytic changed: {'; '.join(bad)}")
                continue
            arm_b_lines.append(
                f">{did}_mpnn{idx} armB_ligandmpnn len={len(seq)} "
                f"base={cat['base']} donor={cat['donor']} wall={wall}\n{seq}"
            )
            kept += 1

        if kept == 0:
            rec.fail("stage2", "every LigandMPNN sequence altered a catalytic residue")
            n_fail += 1
        else:
            rec.add_scores("stage2", {"n_sequences_kept": kept})
            n_done += 1
        rec.write(stage_dir / "designs" / did)
        print(f"         kept {kept}/{max(len(entries) - 1, 0)} sequences")

    out_fa = stage_dir / "armB_ligandmpnn_sequences.fasta"
    out_fa.write_text("\n".join(arm_b_lines) + ("\n" if arm_b_lines else ""))
    print(f"\nstage2: ok={n_done} skipped={n_skip} failed={n_fail}")
    print(f"arm B sequences: {out_fa}  ({len(arm_b_lines)} sequences)")


if __name__ == "__main__":
    main(
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        int(sys.argv[3]) if len(sys.argv) > 3 else 4,
    )
