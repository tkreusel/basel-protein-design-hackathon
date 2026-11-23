"""Write Boltz-2 input YAMLs for both design arms (plan stage 5).

Arm A = RFD4's co-designed sequences.
Arm B = LigandMPNN re-designs of the same backbones.

Both are refolded with 6NT bound so the comparison is like-for-like.

Two choices worth stating:

  msa: empty  -- single-sequence mode (supported, boltz schema.py:1127). De
      novo designs have no database homologues; the organisers' own example
      notes a designed sequence "may have no database homologues and use a
      query-only alignment". Single-sequence mode avoids a 183 GB MMseqs2
      search that would return the query row and nothing else, and avoids
      needing the 28 GB-RAM node. If we later want MSAs for the natural-protein
      control, that is a separate run.

  ligand by CCD code -- Boltz-2 accepts `ccd: 6NT` (verified in the
      organisers' working example), so no SDF is needed here. That answers
      plan 9 question 4 for this code path.

Usage: python -I prepare_boltz_both.py <stage dir> <out dir> [max_per_arm]
"""

from __future__ import annotations

import sys
from pathlib import Path

import yaml

LIGAND_CCD = "6NT"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    if not path.exists():
        return []
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


def main(stage_dir: Path, out_dir: Path, max_per_arm: int = 0) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    arms = {
        "armA": stage_dir / "armA_rfd4_sequences.fasta",
        "armB": stage_dir / "armB_ligandmpnn_sequences.fasta",
    }

    index = []
    for arm, fa in arms.items():
        entries = read_fasta(fa)
        if max_per_arm:
            entries = entries[:max_per_arm]
        print(f"{arm}: {len(entries)} sequence(s) from {fa.name}")
        for hdr, seq in entries:
            tag = hdr.split()[0]
            name = f"{arm}_{tag}"
            spec = {
                "version": 1,
                "sequences": [
                    {"protein": {"id": "A", "sequence": seq, "msa": "empty"}},
                    {"ligand": {"id": "B", "ccd": LIGAND_CCD}},
                ],
            }
            p = out_dir / f"{name}.yaml"
            p.write_text(yaml.safe_dump(spec, sort_keys=False))
            index.append({"name": name, "arm": arm, "design": tag,
                          "length": len(seq), "yaml": str(p), "header": hdr})

    (out_dir / "index.yaml").write_text(yaml.safe_dump(index, sort_keys=False))
    print(f"\nwrote {len(index)} Boltz-2 input(s) to {out_dir}")
    if not index:
        raise SystemExit("no sequences to refold -- did stage 1/2 produce output?")


if __name__ == "__main__":
    main(
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        int(sys.argv[3]) if len(sys.argv) > 3 else 0,
    )
