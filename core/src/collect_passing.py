"""Collect the designs that passed every gate into one obvious place.

Reads a stage manifest, copies each passing design's structure and writes a
ranked summary. Nothing here re-scores anything -- it only gathers what the
gate already decided, so the summary cannot disagree with the manifest.

Ranking is by ligand burial (rel_sasa, lower = better enclosed), because in the
stage1_13 batch catalytic geometry was essentially identical across all designs
(base 2.37-2.73 A) while burial was what actually separated them (0.007-0.54).

Usage: python -I collect_passing.py <stage dir> <out dir>
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kemp.record import read_manifest  # noqa: E402


def main(stage_dir: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    cif_dir = out_dir / "structures_cif"
    pdb_dir = out_dir / "structures_pdb"
    cif_dir.mkdir(exist_ok=True)
    pdb_dir.mkdir(exist_ok=True)

    recs = list(read_manifest(stage_dir / "manifest.jsonl"))
    passing = [
        r for r in recs
        if r.status == "ok" and r.gates.get("stage1_motif") == "pass"
    ]
    passing.sort(key=lambda r: r.scores.get("geom.rel_sasa", 9.9))

    rows = []
    fasta = []
    for rank, r in enumerate(passing, start=1):
        src_cif = Path(r.paths.get("rfd4_cif", ""))
        if src_cif.exists():
            dest = cif_dir / f"{rank:02d}_{r.design_id}.cif"
            shutil.copy2(src_cif, dest)
        src_pdb = Path(r.paths.get("backbone_pdb", ""))
        if src_pdb.exists():
            shutil.copy2(src_pdb, pdb_dir / f"{rank:02d}_{r.design_id}.pdb")

        cat = r.catalytic_residues
        wall = cat["wall"][0] if isinstance(cat.get("wall"), list) else cat.get("wall")
        s = r.scores
        rows.append(dict(
            rank=rank, design_id=r.design_id, length=r.length,
            base_pos=cat.get("base"), donor_pos=cat.get("donor"), wall_pos=wall,
            base_dist=s.get("geom.base_to_abstraction"),
            donor_dist=s.get("geom.donor_to_oxyanion"),
            wall_dist=s.get("geom.wall_stacking"),
            rel_sasa=s.get("geom.rel_sasa"),
            enclosure=s.get("geom.enclosure"),
            n_contacts=s.get("geom.n_contacts"),
            run=r.run_id,
        ))
        fasta.append(
            f">{rank:02d}_{r.design_id} len={r.length} "
            f"base={cat.get('base')} donor={cat.get('donor')} wall={wall} "
            f"rel_sasa={s.get('geom.rel_sasa'):.3f} "
            f"enclosure={s.get('geom.enclosure'):.3f}\n{r.sequence}"
        )

    (out_dir / "sequences.fasta").write_text("\n".join(fasta) + "\n" if fasta else "")

    keys = list(rows[0]) if rows else []
    with open(out_dir / "summary.csv", "w") as fh:
        fh.write(",".join(keys) + "\n")
        for row in rows:
            fh.write(",".join(str(row[k]) for k in keys) + "\n")

    hdr = (f"{'#':>2} {'design_id':13} {'len':>4} {'base':>5} {'donor':>5} "
           f"{'wall':>5} {'rel_sasa':>9} {'encl':>6} {'cont':>5}")
    lines = [
        f"# Designs passing all stage-1 gates, from {stage_dir}",
        f"# {len(passing)} of {len(recs)} designs. Ranked by ligand burial.",
        "# Reference (HG3.17 crystal 7K4Z): base 2.53  donor 3.10  wall 5.53",
        "#                                  rel_sasa 0.011  encl 0.992  cont 66",
        "",
        hdr,
        "-" * len(hdr),
    ]
    for row in rows:
        lines.append(
            f"{row['rank']:>2} {row['design_id']:13} {row['length']:>4} "
            f"{row['base_dist']:>5.2f} {row['donor_dist']:>5.2f} "
            f"{row['wall_dist']:>5.2f} {row['rel_sasa']:>9.3f} "
            f"{row['enclosure']:>6.3f} {row['n_contacts']:>5.0f}"
        )
    text = "\n".join(lines) + "\n"
    (out_dir / "README.txt").write_text(text)
    print(text)
    print(f"structures (cif): {cif_dir}")
    print(f"structures (pdb): {pdb_dir}")
    print(f"sequences:        {out_dir / 'sequences.fasta'}")
    print(f"table:            {out_dir / 'summary.csv'}")


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
