"""Motif location and the LigandMPNN residue sets against the sequence tables of the campaign."""
import pytest

pytest.importorskip("biotite")
from kemp import config, mpnn, structure
from conftest import FIXTURES, read

CASES = [("hgnt__steer_tmol__b01", True), ("hgwt__unsteered__b01", False)]


@pytest.mark.parametrize("run,steered", CASES)
def test_motif_omit_and_shell_sets_match_the_campaign(run, steered):
    seqs = read(run, "sequences_all.csv")
    sel = read(run, "sel.csv") if (FIXTURES / run / "sel.csv").exists() else None
    checked = 0
    for cif in sorted((FIXTURES / run / "designs").glob("*.cif")):
        design = cif.stem
        rows = seqs[seqs.design == design]
        if rows.empty:
            continue
        motif = run.split("__")[0]
        protein, ligand = mpnn.split_design(cif)
        ctx = mpnn.design_context(protein, ligand, motif, steered)
        native = rows[rows.source == "rfd4"].iloc[0]
        assert " ".join(f"A{i}" for i in ctx["motif_ids"]) == native.fixed
        assert len(ctx["omit_ids"]) == native.n_omit
        assert structure.sequence(protein) == native.sequence
        if steered:
            assert len(ctx["shell_ids"]) == rows[rows.variant == "shell"].n_shell_fixed.iloc[0]
        else:
            assert ctx["shell_ids"] is None
        checked += 1
    assert checked >= 1


def test_motif_residues_are_not_in_the_omit_set():
    cif = next((FIXTURES / "hgnt__steer_tmol__b01" / "designs").glob("*.cif"))
    protein, ligand = mpnn.split_design(cif)
    ctx = mpnn.design_context(protein, ligand, "hgnt", True)
    assert not set(ctx["motif_ids"]) & set(ctx["omit_ids"])
    assert set(ctx["motif_ids"]) <= set(ctx["shell_ids"])
