"""Stage-5 report, shortlist flags, composite score, panel rank and final geometry against the campaign tables."""
import numpy as np
import pandas as pd
import pytest

from kemp import final_scoring, flags, report, select_models
from conftest import FIXTURES, read


def test_stage5_report_reproduces_the_campaign():
    run = "hgnt__steer_tmol__b01"
    seqs = read(run, "sequences.csv").set_index("seq_id")
    new = select_models.stage5_report(read(run, "ranking.csv"), seqs).sort_values("seq_id").reset_index(drop=True)
    old = read(run, "stage5_report.csv").sort_values("seq_id").reset_index(drop=True)
    assert list(new.columns) == list(old.columns)
    num = old.select_dtypes("number").columns
    assert np.allclose(new[num].values, old[num].values, equal_nan=True, atol=1e-9)


def test_shortlist_flags_reproduce_the_campaign():
    ev = read("evaluated_inputs.csv")
    expected = read("expected_flags.csv")
    new = flags.add_flags(ev)
    for col in expected.columns[1:]:
        if col.startswith("z_"):
            assert np.allclose(new[col].values, expected[col].values, equal_nan=True, atol=1e-9), col
        else:
            assert (new[col].values == expected[col].values).all(), col


def test_composite_reproduces_final_ranking():
    ev = read("evaluated_inputs.csv")
    ev = ev[~ev.arm.str.contains("sheet")]
    new = final_scoring.composite(final_scoring.add_energy_terms(ev), ipsae_min=0.8, min_metrics=5)
    expected = read("expected_final_scoring.csv")
    assert (new.seq_id.values == expected.seq_id.values).all()
    scored = expected.composite_z.notna().values
    assert np.allclose(new.composite_z.values[scored], expected.composite_z.values[scored], atol=1e-9)
    assert (new.final_rank.values == expected.final_rank.values).all()


def test_composite_orientation_lower_raw_is_better():
    df = pd.DataFrame({"seq_id": list("abc"), "ipsae.tsa": [0.9] * 3})
    for i, c in enumerate(final_scoring.SCORE_COLUMNS):
        df[c] = [1.0, 2.0, 3.0]
    out = final_scoring.composite(df, 0.8, 5)
    assert list(out.seq_id) == ["a", "b", "c"]


def test_panel_rank_reproduces_the_campaign():
    ev = read("evaluated_inputs.csv")
    ranked, _ = report.panel_rank(ev)
    expected = read("expected_panel_rank.csv")
    assert (ranked.seq_id.values == expected.seq_id.values).all()
    assert np.allclose(ranked.panel_rank.values, expected.panel_rank.values, atol=1e-9)


def test_final_geometry_reproduces_the_campaign():
    pytest.importorskip("biotite")
    from kemp import final_selection, structure
    run = FIXTURES / "hgnt__unsteered__b01"
    expected = read("hgnt__unsteered__b01", "expected_geometry.csv")
    seqs = read("hgnt__unsteered__b01", "sequences_all.csv").set_index("seq_id")
    theozyme = structure.load_reference()
    for x in expected.itertuples():
        fixed = [int(t[1:]) for t in str(seqs.loc[x.seq_id].fixed).split()]
        g = final_selection.geometry_metrics(structure.load_cif(run / "designs" / f"{x.design}.cif"),
                                             [structure.load_cif(run / "best_boltz" / f"{x.seq_id}__6NT_model_0.cif")], fixed, x.motif, theozyme)
        for key in ("design_tips_rmsd_theozyme", "boltz_best_tips_rmsd_theozyme", "asp_O_to_N3_A", "gln_NE2_to_N1_A",
                    "n_polar_atoms_3p5_of_N1", "base_CB_within_12A_boltz", "base_CB_within_12A_design"):
            assert g[key] == pytest.approx(getattr(x, key), abs=1e-3), key
