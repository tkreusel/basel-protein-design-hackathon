"""Design gate, refold gate and extended set reproduce the tables of the campaign."""
import pandas as pd

from kemp import boltz_gate, extended_set, filter_pre
from conftest import read, FIXTURES


def test_design_gate_reproduces_the_campaign(run_name):
    new = filter_pre.apply_gate(read(run_name, "design_metrics.csv"), min_keep=30)
    old = read(run_name, "expected_gate.csv")
    for col in old.columns:
        assert (new[col].values == old[col].values).all(), col


def test_fallback_fills_up_to_min_keep_with_tier_3_when_nothing_passes_geometry():
    new = filter_pre.apply_gate(read("hgnt__unsteered_sheet__b51", "design_metrics.csv"), min_keep=30)
    assert (new.gate_tier == 3).sum() > 0 and new.gate_pass.sum() == 30


def test_trp_motif_tolerates_one_polar_atom_near_the_base():
    d = read("hgwt__unsteered__b01", "design_metrics.csv").head(2).copy()
    d["base_sc_polar_4p5"] = [1, 2]
    geometry, sequence = filter_pre.criteria(d)
    assert list(sequence["no_sc_polar"]) == [True, False]
    d["motif"] = "hgnt"
    assert list(filter_pre.criteria(d)[1]["no_sc_polar"]) == [False, False]


def test_refold_gate_reproduces_the_campaign(run_name):
    old = read(run_name, "boltz1_eval.csv")
    derived = [c for c in old if c.startswith("ok_") or c in ("geom_ok_refold", "boltz_gate_strict", "n_crit_ok", "boltz_tier", "boltz_gate")]
    new = boltz_gate.apply_gate(old.drop(columns=derived))
    for col in derived:
        assert (new[col].values == old[col].values).all(), col


def test_refold_fallback_requires_folding_except_for_sheet_runs():
    old = read("hgnt__steer_tmol__b01", "boltz1_eval.csv")
    new = boltz_gate.apply_gate(old.drop(columns=[c for c in old if c.startswith("ok_")]))
    fallback = new[new.boltz_tier == 2]
    assert len(fallback) and (fallback.ok_plddt & fallback.ok_ca).all()


def test_plddt_scale_detection():
    assert boltz_gate.plddt_scale(pd.Series([0.4, 0.9])) == 1.0
    assert boltz_gate.plddt_scale(pd.Series([40.0, 90.0])) == 100.0


def test_extended_set_reproduces_the_campaign():
    expected = read("expected_extended.csv")
    for name in ("hgnt__steer_tmol__b01", "hgwt__unsteered__b01"):
        motif, arm, batch = name.split("__")
        strict = set((FIXTURES / name / "survivors.txt").read_text().split())
        chosen = set(extended_set.select_extended(read(name, "boltz1_eval.csv"), strict).seq_id)
        assert chosen == set(expected[(expected.motif == motif) & (expected.arm == arm) & (expected.batch == batch)].seq_id)
