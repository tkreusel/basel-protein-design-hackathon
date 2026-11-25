"""Arms, seeds, RFD4 overrides and the generated experiment configs against what the campaign ran."""
import pytest

from kemp import config, gen, make_configs
from conftest import FIXTURES, read


def test_seeds_match_the_campaign(run_name):
    d = read(run_name, "design_metrics.csv")
    for (motif, arm, batch, chunk), seed in d.groupby(["motif", "arm", "batch", "chunk"]).seed.first().items():
        assert gen.seed_for(int(batch), arm, motif, chunk_index=int(chunk) // 25) == seed


def test_seed_is_distinct_per_motif_arm_and_batch():
    arms = config.load("arms")["arms"]
    seeds = {gen.seed_for(b, arm, motif) for b in range(1, 11) for arm in arms for motif in ("hgnt", "hgwt")}
    assert len(seeds) == 10 * len(arms) * 2


def test_replicates_and_members():
    arms = config.load("arms")["arms"]
    assert gen.replicates_for("unsteered") * arms["unsteered"]["members"] == 300
    assert gen.replicates_for("steer_tmol") * arms["steer_tmol"]["members"] == 300


def test_unsteered_drops_the_metrics_group_and_steered_selects_the_experiment():
    unsteered = gen.hydra_overrides("hgwt", "unsteered", 3, "/o", 25, 3005, 4)
    assert "~metrics" in unsteered and not any(o.startswith("inference_experiment") for o in unsteered)
    steered = gen.hydra_overrides("hgwt", "steer_tmol_ch18_sheet", 3, "/o", 10, 3055, 4)
    assert "inference_experiment=kc_hgwt_steer_tmol_ch18" in steered
    assert any(o.endswith("kemp4_hg317-tips_with-trp_sheet.yaml") for o in steered)


@pytest.mark.parametrize("motif,arm", [("hgnt", "steer_tmol_ch18"), ("hgwt", "steer_tmol")])
def test_experiment_yaml_is_identical_to_the_one_that_ran(motif, arm):
    expected = (FIXTURES / f"expected_kc_{motif}_{arm}.yaml").read_text().splitlines()
    produced = make_configs.experiment_yaml(motif, arm).splitlines()
    assert produced[:1] == expected[:1] and produced[2:] == expected[2:]     # line 2 is a comment naming the campaign


def test_channel_arm_changes_exactly_the_documented_terms():
    arms = config.load("arms")
    base, band = make_configs.reward_terms(arms, False), make_configs.reward_terms(arms, True)
    assert set(band) - set(base) == {"kemp_channel18.channel_band_error"}
    assert band["kemp_pocket3.ligand_sasa_frac"]["weight"] == -0.25 and base["kemp_pocket3.ligand_sasa_frac"]["weight"] == -0.5
    assert list(band)[-1] == "tmol_lig_ddg.interface_ddg"
