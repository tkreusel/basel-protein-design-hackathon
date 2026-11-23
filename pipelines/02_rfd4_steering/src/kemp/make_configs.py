"""Write the RFD4 experiment configs for the steered arms: rfd4/configs/inference_experiment/kc_<motif>_<arm>.yaml.

usage: make_configs.py [--out DIR]       (default: <pipeline>/rfd4/configs/inference_experiment)

The reward sets, search schedule and arms come from config/arms.yaml. The unsteered arm has no experiment file.
"""
import argparse
from pathlib import Path

from . import config


def reward_terms(arms_cfg, channel_band):
    """Ordered {name: {weight, clip, normalize_by}}; the TMol term always comes last."""
    terms = {k: dict(v) for k, v in arms_cfg["rewards"].items()}
    if channel_band:
        terms.update({k: dict(v) for k, v in arms_cfg["channel_band"]["override"].items()})
        terms.update({k: dict(v) for k, v in arms_cfg["channel_band"]["add"].items()})
    terms.update({k: dict(v) for k, v in arms_cfg["tmol_reward"].items()})
    return terms


def metric_list(arms_cfg, channel_band):
    return arms_cfg["metrics"] + (arms_cfg["channel_band"]["metrics"] if channel_band else [])


def _fmt_value(v):
    if isinstance(v, dict):
        return "{" + ", ".join(f"{k}: {_fmt_value(x)}" for k, x in v.items()) + "}"
    return str(v)


def _fmt_term(name, t):
    return f"{name}: {_fmt_value(t)}"


def experiment_yaml(motif, arm, arms_cfg=None):
    arms_cfg = arms_cfg or config.load("arms")
    arm_def = arms_cfg["arms"][arm]
    assert arm_def["steered"], f"{arm} has no search"
    band = arm_def.get("channel_band", False)
    rewards = ", ".join(_fmt_term(k, v) for k, v in reward_terms(arms_cfg, band).items())
    metrics = ", ".join(metric_list(arms_cfg, band))
    search = arms_cfg["search"]
    lines = []
    for step in search["checkpoints"]:
        clean = "use_clean_prediction: true, " if step <= search["clean_prediction_until"] else ""
        lines.append(f'    "{step}": {{{clean}metrics: [{metrics}], rewards: {{{rewards}}}}}')
    spec = config.motif_spec(motif)
    return (f"# @package _global_\n# campaign experiment / {motif} / {arm}\n"
            f"defaults:\n  - inference_on_user_inputs\n  - _self_\n"
            f"inputs: ${{oc.env:PROJECT_ROOT}}/../inputs/{spec}\n"
            f"run_name_: {motif}_{arm}\n"
            f"search:\n  n_branch: {search['n_branch']}\n  beam_width: {search['beam_width']}\n  checkpoints:\n" + "\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(config.PIPELINE / "rfd4" / "configs" / "inference_experiment"))
    args = ap.parse_args()
    arms_cfg = config.load("arms")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for motif in config.load("motifs")["motifs"]:
        for arm, d in arms_cfg["arms"].items():
            if d["steered"] and not d.get("sheet"):      # the *_sheet arms reuse the experiment of their base arm
                (out / f"kc_{motif}_{arm}.yaml").write_text(experiment_yaml(motif, arm, arms_cfg))
                print("wrote", out / f"kc_{motif}_{arm}.yaml")


if __name__ == "__main__":
    main()
