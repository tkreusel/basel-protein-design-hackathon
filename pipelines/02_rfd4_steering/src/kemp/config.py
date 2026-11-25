"""Paths and configuration.

Everything site-specific comes from environment variables (see config/site.env.example), so the
scripts contain no cluster paths. The pipeline folder is found relative to this file.
"""
import os
from functools import lru_cache
from pathlib import Path

import yaml

PIPELINE = Path(os.environ.get("KEMP_PIPELINE", Path(__file__).resolve().parents[2]))
INPUTS = PIPELINE / "inputs"
CONFIG = PIPELINE / "config"


def env_path(name, default=None):
    value = os.environ.get(name)
    if value:
        return Path(value).expanduser()
    if default is None:
        raise RuntimeError(f"{name} is not set; source config/site.env (see config/site.env.example)")
    return Path(default)


def scratch():
    return env_path("KEMP_SCRATCH")


def campaign_name():
    return os.environ.get("CAMP", "20261122_kemp_c1")


def runs_root(camp=None):
    return scratch() / "runs" / (camp or campaign_name())


def run_dir(motif, arm, batch, camp=None):
    """runs/<camp>/<motif>/<arm>/bNN"""
    return runs_root(camp) / motif / arm / f"b{int(batch):02d}"


def results_dir():
    return env_path("KEMP_RESULTS")


@lru_cache(maxsize=None)
def load(name):
    """Load config/<name>.yaml."""
    with open(CONFIG / f"{name}.yaml") as fh:
        return yaml.safe_load(fh)


def motif_spec(motif, sheet=False):
    """File name of the RFD4 input spec for a motif (the *_sheet variant adds C_SEC=SHEET)."""
    name = load("motifs")["motifs"][motif]["spec"]
    return name[:-len(".yaml")] + "_sheet.yaml" if sheet else name


def arm_family(arm):
    """steer_tmol_sheet -> steer_tmol"""
    return arm[:-len("_sheet")] if arm.endswith("_sheet") else arm


def is_steered(arm):
    return load("arms")["arms"][arm]["steered"]
