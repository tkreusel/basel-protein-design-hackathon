"""Discrimination z-scores and shortlist flags for evaluated sequences (report only; thresholds in config/gates.yaml: shortlist).

z_tsa_vs_sub / z_tsa_vs_prod: difference of the mean Boltz ipTM between the TSA samples and the substrate (product) samples divided by its
standard error over the samples of both ligands. They say whether Boltz places the TSA more confidently than the other ligand for that
sequence; they are not binding energies.
"""
import numpy as np
import pandas as pd

from . import config


def _n(e, lig):
    col = f"{lig}.n"
    return e[col] if col in e else 10


def add_flags(ev, cfg=None):
    cfg = cfg or config.load("gates")["shortlist"]
    ev = ev.copy()
    for lig, tag in (("H5J", "sub"), ("6VP", "prod")):
        se = np.sqrt(ev["6NT.iptm.sd"] ** 2 / _n(ev, "6NT") + ev[f"{lig}.iptm.sd"] ** 2 / _n(ev, lig) + 1e-9)
        ev[f"z_tsa_vs_{tag}"] = (ev["6NT.iptm.mean"] - ev[f"{lig}.iptm.mean"]) / se
    ev["flag_strict_boltz"] = ev.boltz_tier == 1
    ev["flag_placer_signal"] = (ev["placer.apo.frac_pass_core"] >= cfg["placer_apo_frac_min"]) | (ev["placer.holo.frac_pass_core"] >= cfg["placer_holo_frac_min"])
    ev["flag_tsa_over_sub"] = ev.z_tsa_vs_sub > cfg["tsa_over_substrate_z_min"]
    ev["flag_tsa_over_prod"] = ev.z_tsa_vs_prod > cfg["tsa_over_substrate_z_min"]
    ev["flag_shortlist"] = (ev.flag_strict_boltz & ev.flag_placer_signal & (ev["ipsae.tsa"] >= cfg["ipsae_tsa_min"]) & ev.flag_tsa_over_sub
                            & (ev["ros.tsa.dG_sep"] <= cfg["ros_dg_sep_max"]))
    ev["family"] = ev.arm.map(config.arm_family)
    return ev
