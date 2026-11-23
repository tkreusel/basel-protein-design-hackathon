"""Stage 2 (CPU): gate the generated designs before sequence design.

usage: filter_pre.py --metrics gen/metrics.csv --out gate.csv [--designs gen/designs] [--min-keep 30]

Two flag sets are computed for every design, so that arms can be compared under identical filters afterwards:
  pass_geom   ligand and base buried, pocket volume in range, no clashes, motif RMSD small (sequence-independent)
  pass_seq    base environment of the RFD4-native sequence: no polar side chains or cations near the base,
              nitro group H-bonded away from the base, TMol interface energy < 0 (LigandMPNN can repair these)
  pass_core   pass_geom and pass_seq
The decision (gate_core) is pass_core for steered arms and pass_geom for unsteered arms (config/gates.yaml).

Fallback: if fewer than `min_keep` designs of a run (motif x arm x batch) pass the decision, the best remaining designs are added.
gate_tier 1 = passes the arm's gate; 2 = fallback that passes pass_geom; 3 = fallback that does not.
Fallback ranking: one design per lineage first, then more sequence criteria passed, fewer polar side-chain atoms near the base,
lower TMol dG, pocket ratio closer to 0.5.

Also written: lineage (replicate id), lineage_rank, channel_ok, channel18_ok, and <out>_strictness.csv
(per criterion: pass rate and the number of designs that fail only that criterion).
"""
import argparse

import numpy as np
import pandas as pd

from . import config

GEOMETRY_FLAGS = ["lig_sasa", "base_buried", "pocket", "clash", "motif"]
SEQUENCE_FLAGS = ["no_sc_polar", "no_cation", "nitro", "tmol"]


def criteria(d, gates=None, motifs=None):
    """Boolean Series per criterion (NaN counts as failing)."""
    gates = gates or config.load("gates")["design_gate"]
    motifs = motifs or config.load("motifs")["motifs"]
    g, s = gates["geometry"], gates["sequence"]
    polar_allowed = d.motif.map(lambda m: motifs[m]["extra_polar_near_base"])
    geometry = {
        "lig_sasa": d.ligand_sasa_frac <= g["ligand_sasa_frac_max"],
        "base_buried": d.base_sasa_rel <= g["base_sasa_rel_max"],
        "pocket": d.volume_ratio.between(*g["volume_ratio"]),
        "clash": (d.lig_clash_atoms <= g["lig_clash_atoms_max"]) & (d.pp_clash_frac <= g["pp_clash_frac_max"]),
        "motif": (d.motif_base_rmsd <= g["motif_rmsd_max"]) & (d.motif_partner_rmsd <= g["motif_rmsd_max"]),
    }
    sequence = {
        "no_sc_polar": d.base_sc_polar_4p5 <= polar_allowed,
        "no_cation": d.base_cation_6 <= s["base_cation_6_max"],
        "nitro": (d.nitro_hb_far >= s["nitro_hb_far_min"]) & (d.nitro_hb_near_base <= s["nitro_hb_near_base_max"]),
        "tmol": d.tmol_ddg < s["tmol_ddg_max"],
    }
    return geometry, sequence


def apply_gate(d, min_keep=None, gates=None, motifs=None):
    """Add flags, lineage information and the gated/fallback decision to the design table."""
    gates = gates or config.load("gates")["design_gate"]
    min_keep = gates["min_keep_designs"] if min_keep is None else min_keep
    d = d.copy()
    geometry, sequence = criteria(d, gates, motifs)
    for name, ok in {**geometry, **sequence}.items():
        d["ok_" + name] = ok.fillna(False)
    d["pass_geom"] = np.logical_and.reduce([d["ok_" + k] for k in geometry])
    d["pass_seq"] = np.logical_and.reduce([d["ok_" + k] for k in sequence])
    d["pass_core"] = d.pass_geom & d.pass_seq
    d["gate_core"] = np.where(d.arm.str.startswith("unsteered"), d.pass_geom, d.pass_core)

    d["lineage"] = d.design.str.replace(r"_beam\d+_\d+$", "", regex=True)
    ok_cols = [c for c in d if c.startswith("ok_")]
    d["n_ok"] = d[ok_cols].sum(axis=1)
    order = d.sort_values(["n_ok", "tmol_ddg"], ascending=[False, True])
    d["lineage_rank"] = order.groupby("lineage").cumcount().reindex(d.index) + 1
    d["n_seq_ok"] = d[["ok_" + k for k in sequence]].sum(axis=1)
    d["ratio_dev"] = (d.volume_ratio - 0.5).abs()

    d["gate_tier"] = np.where(d.gate_core, 1, 0)
    for _, idx in d.groupby(["motif", "arm", "batch"]).groups.items():
        run = d.loc[idx]
        need = min_keep - int((run.gate_tier == 1).sum())
        for tier, pool in ((2, run[(run.gate_tier == 0) & run.pass_geom]), (3, run[(run.gate_tier == 0) & ~run.pass_geom])):
            if need <= 0:
                break
            pool = pool.assign(not_first=pool.lineage_rank > 1).sort_values(
                ["not_first", "n_seq_ok", "base_sc_polar_4p5", "tmol_ddg", "ratio_dev"], ascending=[True, False, True, True, True])
            take = pool.head(need).index
            d.loc[take, "gate_tier"] = tier
            need -= len(take)
    d["gate_pass"] = d.gate_tier > 0
    d["channel_ok"] = d.bottleneck_radius >= gates["channel_open_min"]
    d["channel18_ok"] = d.bottleneck_radius.between(*gates["channel_band"])
    return d


def strictness(d):
    """Per motif/arm and criterion: pass rate and the number of designs failing only that criterion among all ok_* flags."""
    ok = [c for c in d if c.startswith("ok_")]
    rows = []
    for (motif, arm), x in d.groupby(["motif", "arm"]):
        for c in ok:
            only_this = int(((x[ok].sum(axis=1) == len(ok) - 1) & ~x[c]).sum())
            rows.append(dict(motif=motif, arm=arm, criterion=c[3:], n=len(x), pass_rate=round(x[c].mean(), 3), fails_only_this=only_this))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--designs", help="directory of design CIFs; adds a cif column")
    ap.add_argument("--min-keep", type=int)
    ap.add_argument("--select", help="also write the designs that pass the gate (input table of mpnn.py) to this CSV")
    a = ap.parse_args()

    d = apply_gate(pd.read_csv(a.metrics), a.min_keep)
    if a.designs:
        d["cif"] = a.designs.rstrip("/") + "/" + d.design + ".cif"
    d.to_csv(a.out, index=False)
    if a.select:
        chosen = d[d.gate_pass]
        chosen[["design", "cif", "motif", "arm", "gate_tier", "lineage"]].to_csv(a.select, index=False)
        print("selected for sequence design:", len(chosen), "designs,", chosen.lineage.nunique(), "lineages")

    by = d.groupby(["motif", "arm"])
    summary = by.agg(n=("design", "size"), lineages=("lineage", "nunique"), pass_geom=("pass_geom", "sum"), pass_seq=("pass_seq", "sum"),
                     pass_core=("pass_core", "sum"), gate=("gate_pass", "sum"), ch18=("channel18_ok", "sum"))
    for t in (1, 2, 3):
        summary[f"tier{t}"] = by.gate_tier.apply(lambda x, t=t: int((x == t).sum()))
    print(summary.to_string())
    s = strictness(d)
    s.to_csv(a.out.replace(".csv", "_strictness.csv"), index=False)
    print(s.pivot_table(index="criterion", columns=["motif", "arm"], values="pass_rate").round(2).to_string())


if __name__ == "__main__":
    main()
