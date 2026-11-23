"""Known-geometry recovery -- plan 4.2 acceptance criterion 3.

"The geometric metrics computed on 7K4Z must match values measured directly
from the deposited coordinates. This is a unit test, not a model test."

The reference values come from src/derive_atom_map.py, which selects atoms and
computes distances through an entirely separate code path. So this is a real
cross-check of geometry.py, not a tautology against its own output.

Run:  <env-rfd3 python> -m pytest tests/ -v
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest
from biotite.structure.io.pdbx import CIFFile, get_structure

from kemp import geometry as g

import os
# KEMP_ROOT lets the tests run from a workspace made by scripts/make_workspace.sh (benchmark/ lives there, not in the repo)
ROOT = Path(os.environ["KEMP_ROOT"]) if os.environ.get("KEMP_ROOT") else Path(__file__).resolve().parents[1]
ATOM_MAP = ROOT / "configs" / "atom_map.yaml"
PDB_DIR = ROOT / "benchmark" / "pdb"

# Tolerance: these are the same coordinates through different code, so they
# should agree to within printing precision of the reference table (2 dp).
TOL = 0.015


@pytest.fixture(scope="module")
def amap():
    return g.load_atom_map(ATOM_MAP)


@pytest.fixture(scope="module")
def structures():
    out = {}
    for cif in sorted(PDB_DIR.glob("*.cif")):
        out[cif.stem] = get_structure(CIFFile.read(cif), model=1)
    return out


# (pdb, lig_chain, lig_resid, base_dist, donor_dist, wall_dist)
# Measured by src/derive_atom_map.py over chain A of each structure.
KNOWN = [
    ("7K4Z", "A", 402, 2.51, 3.10, 5.54),   # HG3.17   126,000 M-1s-1
    ("7K4X", "A", 401, 2.49, 2.96, 5.01),   # HG3.7     27,000
    ("8RD5", "A", 401, 2.57, 3.03, 4.79),   # HG3.R5
    ("7K4Q", "A", 401, 2.25, 6.05, 5.02),   # HG3          146  (Lys50)
    ("6SS3", "A", 401, 2.52, 7.37, 3.83),   # HG3.17 Q50K  (ablated)
    ("6SS1", "A", 401, 2.61, 99.00, 5.05),  # HG3.17 Q50A  (no donor at all)
]

BASE_RESID, DONOR_RESID, WALL_RESID = 127, 50, 44


@pytest.mark.parametrize("pdb,lc,lr,base,donor,wall", KNOWN)
def test_known_geometry_recovery(structures, amap, pdb, lc, lr, base, donor, wall):
    """geometry.py reproduces distances measured from deposited coordinates."""
    s = structures[pdb]

    got_base = g.base_to_abstraction(
        s, amap, base_resid=BASE_RESID, lig_chain=lc, lig_resid=lr, chain="A"
    )
    assert got_base == pytest.approx(base, abs=TOL), f"{pdb} base"

    got_donor = g.donor_to_oxyanion(
        s, amap, donor_resid=DONOR_RESID, lig_chain=lc, lig_resid=lr, chain="A"
    )
    assert got_donor == pytest.approx(donor, abs=TOL), f"{pdb} donor"

    got_wall = g.wall_stacking(
        s, amap, wall_resid=WALL_RESID, lig_chain=lc, lig_resid=lr, chain="A"
    )
    assert got_wall == pytest.approx(wall, abs=TOL), f"{pdb} wall"


def test_atom_map_roles_are_the_derived_ones(amap):
    """The empirically derived mapping, pinned so a config edit cannot drift."""
    assert amap.ligand_ccd_id == "6NT"
    assert amap.abstraction_centre == "N3"   # carries the acidic H, near Asp127
    assert amap.oxyanion_centre == "N1"      # the phenolate-O position, near Gln50
    assert amap.nitro_nitrogen == "NO1"
    assert set(amap.nitro_oxygens) == {"O11", "O21"}


def test_atom_map_rejects_degenerate_roles():
    """A config that maps two roles to one atom must not load silently."""
    with pytest.raises(ValueError, match="same atom"):
        g.AtomMap(
            ligand_ccd_id="6NT",
            abstraction_centre="N3",
            oxyanion_centre="N3",
            nitro_nitrogen="NO1",
            nitro_oxygens=("O11", "O21"),
            benzo_ring=(),
            base_tip_atoms=("OD1", "OD2"),
            base_proximal_atom="OD2",
            base_all_atoms=(),
            donating_atom_by_resname={},
            wall_ring_atoms=(),
            raw={},
        )


def test_atom_map_rejects_proximal_not_in_tips():
    with pytest.raises(ValueError, match="not among tip_atoms"):
        g.AtomMap(
            ligand_ccd_id="6NT",
            abstraction_centre="N3",
            oxyanion_centre="N1",
            nitro_nitrogen="NO1",
            nitro_oxygens=("O11", "O21"),
            benzo_ring=(),
            base_tip_atoms=("OD1", "OD2"),
            base_proximal_atom="OG",
            base_all_atoms=(),
            donating_atom_by_resname={},
            wall_ring_atoms=(),
            raw={},
        )


# --- the separation the gate depends on (plan 4.2 criterion 2) -------------


def test_donor_metric_separates_active_from_ablated(structures, amap):
    """Gln50 variants must beat the Lys/Ala ablations on the donor metric.

    This is the half of criterion 2 that needs no physics scorer, so it holds
    even while tmol is unavailable.
    """
    def donor(pdb, lc, lr):
        return g.donor_to_oxyanion(
            structures[pdb], amap, donor_resid=DONOR_RESID,
            lig_chain=lc, lig_resid=lr, chain="A",
        )

    active = [donor("7K4Z", "A", 402), donor("7K4X", "A", 401), donor("8RD5", "A", 401)]
    ablated = [donor("7K4Q", "A", 401), donor("6SS3", "A", 401), donor("6SS1", "A", 401)]

    assert max(active) < min(ablated), (
        f"donor metric does not separate: active={active} ablated={ablated}"
    )
    # And the separation is wide, not marginal.
    assert min(ablated) - max(active) > 2.0


def test_base_metric_does_not_track_activity(structures, amap):
    """Guard rail on a real trap, documented in thresholds.yaml.

    HG3 (146 M-1s-1) has a SHORTER base distance than HG3.17 (126,000). Any
    composite that weights base distance heavily will anti-correlate with
    activity. If this test ever fails, re-read docs/MILESTONE1.md before
    "fixing" it -- the finding, not the test, is what matters.
    """
    def base(pdb, lc, lr):
        return g.base_to_abstraction(
            structures[pdb], amap, base_resid=BASE_RESID,
            lig_chain=lc, lig_resid=lr, chain="A",
        )

    assert base("7K4Q", "A", 401) < base("7K4Z", "A", 402)


# --- primitives ------------------------------------------------------------


def test_kabsch_rmsd_is_zero_for_rigid_motion():
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(8, 3))
    theta = 0.7
    rot = np.array(
        [
            [math.cos(theta), -math.sin(theta), 0.0],
            [math.sin(theta), math.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    moved = pts @ rot.T + np.array([10.0, -4.0, 2.5])
    assert g.kabsch_rmsd(moved, pts) == pytest.approx(0.0, abs=1e-9)


def test_kabsch_rmsd_rejects_reflection():
    """A mirrored motif is NOT a match; the reflection guard must catch it."""
    rng = np.random.default_rng(1)
    pts = rng.normal(size=(10, 3)) * 3.0
    mirrored = pts * np.array([1.0, 1.0, -1.0])
    assert g.kabsch_rmsd(mirrored, pts) > 0.5


def test_kabsch_rmsd_needs_three_atoms():
    with pytest.raises(g.MetricUnavailable, match=">=3 atoms"):
        g.kabsch_rmsd(np.zeros((2, 3)), np.zeros((2, 3)))


def test_kabsch_rmsd_shape_mismatch_raises():
    with pytest.raises(g.MetricUnavailable, match="shape mismatch"):
        g.kabsch_rmsd(np.zeros((5, 3)), np.zeros((4, 3)))


def test_missing_atoms_fail_loud(structures, amap):
    """Plan 0.2: a metric that cannot be computed raises, never returns 0.0."""
    with pytest.raises(g.MetricUnavailable, match="no atoms found"):
        g.base_to_abstraction(
            structures["7K4Z"], amap, base_resid=99999,
            lig_chain="A", lig_resid=402, chain="A",
        )


def test_absent_donor_returns_sentinel_not_error(structures, amap):
    """A 'no oxyanion donor' cell is a design choice, not a failure."""
    assert g.donor_to_oxyanion(
        structures["7K4Z"], amap, donor_resid=None,
        lig_chain="A", lig_resid=402,
    ) == 99.0


def test_specified_but_missing_donor_still_raises(structures, amap):
    """The sentinel must not mask a genuinely broken structure."""
    with pytest.raises(g.MetricUnavailable, match="absent from the structure"):
        g.donor_to_oxyanion(
            structures["7K4Z"], amap, donor_resid=99999,
            lig_chain="A", lig_resid=402,
        )


def test_no_steric_clash_in_natives(structures, amap):
    """The [lit] 1.5 A clash rule must pass on every deposited structure."""
    for pdb, lc, lr, *_ in KNOWN:
        d = g.min_protein_ligand_contact(
            structures[pdb], amap, lig_chain=lc, lig_resid=lr
        )
        assert d > 1.5, f"{pdb} has a {d:.2f} A protein-ligand contact"


def test_ligand_copies_discovered(structures, amap):
    assert g.iter_ligand_copies(structures["7K4Z"], amap) == [("A", 402)]
    assert len(g.iter_ligand_copies(structures["7K4X"], amap)) == 2


# --- the motif gate (regression guard on a real silent failure) -----------
#
# RFD4 produced 10/10 designs whose catalytic residues sat 8.1-23.0 A from the
# ligand while every conditioning column looked correct and nothing warned.
# Only measuring the output caught it. These tests pin the gate that catches
# it, so the pipeline can never again accept a design with no active site.


def _gate(**scores):
    import sys as _sys
    from pathlib import Path as _P
    _sys.path.insert(0, str(_P(__file__).resolve().parents[1] / "src"))
    from extract_designs import apply_motif_gate
    from kemp.record import DesignRecord

    thr = g.load_thresholds(ROOT / "configs" / "thresholds.yaml")
    rec = DesignRecord(design_id="t", run_id="t")
    rec.length = scores.pop("length", 200)
    rec.add_scores("geom", scores)
    return apply_motif_gate(rec, thr)


def test_motif_gate_accepts_native_like_geometry():
    """The HG3.17 crystal values must pass."""
    passed, reasons = _gate(
        base_to_abstraction=2.51, donor_to_oxyanion=3.10,
        wall_stacking=5.54, min_contact=2.51,
    )
    assert passed, reasons


def test_motif_gate_rejects_the_observed_failure():
    """Real measured values from the broken-spec run must be rejected."""
    passed, reasons = _gate(
        base_to_abstraction=13.08, donor_to_oxyanion=11.95,
        wall_stacking=3.90, min_contact=3.90,
    )
    assert not passed
    assert any("base->abstraction" in r for r in reasons)


def test_motif_gate_rejects_a_clash():
    passed, reasons = _gate(
        base_to_abstraction=1.1, donor_to_oxyanion=3.10,
        wall_stacking=5.0, min_contact=1.1,
    )
    assert not passed
    assert any("clash" in r or "contact" in r for r in reasons)


def test_motif_gate_rejects_over_length():
    """The competition cap is an eligibility rule, not a preference."""
    passed, reasons = _gate(
        base_to_abstraction=2.5, donor_to_oxyanion=3.1,
        wall_stacking=5.0, min_contact=2.5, length=261,
    )
    assert not passed
    assert any("length" in r for r in reasons)


# --- ligand burial (the surface-exposure failure mode) -------------------
#
# A 32-design batch satisfied EVERY catalytic distance while 14/24 had the
# ligand on the protein surface (rel_sasa up to 0.54 vs 0.008 in crystals).
# Distance metrics alone cannot see this; these tests pin the burial check.


def test_crystal_ligand_is_buried(structures, amap):
    """The reference crystals define what a real pocket looks like."""
    for pdb, lc, lr in (("7K4Z", "A", 402), ("7K4X", "A", 401)):
        b = g.ligand_burial(structures[pdb], amap, lig_chain=lc, lig_resid=lr)
        assert b["rel_sasa"] < 0.05, f"{pdb} rel_sasa {b['rel_sasa']}"
        assert b["enclosure"] > 0.95, f"{pdb} enclosure {b['enclosure']}"
        assert b["n_contacts"] >= 60, f"{pdb} contacts {b['n_contacts']}"


def test_crystals_pass_the_burial_thresholds(structures, amap):
    """Our cut-offs must not reject the real thing."""
    thr = g.load_thresholds(ROOT / "configs" / "thresholds.yaml")
    bur = thr["geometry"]["ligand_burial"]
    for pdb, lc, lr in (("7K4Z", "A", 402), ("4BS0", "A", 1305)):
        b = g.ligand_burial(structures[pdb], amap, lig_chain=lc, lig_resid=lr)
        assert b["rel_sasa"] <= bur["rel_sasa_pass_max"]
        assert b["enclosure"] >= bur["enclosure_pass_min"]
        assert b["n_contacts"] >= bur["n_contacts_pass_min"]


def test_burial_gate_rejects_surface_ligand():
    """Real measured values from a surface-exposed design must be rejected."""
    passed, reasons = _gate(
        base_to_abstraction=2.45, donor_to_oxyanion=2.94, wall_stacking=5.46,
        min_contact=2.45, rel_sasa=0.539, enclosure=0.527, n_contacts=14.0,
    )
    assert not passed
    assert any("rel_sasa" in r or "enclosure" in r or "contacts" in r
               for r in reasons), reasons


def test_burial_gate_accepts_buried_design():
    """And a genuinely buried design must pass."""
    passed, reasons = _gate(
        base_to_abstraction=2.57, donor_to_oxyanion=3.07, wall_stacking=5.47,
        min_contact=2.57, rel_sasa=0.007, enclosure=1.000, n_contacts=52.0,
    )
    assert passed, reasons


def test_donor_floor_rejects_strained_contact():
    """A 2.41 A donor-oxyanion contact is strained, not a better H-bond.

    Real value from burial-conditioned run job21, design 47ec9b0b4620. No
    catalytically active crystal goes below 2.96 A. This passed the original
    generic 2.40 A clash floor by 0.01 A.
    """
    passed, reasons = _gate(
        base_to_abstraction=2.50, donor_to_oxyanion=2.41, wall_stacking=5.45,
        min_contact=2.41, rel_sasa=0.007, enclosure=1.000, n_contacts=52.0,
    )
    assert not passed
    assert any("clash" in r or "hard_fail_min" in r for r in reasons), reasons


def test_donor_floor_still_accepts_crystal_values():
    """The tightened floor must not reject the real active geometries."""
    for d in (2.96, 2.98, 3.02, 3.03, 3.10):
        passed, reasons = _gate(
            base_to_abstraction=2.53, donor_to_oxyanion=d, wall_stacking=5.53,
            min_contact=2.53, rel_sasa=0.011, enclosure=0.992, n_contacts=66.0,
        )
        assert passed, f"donor {d} rejected: {reasons}"
