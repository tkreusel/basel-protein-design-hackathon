"""Theozyme constraints: per-constraint deviations + (optionally) a tmol constraint term.

Two sources of constraints, same downstream format:
  1. config/theozyme_reference.yaml  - measured from the crystal TSA complexes (default)
  2. a Rosetta enzdes-style .cst file (parse_enzdes_cst) - YOUR theozyme; not available when this
     was written, so the parser follows the documented block format and must be confirmed against
     the real file (see PIPELINE.md, open question 1).

Reporting: for each constraint  dev = value - target (circularly wrapped for angles),
z = dev/tol, pass = |z| <= 1; plus n_pass, rms_z and a flat-bottomed penalty sum(max(0,|z|-1)^2).
The same flat-bottom form is what the tmol constraint term uses, so "constraint energy" and
"deviation report" tell one story.
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from . import geometry as G
from . import ligands as L


@dataclass
class Constraint:
    key: str           # metric key, e.g. 'geo.d_O_X', or 'cst:<block>:<dof>' for .cst files
    target: float
    tol: float
    kind: str = "linear"   # linear | circular
    force_k: float = 1.0
    atoms: tuple | None = None   # symbolic atoms for the tmol term, see ATOM_SPEC

    @property
    def short(self) -> str:
        return self.key.replace("geo.", "").replace(":", "_")


# symbolic atom definitions for the tmol constraint term (role names resolved per pose)
#  ('P', role) -> protein atom of site['base'] ; ('L', role) -> ligand atom by ROLES key / name
ATOM_SPEC = {
    "geo.d_O_X":          (("P", "O"), ("L", "abstraction")),
    "geo.d_O_H":          (("P", "O"), ("L", "abstraction_H")),
    "geo.ang_O_H_X":      (("P", "O"), ("L", "abstraction_H"), ("L", "abstraction")),
    "geo.ang_C_O_H":      (("P", "C"), ("P", "O"), ("L", "abstraction_H")),
    "geo.tors_syn":       (("P", "O_other"), ("P", "C"), ("P", "O"), ("L", "abstraction_H")),
    "geo.tors_C_O_H_X":   (("P", "C"), ("P", "O"), ("L", "abstraction_H"), ("L", "abstraction")),
    "geo.tors_O_H_X_C3A": (("P", "O"), ("L", "abstraction_H"), ("L", "abstraction"), ("L", "ring_attached")),
}


def load_reference(path: Path, group: str = "hg3_lineage", drop: tuple = ()) -> list[Constraint]:
    ref = yaml.safe_load(Path(path).read_text())["groups"][group]["constraints"]
    out = []
    for k, s in ref.items():
        if k in drop:
            continue
        out.append(Constraint(k, s["mean"], s["tol"], s["kind"], s.get("force_k", 1.0), ATOM_SPEC.get(k)))
    return out


def evaluate(metrics: dict, constraints: list[Constraint], prefix: str = "cst.") -> dict:
    out, zs, pen, npass, n = {}, [], 0.0, 0, 0
    for c in constraints:
        v = metrics.get(c.key)
        if v is None or not np.isfinite(v):
            out[f"{prefix}{c.short}.dev"] = float("nan")
            continue
        dev = G.wrap180(v - c.target) if c.kind == "circular" else v - c.target
        z = dev / c.tol
        out[f"{prefix}{c.short}.dev"] = dev
        out[f"{prefix}{c.short}.z"] = z
        out[f"{prefix}{c.short}.pass"] = float(abs(z) <= 1.0)
        zs.append(z)
        pen += max(0.0, abs(z) - 1.0) ** 2
        npass += abs(z) <= 1.0
        n += 1
    out[f"{prefix}n_total"] = float(n)
    out[f"{prefix}n_pass"] = float(npass)
    out[f"{prefix}frac_pass"] = npass / n if n else float("nan")
    out[f"{prefix}rms_z"] = float(np.sqrt(np.mean(np.square(zs)))) if zs else float("nan")
    out[f"{prefix}penalty"] = float(pen)
    return out


# ------------------------------------------------------------------ tmol constraint term
def _flat(dev, tol, scale=1.0):
    return scale * (torch_relu(dev.abs() - tol) / tol) ** 2


def torch_relu(x):
    import torch
    return torch.clamp(x, min=0.0)


def flat_distance(atoms, params):  # params: target, tol, scale
    import torch
    d = torch.linalg.norm(atoms[:, 0] - atoms[:, 1], dim=-1)
    return _flat(d - params[:, 0], params[:, 1], params[:, 2])


def _angle(a, b, c):
    import torch
    v1, v2 = a - b, c - b
    cross = torch.linalg.norm(torch.cross(v1, v2, dim=-1), dim=-1)
    dot = (v1 * v2).sum(-1)
    return torch.atan2(cross, dot)


def flat_angle(atoms, params):  # radians
    return _flat(_angle(atoms[:, 0], atoms[:, 1], atoms[:, 2]) - params[:, 0], params[:, 1], params[:, 2])


def _torsion(a, b, c, d):
    import torch
    b0, b1, b2 = a - b, c - b, d - c
    b1n = b1 / torch.linalg.norm(b1, dim=-1, keepdim=True)
    v = b0 - (b0 * b1n).sum(-1, keepdim=True) * b1n
    w = b2 - (b2 * b1n).sum(-1, keepdim=True) * b1n
    x = (v * w).sum(-1)
    y = (torch.cross(b1n, v, dim=-1) * w).sum(-1)
    return torch.atan2(y, x)


def flat_torsion(atoms, params):
    import torch
    t = _torsion(atoms[:, 0], atoms[:, 1], atoms[:, 2], atoms[:, 3])
    d = torch.remainder(t - params[:, 0] + math.pi, 2 * math.pi) - math.pi
    return _flat(d, params[:, 1], params[:, 2])


def add_to_pose(pose, site, constraints: list[Constraint], block_map: dict, lig_resname: str = "6NT",
                scale: float = 1.0):
    """Attach flat-bottomed constraints (those with ATOM_SPEC) to the pose's ConstraintSet.

    Returns (pose, list of constraint keys applied). Constraints whose atoms cannot be resolved
    (e.g. no H on the ligand, no base found) are skipped and reported, never silently dropped.
    """
    import torch
    from . import tmol_core as TC

    roles = L.ROLES[lig_resname]
    base = site.get("base")
    if base is None:
        return pose, [], [c.key for c in constraints if c.atoms]
    bres, lblk = block_map["res"][base["resid"]], block_map["lig"]

    def resolve(a):
        kind, role = a
        if kind == "P":
            return bres, pose.block_type(0, bres).atom_to_idx[base[role]]
        name = roles.get(role, role)
        if name is None:
            raise KeyError(role)
        return lblk, pose.block_type(0, lblk).atom_to_idx[name]

    groups = {2: ([], [], flat_distance), 3: ([], [], flat_angle), 4: ([], [], flat_torsion)}
    applied, skipped = [], []
    for c in constraints:
        if not c.atoms:
            continue
        try:
            idx = [resolve(a) for a in c.atoms]
        except KeyError:
            skipped.append(c.key)
            continue
        ats, prm, _ = groups[len(c.atoms)]
        ats.append([[0, b, a] for b, a in idx])
        if len(c.atoms) == 2:
            prm.append([c.target, c.tol, scale])
        else:  # degrees -> radians
            prm.append([math.radians(c.target), math.radians(c.tol), scale])
        applied.append(c.key)
    for n, (ats, prm, fn) in groups.items():
        if ats:
            pose = TC.with_constraints(pose, fn, torch.tensor(ats, dtype=torch.int32),
                                       torch.tensor(prm, dtype=torch.float32))
    return pose, applied, skipped


# ------------------------------------------------------------------ Rosetta enzdes-style .cst
CST_TYPES_B1 = {  # Rosetta atom_type -> residue -> protein atoms (subset; extend when the real file is known)
    "OOC": {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]},
    "Nhis": {"HIS": ["ND1", "NE2"]},
    "NH2O": {"ASN": ["ND2"], "GLN": ["NE2"]},
    "OH": {"SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"]},
    "Nlys": {"LYS": ["NZ"]},
    "Narg": {"ARG": ["NE", "NH1", "NH2"]},
}
ONE2THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY",
             "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
             "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
DOFS = ["distanceAB", "angle_A", "angle_B", "torsion_A", "torsion_AB", "torsion_B"]


def parse_enzdes_cst(text: str) -> list[dict]:
    """Parse CST::BEGIN ... CST::END blocks.

    Returns per block: {'ligand': {'atoms': [A1,A2,A3], 'res3': [...]},
                        'protein': {'atoms': [B1,B2,B3] | None, 'types': [...] | None, 'res3': [...]},
                        'dofs': {name: (target, tol, k)}}
    """
    blocks = []
    for raw in re.findall(r"CST::BEGIN(.*?)CST::END", text, flags=re.S):
        b = {"ligand": {"atoms": None, "res3": []}, "protein": {"atoms": None, "types": None, "res3": []}, "dofs": {}}
        for line in raw.splitlines():
            line = line.split("#")[0].strip()
            m = re.match(r"TEMPLATE::\s*ATOM_MAP:\s*(\d)\s+(atom_name|atom_type|residue3|residue1):\s*(.*)", line)
            if m:
                which, field, rest = int(m.group(1)), m.group(2), m.group(3)
                vals = [v for v in re.split(r"[,\s]+", rest.strip()) if v]
                tgt = b["ligand"] if which == 1 else b["protein"]
                if field == "atom_name":
                    tgt["atoms"] = vals
                elif field == "atom_type":
                    tgt["types"] = vals
                elif field == "residue3":
                    tgt["res3"] += vals
                else:
                    tgt["res3"] += [ONE2THREE[c] for c in "".join(vals)]
                continue
            m = re.match(r"CONSTRAINT::\s*(distanceAB|angle_A|angle_B|torsion_A|torsion_AB|torsion_B):\s*(.*)", line)
            if m:
                f = [float(x) for x in m.group(2).split()]
                b["dofs"][m.group(1)] = (f[0], f[1], f[2] if len(f) > 2 else 1.0)
        blocks.append(b)
    return blocks


def evaluate_enzdes_cst(rep, blocks: list[dict], alias: dict | None = None) -> dict:
    """Per-constraint deviations of `rep` (replicate, ideally with H) against a parsed .cst.

    Protein residue per block = nearest allowed residue to ligand atom A1 (the .cst has no residue
    numbers). `alias` maps CST ligand atom names -> canonical names (e.g. {'C3': 'N3'}).
    """
    alias = alias or {}
    ly = G.lig_xyz(rep)
    P = G.Atoms.of(rep.protein)
    out: dict = {}

    def lig_atom(n):
        n2 = alias.get(n, n)
        if n2 not in ly:
            raise KeyError(f"ligand atom '{n}' (alias '{n2}') not in {sorted(ly)}; provide alias map")
        return ly[n2]

    n_pass = n_tot = 0
    for bi, b in enumerate(blocks, start=1):
        A = [lig_atom(n) for n in b["ligand"]["atoms"][:3]]
        pr = b["protein"]
        allowed = pr["res3"] or list(ONE2THREE.values())
        best = None
        for rid in sorted({r for r, _ in P._idx}):
            rn = P.resname(rid)
            if rn not in allowed:
                continue
            if pr["atoms"]:
                cand = [pr["atoms"]]
                if not all(P.xyz(rid, a) is not None for a in pr["atoms"]):
                    continue
            else:  # atom_type for B1; B2, B3 unsupported here
                t = pr["types"][0]
                if t not in CST_TYPES_B1 or rn not in CST_TYPES_B1[t]:
                    continue
                raise NotImplementedError("atom_type-only protein atoms need a residue-specific B2/B3 table; "
                                          "share the real .cst so it can be added")
            B = [P.xyz(rid, a) for a in cand[0]]
            d = G.dist(A[0], B[0])
            if best is None or d < best[0]:
                best = (d, rid, rn, B)
        if best is None:
            out[f"cst{bi}.found"] = 0.0
            continue
        _, rid, rn, B = best
        out[f"cst{bi}.resid"] = float(rid)
        vals = {"distanceAB": G.dist(A[0], B[0]),
                "angle_A": G.angle(A[1], A[0], B[0]) if len(A) > 1 else None,
                "angle_B": G.angle(A[0], B[0], B[1]) if len(B) > 1 else None,
                "torsion_A": G.dihedral(A[2], A[1], A[0], B[0]) if len(A) > 2 else None,
                "torsion_AB": G.dihedral(A[1], A[0], B[0], B[1]) if len(A) > 1 and len(B) > 1 else None,
                "torsion_B": G.dihedral(A[0], B[0], B[1], B[2]) if len(B) > 2 else None}
        for name, (t, tol, k) in b["dofs"].items():
            v = vals[name]
            if v is None:
                continue
            dev = G.wrap180(v - t) if name.startswith(("angle", "torsion")) else v - t
            out[f"cst{bi}.{name}.value"] = v
            out[f"cst{bi}.{name}.dev"] = dev
            out[f"cst{bi}.{name}.z"] = dev / tol
            ok = abs(dev) <= tol
            out[f"cst{bi}.{name}.pass"] = float(ok)
            n_pass += ok
            n_tot += 1
    out["cst_file.n_total"], out["cst_file.n_pass"] = float(n_tot), float(n_pass)
    return out
