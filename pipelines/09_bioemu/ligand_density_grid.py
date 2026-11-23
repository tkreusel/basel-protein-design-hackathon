#!/usr/bin/env python3
"""
ligand_density_grid.py - compare ligand density across many BioEmu runs.

    python ligand_density_grid.py runs/
    python ligand_density_grid.py runs/ --view side --ncols 4 -o grid.png

`runs/` is a folder of run folders, each holding a topology (topology.cif or
.pdb) and a trajectory (samples.xtc). Sub-folders are searched recursively.

Every run is superposed on the SAME reference protein (the first run, or
--ref-run NAME) and drawn in the same orientation on the same axes, so the
panels can be compared directly. Shaded bands hold 50 / 80 / 95 % of that run's
ligand positions; the grey trace is that run's own mean Calpha backbone.

Panels are ordered by ligand RMSF, lowest (most tightly held) first, so the
grid reads in the same order as ligand_rmsf.py.

Requires: pip install mdtraj gemmi scipy matplotlib numpy
"""
import argparse, os, warnings
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

ap = argparse.ArgumentParser()
ap.add_argument("folder", help="folder containing one sub-folder per run")
ap.add_argument("--ligand", help="residue name (e.g. 6NT) or mdtraj selection")
ap.add_argument("--view", default="front", choices=["front", "side"])
ap.add_argument("--ncols", type=int, help="panels per row (default: about square)")
ap.add_argument("--ref-run", help="run whose orientation everything is aligned to "
                                  "(default: the first one found)")
ap.add_argument("--sort", default="rmsf", choices=["rmsf", "name"])
ap.add_argument("--max-runs", type=int, default=40)
ap.add_argument("-o", "--out", default="ligand_density_grid.png")
a = ap.parse_args()

SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
BANDS = ["#86b6ef", "#2a78d6", "#104281"]            # 95 %, 80 %, 50 %
SILHOUETTE = "#edebe3"
DEPTH = LinearSegmentedColormap.from_list("depth", ["#bbb8ad", "#4a4843"])


# --- mmCIF -> mdtraj Topology (mdtraj cannot read mmCIF itself)
def cif_topology(path):
    import gemmi

    def element_of(symbol, atom_name):
        for s in (symbol, (atom_name or "").lstrip("0123456789")[:1]):
            if s:
                s = s.capitalize() if len(s) > 1 else s.upper()
                try:
                    return md.element.get_by_symbol(s)
                except KeyError:
                    pass
        return md.element.virtual

    doc = gemmi.cif.read(str(path))
    for block in doc:
        def col(tag, default=None):
            c = block.find_values("_atom_site." + tag)
            if c is None or len(c) == 0:
                return default
            return [gemmi.cif.as_string(v) for v in c]

        x = col("Cartn_x")
        if not x:
            continue
        n = len(x)
        name = col("label_atom_id") or col("auth_atom_id") or [""] * n
        comp = col("label_comp_id") or col("auth_comp_id") or ["UNK"] * n
        asym = col("label_asym_id") or col("auth_asym_id") or ["A"] * n
        seq = col("auth_seq_id") or col("label_seq_id") or [""] * n
        ins = col("pdbx_PDB_ins_code", [""] * n)
        sym = col("type_symbol", [""] * n)
        model = col("pdbx_PDB_model_num")
        rows = range(n)
        if model:
            rows = [i for i in rows if model[i] == model[0]]
        top = md.Topology()
        chain = residue = None
        chain_key = res_key = object()
        counter = 0
        for i in rows:
            if asym[i] != chain_key:
                chain = top.add_chain()
                chain_key, res_key = asym[i], object()
            key = (asym[i], seq[i], ins[i], comp[i])
            if key != res_key:
                counter += 1
                try:
                    num = int(seq[i])
                except (TypeError, ValueError):
                    num = counter
                residue = top.add_residue(comp[i] or "UNK", chain, resSeq=num)
                res_key = key
            top.add_atom(name[i] or f"X{i}", element_of(sym[i], name[i]), residue)
        if top.n_atoms:
            return top
    raise ValueError(f"no atom coordinates in {path}")


def find_runs(root):
    out = []
    for dirpath, _, files in os.walk(root):
        xtcs = sorted(f for f in files if f.endswith((".xtc", ".dcd", ".trr")))
        tops = ([f for f in files if f in ("topology.cif", "topology.pdb")]
                or sorted(f for f in files if f.endswith((".cif", ".pdb"))))
        if not xtcs or not tops:
            continue
        xtc = "samples.xtc" if "samples.xtc" in xtcs else xtcs[0]
        rel = os.path.relpath(dirpath, root)
        out.append((rel if rel != "." else os.path.basename(os.path.abspath(root)),
                    os.path.join(dirpath, xtc), os.path.join(dirpath, tops[0])))
    return sorted(out)


def ligand_atoms(topology, spec):
    if spec and " " not in spec.strip():
        idx = [x.index for x in topology.atoms
               if x.residue.name == spec.strip() and x.element.symbol != "H"]
        if not idx:
            raise ValueError(f"no residue named {spec}")
        return np.array(idx)
    if spec:
        idx = topology.select(f"({spec}) and not element H")
        if len(idx) == 0:
            raise ValueError(f"selection {spec!r} matched no atoms")
        return idx
    res = [r for r in topology.residues if not r.is_protein and not r.is_water]
    if not res:
        raise ValueError("no non-protein residue found (pass --ligand)")
    r = max(res, key=lambda r: r.n_atoms)
    return np.array([x.index for x in r.atoms if x.element.symbol != "H"])


def smooth_path(P, up=5, passes=2):
    if len(P) < 4:
        return P
    g = np.linspace(0, len(P) - 1, (len(P) - 1) * up + 1)
    Q = np.column_stack([np.interp(g, np.arange(len(P)), P[:, c]) for c in range(P.shape[1])])
    w = max(2, up // 2)
    kern = np.ones(2 * w + 1) / (2 * w + 1)
    for _ in range(passes):
        Q = np.column_stack([np.convolve(np.pad(Q[:, c], (w, w), mode="edge"), kern, "valid")
                             for c in range(Q.shape[1])])
    return Q



def clip(x, n=24):
    return x if len(x) <= n else x[:n - 1] + "…"

def trim_shared(names):
    """drop the prefix/suffix every run name shares, so only the difference shows"""
    if len(names) < 2:
        return list(names), "", ""
    pre = os.path.commonprefix(names)
    suf = os.path.commonprefix([n[::-1] for n in names])[::-1]
    for cut_suf in (suf, ""):
        for cut_pre in (pre, ""):
            out = [n[len(cut_pre):len(n) - len(cut_suf) or None] for n in names]
            if all(o.strip(" _-.") for o in out) and len(set(out)) == len(set(names)):
                return out, cut_pre, cut_suf
    return list(names), "", ""


def ca_key(topology, idx):
    return [(topology.atom(i).residue.chain.index, topology.atom(i).residue.resSeq)
            for i in idx]


# --- load every run
runs = find_runs(a.folder)
if not runs:
    raise SystemExit(f"No runs found under {a.folder}")
print(f"Found {len(runs)} run folder(s) under {a.folder}\n")

loaded = []
for name, xtc, top_path in runs:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            top = (cif_topology(top_path) if top_path.lower().endswith((".cif", ".mmcif"))
                   else md.load_topology(top_path))
            t = md.load(xtc, top=top)
        lig = ligand_atoms(t.topology, a.ligand)
        ca = t.topology.select("protein and name CA")
        if t.n_frames < 2 or len(ca) < 3:
            raise ValueError(f"{t.n_frames} frames, {len(ca)} Calpha")
        loaded.append(dict(name=name, t=t, lig=lig, ca=ca,
                           ligname=t.topology.atom(lig[0]).residue.name))
        print(f"  {name:42s} {t.n_frames:4d} samples, ligand "
              f"{loaded[-1]['ligname']} ({len(lig)} atoms)")
    except Exception as e:
        print(f"  {name:42s} SKIPPED - {e}")
if not loaded:
    raise SystemExit("No run could be analysed.")
if len(loaded) > a.max_runs:
    print(f"\nKeeping the first {a.max_runs} runs (--max-runs to change)")
    loaded = loaded[:a.max_runs]

# --- common reference frame: superpose every run on one reference protein
ref = next((r for r in loaded if r["name"] == a.ref_run), loaded[0])
if a.ref_run and ref["name"] != a.ref_run:
    print(f"\nwarning: --ref-run {a.ref_run!r} not found, using {ref['name']!r}")
ref_mean = md.Trajectory(ref["t"].xyz.mean(0, keepdims=True), ref["t"].topology)
ref_ca, ref_keys = ref["ca"], ca_key(ref["t"].topology, ref["ca"])
key_pos = {k: i for i, k in enumerate(ref_keys)}

for r in loaded:
    keys = ca_key(r["t"].topology, r["ca"])
    if len(keys) == len(ref_keys):
        sel, rsel = r["ca"], ref_ca                       # same protein, 1:1
    else:                                                 # match by chain + residue number
        pairs = [(c, ref_ca[key_pos[k]]) for c, k in zip(r["ca"], keys) if k in key_pos]
        if len(pairs) < 3:
            print(f"  note: {r['name']} shares no numbering with the reference; "
                  f"drawn in its own frame")
            r["t"].superpose(r["t"], 0, atom_indices=r["ca"])
            continue
        sel, rsel = (np.array(x) for x in zip(*pairs))
    r["t"].superpose(ref_mean, 0, atom_indices=sel, ref_atom_indices=rsel)

# principal axes of the reference, so all panels share an orientation
P = ref_mean.xyz[0, ref_ca] * 10.0
centre = P.mean(0)
_, _, axes3 = np.linalg.svd(P - centre)
i, j, k = (0, 1, 2) if a.view == "front" else (0, 2, 1)

for r in loaded:
    xyz = r["t"].xyz * 10.0
    r["L"] = ((xyz[:, r["lig"]].reshape(-1, 3) - centre) @ axes3.T)
    r["trace"] = smooth_path((xyz[:, r["ca"]].mean(0) - centre) @ axes3.T)
    Lp = xyz[:, r["lig"]]                                  # already Angstrom
    r["rmsf"] = float(np.sqrt(((Lp - Lp.mean(0)) ** 2).sum(-1).mean(0)).mean())

loaded.sort(key=(lambda r: r["rmsf"]) if a.sort == "rmsf" else (lambda r: r["name"]))

# --- shared limits so every panel is on the same scale
allpts = np.vstack([np.vstack([r["trace"][:, [i, j]], r["L"][:, [i, j]]]) for r in loaded])
lo, hi = allpts.min(0) - 5, allpts.max(0) + 5
span = (hi - lo).max()                                     # square panels, equal scale
mid = (hi + lo) / 2
lo, hi = mid - span / 2, mid + span / 2
gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], 200), np.linspace(lo[1], hi[1], 200))
grid = np.vstack([gx.ravel(), gy.ravel()])

# --- plot
short_names, cut_pre, cut_suf = trim_shared([r["name"] for r in loaded])
n = len(loaded)
ncols = a.ncols or min(5, max(1, int(np.ceil(np.sqrt(n)))))
nrows = int(np.ceil(n / ncols))
panel = 2.7
sub = (f"{n} runs, all superposed on '{clip(ref['name'], 34)}' and drawn on the same "
       f"axes ({span:.0f} × {span:.0f} Å).\nOrdered by ligand RMSF, tightest first.")
if cut_pre or cut_suf:
    sub += f" Names shown without the shared '{clip(cut_pre)}…{clip(cut_suf)}'."
nlines = sub.count("\n") + 1
header = 0.62 + 0.21 * nlines + 0.52
figh = panel * nrows + header + 0.3
fig, axs = plt.subplots(nrows, ncols, figsize=(panel * ncols, figh), squeeze=False)
fig.patch.set_facecolor(SURF)
fig.subplots_adjust(left=0.012, right=0.988, top=1 - header / figh, bottom=0.012,
                    wspace=0.06, hspace=0.16)

for ax in axs.ravel():
    ax.set_facecolor(SURF)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

for p, r in enumerate(loaded):
    ax = axs[p // ncols][p % ncols]
    tr = r["trace"]

    pk = gaussian_kde(tr[:, [i, j]].T)
    pk.set_bandwidth(pk.factor * 0.45)
    pd = pk(grid).reshape(gx.shape)
    lvl = np.quantile(pk(tr[:, [i, j]].T), 0.02)
    ax.contourf(gx, gy, pd, levels=[lvl, pd.max() * 1.01], colors=[SILHOUETTE], zorder=0.5)

    pts = r["L"][:, [i, j]].T
    try:
        kde = gaussian_kde(pts)
        dens = kde(grid).reshape(gx.shape)
        at = kde(pts)
        levels = [np.quantile(at, 1 - q) for q in (0.95, 0.80, 0.50)] + [dens.max() * 1.01]
        if np.all(np.diff(levels) > 0):
            ax.contourf(gx, gy, dens, levels=levels, colors=BANDS, alpha=0.88, zorder=1)
            ax.contour(gx, gy, dens, levels=levels[:3], colors=SURF, linewidths=0.6,
                       zorder=1.3)
        else:
            raise np.linalg.LinAlgError
    except (np.linalg.LinAlgError, ValueError):            # ligand barely moves
        ax.scatter(*pts, s=6, color=BANDS[2], lw=0, zorder=2)

    xy, depth = tr[:, [i, j]], tr[:, k]
    ax.plot(xy[:, 0], xy[:, 1], color=SURF, lw=2.8, zorder=2, alpha=0.9,
            solid_capstyle="round", solid_joinstyle="round")
    u = (0.5 * (depth[:-1] + depth[1:]) - depth.min()) / (np.ptp(depth) + 1e-9)
    ax.add_collection(LineCollection(np.stack([xy[:-1], xy[1:]], axis=1), colors=DEPTH(u),
                                     linewidths=0.6 + 1.1 * u, zorder=2.5,
                                     capstyle="round", joinstyle="round"))

    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
    label = short_names[p] if len(short_names[p]) <= 34 else \
        short_names[p][:16] + "…" + short_names[p][-17:]
    ax.set_title(label, loc="left", fontsize=8.5, fontweight="bold", color=INK, pad=3)
    ax.text(0.015, 0.985, f"RMSF {r['rmsf']:.2f} Å · {r['t'].n_frames} samples"
            + (f" · {r['ligname']}" if len({x['ligname'] for x in loaded}) > 1 else ""),
            transform=ax.transAxes, va="top", fontsize=7.5, color=INK2)

for p in range(n, nrows * ncols):
    axs[p // ncols][p % ncols].set_visible(False)

fig.suptitle(f"Where the ligand sits in each run — {a.view} view", x=0.012, ha="left",
             fontsize=13, fontweight="bold", color=INK, y=1 - 0.33 / figh)
fig.text(0.012, 1 - 0.60 / figh, sub, color=INK2, fontsize=9, va="top", linespacing=1.6)
fig.legend(handles=[Patch(color=BANDS[2], label="50% of ligand positions"),
                    Patch(color=BANDS[1], label="80%"), Patch(color=BANDS[0], label="95%"),
                    plt.Line2D([], [], color="#4a4843", lw=1.8, label="Cα trace (mean)")],
           loc="center left",
           bbox_to_anchor=(0.012, 1 - (0.62 + 0.21 * nlines + 0.22) / figh),
           ncol=4, frameon=False, fontsize=8.5, labelcolor=INK2, handlelength=1.3)
fig.savefig(a.out, dpi=200)
print(f"\nSaved {a.out}")
