#!/usr/bin/env python3
"""
ligand_rmsf.py - compare how much the ligand moves, across many BioEmu runs.

    python ligand_rmsf.py runs/
    python ligand_rmsf.py runs/ --ligand 6NT -o rmsf.png

`runs/` is a folder of run folders, each holding a topology (topology.cif or
.pdb) and a trajectory (samples.xtc). Sub-folders are searched recursively.

For each run every sample is superposed on the protein, so the number measures
how much the ligand moves *relative to the protein*, not overall drift. Then

    RMSF = for each ligand heavy atom, the root-mean-square distance from that
           atom's own average position, averaged over the atoms.

One bar per run, lowest RMSF first (the most tightly held ligand at the top).
A thin line on each bar spans the lowest to highest single-atom value, so you
can see whether the whole ligand moves together or just part of it.

Writes the plot and a .csv of the same numbers.

Requires: pip install mdtraj gemmi matplotlib numpy
"""
import argparse, csv, os, warnings
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("folder", help="folder containing one sub-folder per run")
ap.add_argument("--ligand", help="residue name (e.g. 6NT) or mdtraj selection; "
                                 "default: largest non-protein residue in each run")
ap.add_argument("--align", default="protein and name CA",
                help='atoms to superpose on (default "protein and name CA")')
ap.add_argument("-o", "--out", default="ligand_rmsf.png")
a = ap.parse_args()

SURF, INK, INK2, MUTED, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BAR, RANGE = "#2a78d6", "#0d366b"


# --- read an mmCIF _atom_site table into an mdtraj Topology (mdtraj can't do mmCIF)
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
    """every folder below root that holds both a trajectory and a topology"""
    out = []
    for dirpath, _, files in os.walk(root):
        xtcs = sorted(f for f in files if f.endswith((".xtc", ".dcd", ".trr")))
        tops = ([f for f in files if f in ("topology.cif", "topology.pdb")]
                or sorted(f for f in files if f.endswith((".cif", ".pdb"))))
        if not xtcs or not tops:
            continue
        xtc = "samples.xtc" if "samples.xtc" in xtcs else xtcs[0]
        name = os.path.relpath(dirpath, root)
        out.append((name if name != "." else os.path.basename(os.path.abspath(root)),
                    os.path.join(dirpath, xtc), os.path.join(dirpath, tops[0])))
    return sorted(out)


def ligand_atoms(topology, spec):
    if spec and " " not in spec.strip():
        idx = [x.index for x in topology.atoms if x.residue.name == spec.strip()
               and x.element.symbol != "H"]
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


runs = find_runs(a.folder)
if not runs:
    raise SystemExit(f"No runs found under {a.folder} "
                     "(each run folder needs a .xtc and a topology .cif/.pdb)")
print(f"Found {len(runs)} run folder(s) under {a.folder}\n")

rows = []
for name, xtc, top_path in runs:
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            top = (cif_topology(top_path) if top_path.lower().endswith((".cif", ".mmcif"))
                   else md.load_topology(top_path))
            t = md.load(xtc, top=top)
        lig = ligand_atoms(t.topology, a.ligand)
        if t.n_frames < 2:
            raise ValueError(f"only {t.n_frames} frame(s)")
        sel = t.topology.select(a.align)
        if len(sel) < 3:
            raise ValueError(f"alignment selection matched {len(sel)} atoms")
        t.superpose(t, 0, atom_indices=sel)                     # remove global motion
        L = t.xyz[:, lig] * 10.0                                # Angstrom
        per_atom = np.sqrt(((L - L.mean(0)) ** 2).sum(-1).mean(0))
        rows.append(dict(name=name, rmsf=float(per_atom.mean()),
                         lo=float(per_atom.min()), hi=float(per_atom.max()),
                         n=t.n_frames, lig=t.topology.atom(lig[0]).residue.name,
                         natoms=len(lig)))
        print(f"  {name:40s} RMSF {per_atom.mean():5.2f} A   "
              f"({rows[-1]['lig']}, {len(lig)} atoms, {t.n_frames} samples)")
    except Exception as e:
        print(f"  {name:40s} SKIPPED - {e}")

if not rows:
    raise SystemExit("No run could be analysed.")
rows.sort(key=lambda r: r["rmsf"])                              # lowest RMSF first

# --- plot
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF,
                     "savefig.facecolor": SURF, "font.size": 9, "text.color": INK,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
                     "axes.edgecolor": "#c3c2b7", "axes.spines.top": False,
                     "axes.spines.right": False, "axes.spines.left": False})
many_ligands = len({r["lig"] for r in rows}) > 1


def shorten(name, keep=46):
    return name if len(name) <= keep else name[:keep // 2 - 1] + "…" + name[-keep // 2:]


labels = [shorten(r["name"]) + (f"  ({r['lig']})" if many_ligands else "") for r in rows]
h = 1.9 + 0.34 * len(rows)
lab_in = 0.062 * max(len(l) for l in labels) + 0.25      # room for the names
plot_in, val_in = 6.5, 0.9                               # plot area, value labels
width = lab_in + plot_in + val_in
fig, ax = plt.subplots(figsize=(width, h))
fig.subplots_adjust(left=lab_in / width, right=(lab_in + plot_in + 0.55) / width,
                    top=1 - 1.25 / h, bottom=0.75 / h)
y = np.arange(len(rows))
vals = np.array([r["rmsf"] for r in rows])
ax.barh(y, vals, height=0.62, color=BAR, edgecolor=SURF, linewidth=1.5, zorder=2)
for k, r in enumerate(rows):
    ax.plot([r["lo"], r["hi"]], [k, k], color=RANGE, lw=1.6, zorder=3,
            solid_capstyle="butt", alpha=0.85)
    ax.annotate(f"{r['rmsf']:.2f}", (max(r["hi"], r["rmsf"]), k), xytext=(6, 0),
                textcoords="offset points", va="center", fontsize=8.5, color=INK)
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8.5)
ax.invert_yaxis()
ax.tick_params(axis="y", length=0)
ax.set_xlim(0, max(max(r['hi'] for r in rows), vals.max()) * 1.18)
ax.set_xlabel("Ligand RMSF (Å)   —   lower means the ligand is held more tightly")
ax.grid(axis="x", color=GRID, lw=0.6)
ax.set_axisbelow(True)
x0 = 0.35 / width
fig.suptitle("How much the ligand moves in each run", x=x0, ha="left",
             fontsize=13, fontweight="bold", color=INK, y=1 - 0.3 / h)
fig.text(x0, 1 - 0.80 / h, f"Heavy-atom RMSF of {len(rows)} runs, after superposing every "
         "sample on the protein. Bar = mean over the ligand's atoms,\ndark line = range "
         "from its least to its most mobile atom. Sorted lowest first.",
         color=INK2, fontsize=9, va="center", linespacing=1.5)
fig.savefig(a.out, dpi=200)

csv_path = os.path.splitext(a.out)[0] + ".csv"
with open(csv_path, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["run", "ligand", "ligand_heavy_atoms", "n_samples",
                "rmsf_A", "min_atom_rmsf_A", "max_atom_rmsf_A"])
    for r in rows:
        w.writerow([r["name"], r["lig"], r["natoms"], r["n"],
                    f"{r['rmsf']:.3f}", f"{r['lo']:.3f}", f"{r['hi']:.3f}"])
print(f"\nSaved {a.out} and {csv_path}")
