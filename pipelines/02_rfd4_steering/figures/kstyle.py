"""Shared style for the campaign figures (colours fixed per arm).

Inputs: $KEMP_RESULTS/tables/*.csv (written by kemp.export_results). Outputs: $KEMP_RESULTS/figures/YYYYMMDD_fNN_<name>.png
"""
import os, textwrap, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES = os.environ["KEMP_RESULTS"]; T = os.path.join(RES, "tables"); OUT = os.path.join(RES, "figures"); DATE = os.environ.get("FIG_DATE", "20261122")
os.makedirs(OUT, exist_ok=True)
INK, NOTE, SPINE = "#1f2933", "#52606d", "#52606d"
ARMS = ["unsteered", "steer_tmol", "steer_tmol_ch18"]
COL = {"unsteered": "#8a9099", "steer_tmol": "#3f7fe0", "steer_tmol_ch18": "#e2661a"}            # baseline grey, steered blue, steered + channel orange (fixed in every figure)
LAB = {"unsteered": "unsteered", "steer_tmol": "steered", "steer_tmol_ch18": "steered + channel 1.8 A"}
MOT = {"hgnt": "o", "hgwt": "s"}; MOTLAB = {"hgnt": "hgnt (Gln50, Asp127)", "hgwt": "hgwt (+ Trp44)"}
MM = 25.4
plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6.5, "font.family": "DejaVu Sans", "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK, "axes.edgecolor": SPINE, "axes.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "pdf.fonttype": 42, "savefig.dpi": 300})
def fig(w_mm=183, h_mm=90, **kw): return plt.subplots(figsize=(w_mm / MM, h_mm / MM), **kw)
def footnote(f, text, y=0.005): f.text(0.01, y, text, fontsize=6, color=NOTE, ha="left", va="bottom")
def save(f, n, name): 
    for ext in ("png",): f.savefig(os.path.join(OUT, f"{DATE}_f{n:02d}_{name}.{ext}"), bbox_inches="tight", pad_inches=0.04)
    plt.close(f); print("saved", f"{DATE}_f{n:02d}_{name}")
def ev(include_sheet=False):
    e = pd.read_csv(os.path.join(T, "evaluated_sequences_full_suite.csv"))
    if not include_sheet: e = e[~e.arm.str.contains("sheet")]
    return e
def designs(include_sheet=False):
    d = pd.read_csv(os.path.join(T, "all_designs_metrics.csv")); d["batch"] = ["b%02d" % int(b) for b in d.batch]
    if not include_sheet: d = d[~d.arm.str.contains("sheet")]
    return d
def cands(include_sheet=False):
    c = pd.read_csv(os.path.join(T, "candidates.csv"))
    if not include_sheet: c = c[~c.arm.str.contains("sheet")]
    return c
def legend_mono(ax, handles, labels, title=None, loc="center left", anchor=(1.01, 0.5)):
    return ax.legend(handles, labels, prop={"family": "monospace", "size": 6.5}, frameon=False, loc=loc, bbox_to_anchor=anchor, title=title, title_fontsize=6.5, handlelength=1.6)
def arm_n(df, col="seq_id"): return {a: int((df.arm == a).sum()) for a in ARMS}
def hist_panels(axs, df, col, label, bins=None, motifs=("hgnt", "hgwt"), zero=False, xlim=None):
    """overlaid step histograms (density) per arm; one axis per motif."""
    vals = df[col].dropna()
    if bins is None: bins = np.linspace(vals.min(), vals.max(), 14) if xlim is None else np.linspace(xlim[0], xlim[1], 14)
    for ax, m in zip(np.atleast_1d(axs), motifs):
        for a in ARMS:
            v = df[(df.motif == m) & (df.arm == a)][col].dropna()
            if len(v) < 2: continue
            ax.hist(v, bins=bins, histtype="stepfilled", color=COL[a], alpha=0.16, lw=0, density=True)
            ax.hist(v, bins=bins, histtype="step", color=COL[a], lw=1.3 if a != "unsteered" else 1.0, ls="-" if a != "unsteered" else "--", density=True)
            ax.plot([v.median()], [0], marker="v", ms=4, color=COL[a], mec="white", mew=0.5, clip_on=False, zorder=5, transform=ax.get_xaxis_transform() if False else ax.transData)
        if zero: ax.axvline(0, color=INK, lw=0.6, ls=":")
        ax.set_xlabel(textwrap.fill(label, 17), fontsize=6.5); ax.set_yticks([])
        if xlim: ax.set_xlim(*xlim)
def arm_legend(ax, df, extra=""):
    hs = [plt.Line2D([], [], color=COL[a], lw=1.4, ls="--" if a == "unsteered" else "-") for a in ARMS]
    nm = lambda a: ", ".join(f"{m}:{int(((df.arm == a) & (df.motif == m)).sum())}" for m in ("hgnt", "hgwt"))
    return legend_mono(ax, hs, [f"{LAB[a]:<22} n({nm(a)})" for a in ARMS], title="evaluated sequences" + extra)
def scat(ax, df, x, y, size=16, jitter=0.0, ring=None, seed=1):
    rng = np.random.default_rng(seed)
    for a in ARMS:
        for m in ("hgnt", "hgwt"):
            d = df[(df.arm == a) & (df.motif == m)]
            if len(d) == 0: continue
            jx = rng.normal(0, jitter, len(d)) if jitter else 0; jy = rng.normal(0, jitter, len(d)) if jitter else 0
            ax.scatter(d[x] + jx, d[y] + jy, s=size, marker=MOT[m], facecolor="white" if a == "unsteered" else COL[a], edgecolor=COL[a] if a == "unsteered" else "white", linewidth=0.9 if a == "unsteered" else 0.4, zorder=3, alpha=0.95)
    if ring is not None and ring.any():
        r = df[ring]; ax.scatter(r[x], r[y], s=size * 3.2, marker="o", facecolor="none", edgecolor=INK, linewidth=0.8, zorder=4)
def scat_legend(ax, df, ring_label=None):
    hs = [plt.Line2D([], [], marker="o", ls="", ms=5, mfc="white" if a == "unsteered" else COL[a], mec=COL[a] if a == "unsteered" else "white") for a in ARMS]
    hs += [plt.Line2D([], [], marker=MOT[m], ls="", ms=5, mfc=NOTE, mec="white") for m in MOT]
    lab = [f"{LAB[a]:<22} n={int((df.arm == a).sum())}" for a in ARMS] + [MOTLAB[m] for m in MOT]
    if ring_label: hs.append(plt.Line2D([], [], marker="o", ls="", ms=8, mfc="none", mec=INK)); lab.append(ring_label)
    return legend_mono(ax, hs, lab)
