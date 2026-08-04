#!/usr/bin/env python3
"""
Figure 2F-H: D2 spiny projection neuron subclustering.

F  the two states on the UMAP, restricted to the D2 region
G  abundance of each state per genotype, against the analysed scope
H  transferred pseudotime within each state, by genotype

Reads the tables written by 12_d2_subclustering.R. Also writes the Supplementary
Figure S2F marker heatmap, since it comes from the same source table.

Usage:  python plot_figure2FGH.py <table_dir> [out_dir]
"""
import csv, os, sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

TBL = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
os.makedirs(OUT, exist_ok=True)

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.labelsize": 7.5,
})

SC    = ["precursor-like", "differentiating"]
CS    = {"precursor-like": "#C0453B", "differentiating": "#5B7FA6"}
WT_C  = "0.72"
LAB   = {"precursor-like": "precursor-\nlike", "differentiating": "differ-\nentiating"}
SCOPE = {"WT": 8444, "MUT": 10852}


def read_csv(path):
    """Read a CSV whose leading '#' lines are provenance notes."""
    lines = [l for l in open(path).read().split("\n") if l and not l.lstrip('"').startswith("#")]
    return list(csv.DictReader(lines))


cells = read_csv(os.path.join(TBL, "d2_subcluster_cells.csv"))
bg    = read_csv(os.path.join(TBL, "d2_umap_background.csv"))
mk    = read_csv(os.path.join(TBL, "d2_subcluster_markers.csv"))
stats = read_csv(os.path.join(TBL, "Table_S12_D2_subclustering.csv"))
ST    = {r["subcluster"]: r for r in stats}

pt_of = lambda s, g: np.array([float(r["pt"]) for r in cells
                               if r["subcluster"] == s and r["geno"] == g
                               and r["pt"] not in ("", "NA")])
xy_of = lambda s: (np.array([float(r["UMAP1"]) for r in cells if r["subcluster"] == s]),
                   np.array([float(r["UMAP2"]) for r in cells if r["subcluster"] == s]))

# Zoom window: 1st-99th percentile of D2 cells plus a small margin, so the two
# states separate visually instead of overlapping in the full embedding.
ax_ = np.array([float(r["UMAP1"]) for r in cells])
ay_ = np.array([float(r["UMAP2"]) for r in cells])
PAD = 0.35
ZX = (np.percentile(ax_, 1) - PAD, np.percentile(ax_, 99) + PAD)
ZY = (np.percentile(ay_, 1) - PAD, np.percentile(ay_, 99) + PAD)


def fmt_p(p):
    """One-significant-figure scientific notation, so the label tracks the table."""
    if p >= 0.05:
        return "n.s."
    exp = int(np.floor(np.log10(p)))
    mant = p / 10 ** exp
    return f"$P$ = {mant:.0f} \u00d7 10$^{{{exp}}}$"


def panel_F():
    fig, ax = plt.subplots(figsize=(2.9, 2.8))
    bx = [(float(r["UMAP1"]), float(r["UMAP2"])) for r in bg
          if ZX[0] < float(r["UMAP1"]) < ZX[1] and ZY[0] < float(r["UMAP2"]) < ZY[1]]
    ax.scatter([p[0] for p in bx], [p[1] for p in bx], s=4.5, c="0.90", lw=0, rasterized=True)
    for s in ["differentiating", "precursor-like"]:
        x, y = xy_of(s)
        ax.scatter(x, y, s=8.0, c=CS[s], lw=0, alpha=0.85, rasterized=True)
    px, py = xy_of("precursor-like")
    dx, dy = xy_of("differentiating")
    ax.text(np.percentile(px, 92), np.percentile(py, 94) + 0.20, "precursor-like",
            fontsize=7, color="white", ha="center", va="bottom",
            bbox=dict(boxstyle="round,pad=0.24", fc=CS["precursor-like"], ec="none"))
    ax.text(np.percentile(dx, 30), np.percentile(dy, 4) - 0.10, "differentiating",
            fontsize=7, color="white", ha="center", va="top",
            bbox=dict(boxstyle="round,pad=0.24", fc=CS["differentiating"], ec="none"))
    ax.set_xlim(*ZX); ax.set_ylim(ZY[0] - 0.10, ZY[1])
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlabel("UMAP 1 (zoom on D2)", fontsize=7, labelpad=2)
    ax.set_ylabel("UMAP 2", fontsize=7, labelpad=1)
    ax.set_title("D2 SPNs contain two states", fontsize=7.5, loc="left", pad=6)
    fig.savefig(os.path.join(OUT, "panel_F_D2_subclusters_umap.pdf"), bbox_inches="tight")
    plt.close(fig)


def panel_G():
    fig, ax = plt.subplots(figsize=(2.7, 2.8))
    wd = 0.34
    for k, s in enumerate(SC):
        r = ST[s]
        pw, pm = float(r["pct_wt"]), float(r["pct_mut"])
        ratio, p = float(r["abundance_ratio"]), float(r["abundance_p"])
        ax.bar(k - wd / 2, pw, wd, color=WT_C, lw=0)
        ax.bar(k + wd / 2, pm, wd, color=CS[s], lw=0)
        for xo, lb in [(-wd / 2, "wt"), (wd / 2, "mut")]:
            ax.text(k + xo, -0.18, lb, ha="center", va="top", fontsize=6, color="0.4")
        top = max(pw, pm)
        ax.text(k, top + 0.16, f"{ratio:.1f}\u00d7", ha="center", fontsize=7,
                color=CS[s] if p < 0.05 else "0.45",
                fontweight="bold" if p < 0.05 else "normal")
        ax.text(k, top + 0.55, fmt_p(p), ha="center", fontsize=6, color="0.45")
    ax.set_xticks([0, 1]); ax.set_xticklabels([LAB[s] for s in SC], fontsize=7)
    ax.tick_params(axis="x", pad=10)
    ax.set_ylabel("% of analysed nuclei", fontsize=7.5); ax.set_ylim(0, 7.4)
    ax.set_title("Only the precursor pool expands", fontsize=7.5, loc="left", pad=6)
    fig.savefig(os.path.join(OUT, "panel_G_D2_subcluster_abundance.pdf"), bbox_inches="tight")
    plt.close(fig)


def panel_H():
    fig, ax = plt.subplots(figsize=(3.3, 2.8))
    for k, s in enumerate(SC):
        for j, g in enumerate(["WT", "MUT"]):
            x = k * 1.25 + j * 0.42
            v = pt_of(s, g)
            parts = ax.violinplot([v], positions=[x], widths=0.38, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(CS[s] if g == "MUT" else WT_C)
                b.set_alpha(0.8); b.set_edgecolor("none")
            ax.hlines(np.median(v), x - 0.16, x + 0.16, color="black", lw=1.0, zorder=3)
            ax.text(x, 31.4, "wt" if g == "WT" else "mut",
                    ha="center", va="bottom", fontsize=6, color="0.4")
        p = float(ST[s]["pt_p"])
        ytop = max(np.percentile(pt_of(s, g), 75) for g in ("WT", "MUT")) + 4.2
        ax.plot([k * 1.25, k * 1.25, k * 1.25 + 0.42, k * 1.25 + 0.42],
                [ytop, ytop + 1.1, ytop + 1.1, ytop], lw=0.6, c="0.4")
        ax.text(k * 1.25 + 0.21, ytop + 2.0, fmt_p(p),
                ha="center", fontsize=6, color="0.35" if p >= 0.05 else "black")
    ax.set_xticks([0.21, 1.46]); ax.set_xticklabels([LAB[s] for s in SC], fontsize=7)
    ax.tick_params(axis="x", pad=10)
    ax.set_ylabel("Transferred pseudotime", fontsize=7.5)
    ax.set_xlim(-0.35, 2.05); ax.set_ylim(30, 70)
    ax.set_title("Differentiating cells are delayed in mutant", fontsize=7.5, loc="left", pad=6)
    fig.savefig(os.path.join(OUT, "panel_H_D2_within_state_maturation.pdf"), bbox_inches="tight")
    plt.close(fig)


def panel_S2F():
    """Supplementary S2F: marker detection per state, SPN block above LGE block."""
    fig, ax = plt.subplots(figsize=(2.3, 3.0))
    spn = [r for r in mk if r["group"] == "SPN differentiation"]
    lge = [r for r in mk if r["group"] == "LGE precursor"]
    rows = (lge + spn)[::-1]
    M = np.array([[float(r["prec"]), float(r["diff"])] for r in rows])
    ax.imshow(M, cmap="Reds", vmin=0, vmax=100, aspect="auto")
    for i in range(len(rows)):
        for j in range(2):
            ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=6,
                    color="white" if M[i, j] > 55 else "0.15")
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels([f"$\\it{{{r['gene']}}}$" for r in rows], fontsize=6.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels([LAB[s] for s in SC], fontsize=6.5)
    ns = len(spn)
    ax.axhline(ns - 0.5, color="0.3", lw=0.8)
    ax.text(1.62, (ns - 1) / 2, "SPN\ndifferentiation", fontsize=6, va="center",
            ha="left", color="0.3")
    ax.text(1.62, ns + (len(rows) - ns - 1) / 2, "LGE\nprecursor", fontsize=6,
            va="center", ha="left", color="0.3")
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("% cells expressing", fontsize=7, loc="left", pad=6)
    fig.savefig(os.path.join(OUT, "panel_F_D2_subcluster_markers.pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    panel_F(); panel_G(); panel_H(); panel_S2F()
    print("wrote Figure 2F-H and Supplementary S2F to", os.path.abspath(OUT))
