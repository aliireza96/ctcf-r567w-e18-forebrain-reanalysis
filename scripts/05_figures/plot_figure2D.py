#!/usr/bin/env python3
"""
Figure 2D - shift in median transferred pseudotime, mutant versus wild type.

Reads tables/Table_S2_pseudotime_shift_FDR.csv and reproduces the installed panel.

Points are FILLED when the Benjamini-Hochberg q across the nine classes is below
0.05, open otherwise. No significance markers are drawn inside the panel; the
legend states the criterion. The horizontal bar is the equal-N resampling interval
carried over from the original analysis, which bounds sampling-depth instability
and is not a confidence interval for a genotype effect.

Usage: python plot_figure2D.py [FIG2_DIR] [OUT_DIR]
"""
import os, sys, csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

_args = [a for a in sys.argv[1:] if not a.startswith("-")]
FIG2 = _args[0] if len(_args) > 0 else ".."
OUT  = _args[1] if len(_args) > 1 else "."

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})
GCOL = {"RG root (dorsal)": "#1B7A6B", "LGE root": "#C0453B",
        "MGE root": "#5B4FA8", "CGE root": "#C08A2E"}


def read_csv(p):
    with open(p) as f:
        lines = [ln for ln in f.read().split("\n")
                 if ln and not ln.lstrip('"').startswith("#")]
    return list(csv.DictReader(lines))


def main():
    rows = read_csv(os.path.join(FIG2, "tables", "Table_S2_pseudotime_shift_FDR.csv"))
    n = len(rows)

    fig, ax = plt.subplots(figsize=(6.6, 3.3))
    for k, r in enumerate(rows):
        y = n - 1 - k
        col = GCOL[r["root"]]
        if r["ci_low"]:
            ax.plot([float(r["ci_low"]), float(r["ci_high"])], [y, y],
                    color=col, lw=1.4, solid_capstyle="butt", zorder=2)
        filled = float(r["q_value_BH"]) < 0.05
        ax.plot(float(r["delta_median_pseudotime"]), y, "o", ms=6,
                mfc=col if filled else "white", mec=col, mew=1.3, zorder=4)
        ax.text(2.45, y, f"{r['n_wt']} / {r['n_mut']}", fontsize=5.8,
                va="center", ha="right", color="0.45")
    ax.axvline(0, color="0.55", lw=0.8, ls="--", zorder=1)

    prev = None
    for k, r in enumerate(rows):
        if r["root"] != prev and k > 0:
            ax.axhline(n - 0.5 - k, color="0.88", lw=0.6, zorder=0)
        prev = r["root"]

    ax.set_yticks(range(n))
    ax.set_yticklabels([r["lineage"] for r in rows][::-1], fontsize=6.8, color="black")
    ax.set_xlabel("\u0394 median pseudotime (mutant \u2212 wild type)", fontsize=7.2)
    ax.set_xlim(-4.5, 2.6)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xticks([-4, -3, -2, -1, 0, 1, 2])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.text(2.45, n - 0.45, "n WT / MUT", fontsize=5.8, ha="right",
            va="center", color="0.45")

    # root bars in axes coordinates, so they clear the cell-class labels
    fig.canvas.draw()
    inv = ax.transAxes.inverted()
    fy = lambda v: inv.transform(ax.transData.transform((0, v)))[1]
    for g in dict.fromkeys(r["root"] for r in rows):
        idx = [n - 1 - k for k, r in enumerate(rows) if r["root"] == g]
        ax.plot([-0.30, -0.30], [fy(min(idx) - 0.42), fy(max(idx) + 0.42)],
                color=GCOL[g], lw=2.2, transform=ax.transAxes, clip_on=False, zorder=5)
        ax.text(-0.315, fy(float(np.mean(idx))), g, fontsize=6.1, color=GCOL[g],
                ha="right", va="center", transform=ax.transAxes, clip_on=False)

    fig.savefig(os.path.join(OUT, "panel_D_delta_median_FDR.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("Figure 2D written to", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
