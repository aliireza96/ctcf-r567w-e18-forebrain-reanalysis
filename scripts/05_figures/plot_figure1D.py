#!/usr/bin/env python3
"""
Figure 1D - cell-type composition, mutant versus wild type.

No significance markers are drawn inside the panel. A point is FILLED when it passes
BOTH the FDR test (BH q < 0.05 across the 12 cell types) and the equal-N resampling
check (|delta| >= 1 percentage point and sign consistency >= 0.95); open otherwise.
The figure legend states the criteria and the per-cell-type q values are in
Table_S1_composition_FDR.csv.

Reads tables/Table_S1_composition_FDR.csv and reproduces the installed panel.
That table is produced by compute_figure1D_stats.py from the per-cell-type nucleus
counts in E18p5_clean/results/tables/04e_composition_descriptive_celltype_figure.csv.

Usage: python plot_figure1D.py [FIG1_DIR] [OUT_DIR]
"""
import os, sys, csv
import matplotlib as mpl
import matplotlib.pyplot as plt

_args = [a for a in sys.argv[1:] if not a.startswith("-")]
FIG1 = _args[0] if len(_args) > 0 else ".."
OUT  = _args[1] if len(_args) > 1 else "."

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})
UP, DN = "#C0453B", "#3A6EA5"
NICE = {"Deep layer EN": "Deep-layer EN", "Upper layer EN": "Upper-layer EN",
        "Immature Astrocytes": "Immature astrocytes", "MGE-IN": "MGE interneurons",
        "Migrating CGE-derived IN": "CGE interneurons", "LGE-IN prec": "LGE-IN precursors",
        "ExtendedAmygdala_GABA": "Extended amygdala GABA", "Cajal-Retzius": "Cajal-Retzius",
        "OPC": "OPC", "SPN-D1": "SPN-D1", "SPN-D2": "SPN-D2", "Cycling RG": "Cycling RG"}


def read_csv(p):
    with open(p) as f:
        lines = [ln for ln in f.read().split("\n")
                 if ln and not ln.lstrip('"').startswith("#")]
    return list(csv.DictReader(lines))


def main():
    rows = read_csv(os.path.join(FIG1, "tables", "Table_S1_composition_FDR.csv"))
    rows = sorted(rows, key=lambda r: float(r["delta_pct_points"]))

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    for k, r in enumerate(rows):
        d = float(r["delta_pct_points"])
        q = float(r["q_value_BH"])
        lo, hi = float(r["ci_low_pct"]), float(r["ci_high_pct"])
        col = UP if d > 0 else DN
        # pale bar: equal-N resampling interval (sampling stability, NOT a genotype CI)
        ax.plot([lo, hi], [k, k], color=col, lw=0.9, solid_capstyle="butt", alpha=0.5, zorder=2)
        ax.plot([0, d], [k, k], color=col, lw=1.5, solid_capstyle="butt", zorder=3)
        # FILLED requires BOTH criteria: FDR q<0.05 AND the equal-N resampling check
        # (|delta| >= 1 pp and sign consistency >= 0.95). Two cell types pass FDR but
        # fail resampling (Cajal-Retzius, SPN-D2) and are drawn open with their stars.
        filled = q < 0.05 and r["passes_resampling"] == "yes"
        ax.plot(d, k, "o", ms=5.2, mfc=col if filled else "white", mec=col, mew=1.1, zorder=4)
    ax.axvline(0, color="0.55", lw=0.8, ls="--", zorder=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([NICE.get(r["celltype"], r["celltype"]) for r in rows],
                       fontsize=6.8, color="black")
    ax.set_xlabel("Difference in percentage of nuclei (mutant \u2212 wild type)", fontsize=7.2)
    ax.set_xlim(-8.4, 7.0)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.set_xticks([-8, -6, -4, -2, 0, 2, 4, 6])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.savefig(os.path.join(OUT, "panel_D_composition_FDR.pdf"), bbox_inches="tight")
    plt.close(fig)
    print("Figure 1D written to", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
