#!/usr/bin/env python3
"""
Figure 3 - cell-type-resolved differential expression. Plotting script.

Reads the source CSVs in the Figure_3 folder and reproduces all four installed
panels. Those CSVs come from compute_figure3_tables.py, which runs the permutation
battery and GO enrichment against the per-cell-class differential-expression tables
under E18p5_clean/results/tables/. This script does only the plotting, so the figure
can be regenerated without re-running the battery.

Input CSVs (in FIG3_DIR/tables/, default "../tables"):
  Table_S4_GO_terms.csv          panel A: curated GO terms with fold enrichment and q
  Table_S3_gene_set_battery.csv  panels B-D: 10 gene sets x 9 classes, mean log2FC,
                                 z vs random, nominal perm P, and BH q across the
                                 original 108 tests. DIRECT GO annotation, not the
                                 GOALL closure - see the CSV header.

Output: panel_A_GO_merged.pdf, panel_B_guidance_adhesion.pdf,
        panel_C_synaptic_function.pdf, panel_D_proliferation_stress.pdf

Usage: python plot_figure3.py [FIG3_DIR] [OUT_DIR]
"""
import os, sys, csv
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

_args = [a for a in sys.argv[1:] if not a.startswith("-")]
FIG3 = _args[0] if len(_args) > 0 else ".."
OUT  = _args[1] if len(_args) > 1 else "."

mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})
UPC, DNC = "#C0453B", "#3A6EA5"   # up / down in mutant
CLASS_ORDER = ["Cycling RG", "Immature astrocytes", "Deep-layer EN", "Upper-layer EN",
               "SPN-D1", "SPN-D2", "LGE-IN precursors", "MGE interneurons", "CGE interneurons"]
# panel B row-category bands (left) and column-group headers (top)

def read_csv(p):
    """Read a shipped source CSV, skipping leading '#' provenance comments."""
    with open(p) as f:
        lines = [ln for ln in f.read().split("\n")
                 if ln and not ln.lstrip('"').startswith("#")]
    return list(csv.DictReader(lines))

# ---------- 3A: curated GO enrichment, mirrored ----------
GOTERMS = os.path.join(FIG3, "tables", "Table_S4_GO_terms.csv")


def panel_A():
    """Panel A: GO enrichment of the global DE result, redundant terms collapsed.

    Mirrored axis: DOWN-in-mutant terms extend left from zero, UP-in-mutant terms
    extend right, so direction reads without a colour key. Fold enrichment is
    plotted as magnitude on both sides; the left axis is labelled with positive
    values because it is an enrichment, not a negative quantity.
    """
    rows = [r for r in read_csv(GOTERMS) if r["figure"] == "Figure 3A"]
    down = [r for r in rows if r["direction"].upper().startswith("D")]
    up   = [r for r in rows if r["direction"].upper().startswith("U")]
    down = sorted(down, key=lambda r: float(r["fold_enrichment"]))
    up   = sorted(up,   key=lambda r: float(r["fold_enrichment"]))
    ordered = down + up

    fig, ax = plt.subplots(figsize=(5.0, 2.6))
    for k, r in enumerate(ordered):
        fe = float(r["fold_enrichment"])
        isdown = r["direction"].upper().startswith("D")
        x = -fe if isdown else fe
        col = DNC if isdown else UPC
        ax.plot([0, x], [k, k], color=col, lw=0.9, solid_capstyle="butt", zorder=2)
        ax.plot(x, k, "o", ms=3.6, color=col, zorder=3)
    ax.axvline(0, color="0.3", lw=0.8, zorder=1)
    ax.axhline(len(down) - 0.5, color="0.85", lw=0.6, zorder=1)
    ax.set_yticks(range(len(ordered)))
    ax.set_yticklabels([r["term_shown"] for r in ordered], fontsize=6.2, color="black")
    lim = max(float(r["fold_enrichment"]) for r in ordered) * 1.14
    ax.set_xlim(-lim, lim)
    ticks = [t for t in range(0, int(lim) + 1, 2)]
    ax.set_xticks([-t for t in ticks[::-1]] + ticks[1:])
    ax.set_xticklabels([str(t) for t in ticks[::-1]] + [str(t) for t in ticks[1:]],
                       fontsize=6.2)
    ax.set_xlabel("Fold enrichment", fontsize=7)
    ax.text(-lim * 0.97, len(ordered) - 0.35, "LOWER in mutant", fontsize=6.2,
            color=DNC, ha="left", va="bottom")
    ax.text(lim * 0.97, len(ordered) - 0.35, "HIGHER in mutant", fontsize=6.2,
            color=UPC, ha="right", va="bottom")
    ax.set_ylim(-0.7, len(ordered) + 0.1)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.savefig(os.path.join(OUT, "panel_A_GO_merged.pdf"), bbox_inches="tight")
    plt.close(fig)


BATTERY = os.path.join(FIG3, "tables", "Table_S3_gene_set_battery.csv")

PANEL_FILE = {"B": "panel_B_guidance_adhesion.pdf",
              "C": "panel_C_synaptic_function.pdf",
              "D": "panel_D_proliferation_stress.pdf"}
PANEL_HEIGHT = {"B": 1.9, "C": 1.9, "D": 2.8}
LONG_NAME = {"Ionotropic glutamate receptor sig": "Ionotropic glutamate receptor signalling"}


def panels_BCD():
    """Panels B, C and D: gene-set shifts within each cell class.

    Panels B and C carry a small white dot on cells reaching the NOMINAL
    permutation threshold P < 0.05. Panel D carries no markers: it has one
    nominal cell against 1.8 expected by chance across its 36 tests, fewer
    than chance, so marking it would imply a signal the data do not support.
    No cell in the battery survives FDR correction (lowest q = 0.076); the
    dots mark where the panel A enrichment localises, not independent claims.
    """
    rows = [r for r in read_csv(BATTERY) if not r["panel"].startswith("#")]
    classes = list(dict.fromkeys(r["cell_type"] for r in rows))
    ci = {c: k for k, c in enumerate(classes)}
    groups = [("Dorsal", 0, 3), ("LGE-derived", 4, 6), ("Cortical IN", 7, 8)]

    for pid, figh, mark in [("B", 2.1, True), ("C", 2.1, True), ("D", 2.4, False)]:
        sets_ = list(dict.fromkeys(r["gene_set"] for r in rows if r["panel"] == pid))
        gi = {s: k for k, s in enumerate(sets_)}
        Z = np.full((len(sets_), len(classes)), np.nan)
        P = np.full_like(Z, np.nan)
        for r in rows:
            if r["panel"] != pid:
                continue
            Z[gi[r["gene_set"]], ci[r["cell_type"]]] = float(r["z_vs_random"])
            P[gi[r["gene_set"]], ci[r["cell_type"]]] = float(r["perm_p"])

        fig, ax = plt.subplots(figsize=(5.4, figh))
        # fixed colour scale, shared by all three panels so they are comparable
        im = ax.imshow(Z, cmap="RdBu_r", vmin=-4.3, vmax=4.3, aspect="auto")
        if mark:
            for a in range(Z.shape[0]):
                for b in range(Z.shape[1]):
                    if P[a, b] < 0.05:
                        ax.plot(b, a, "o", ms=2.2, mfc="white", mec="none", zorder=3)
        ax.set_xticks(range(len(classes)))
        ax.set_xticklabels(classes, rotation=40, ha="right", fontsize=6.0)
        ax.set_yticks(range(len(sets_)))
        ax.set_yticklabels(sets_, fontsize=6.6, color="black")
        if pid in ("B", "C"):
            for nm, a, b in groups:
                ax.plot([a - 0.42, b + 0.42], [-0.72, -0.72], color="0.35",
                        lw=1.0, clip_on=False)
                ax.text((a + b) / 2, -0.95, nm, ha="center", va="bottom",
                        fontsize=6.0, color="0.3", clip_on=False)
        cb = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.02)
        cb.set_label("Shift vs random gene sets (z)", fontsize=6.6)
        cb.ax.tick_params(labelsize=6)
        fig.savefig(os.path.join(OUT, PANEL_FILE[pid]), bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    panel_A(); panels_BCD()
    print("Figure 3 panels A-D written to", os.path.abspath(OUT))
