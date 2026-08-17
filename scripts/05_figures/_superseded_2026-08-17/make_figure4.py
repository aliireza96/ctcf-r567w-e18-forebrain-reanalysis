#!/usr/bin/env python3
"""
make_figure4.py

Regenerates the four installed panels of Figure 4 from the shipped source CSVs.

  panel_A_global_occupancy_loss.pdf   CTCF signal change by MACS2 peak class
  panel_B_umotif_gradient.pdf         U-motif content vs occupancy-change decile
  panel_C_lost_site_location.pdf      genomic-position composition per decile
  panel_D_cpcdh_locus.pdf             cPcdh locus, three assays on one axis

Reads only from the CSVs in the parent Figure_4 folder plus, for panel A's
per-peak distributions, figure4_normalised_peak_signal.csv. Writes nothing back
to those files.

Usage:  python make_figure4.py [FIGURE4_DIR] [OUT_DIR]

Requires: matplotlib, numpy.  No seaborn, no project code.

PROVENANCE OF THE GENOMIC UNIT AND SIGNAL VALUES
------------------------------------------------
Panels A-C use 45,021 disjoint components obtained after collapsing exact-
coordinate MACS2 multi-summit records and reducing all overlapping WT and
mutant intervals. They do not treat the 47,052 source records as independent
sites.

All signal in this figure comes from the GEO spike-in NORMALISED bigwigs
(GSM6614265 wild type, GSM6614266 homozygous), NOT from the narrowPeak
signalValue column. The paper applies its HEK293T spike-in scale factor when
building bigwigs; peak calling is run with default MACS2 on the mouse reads
after the mixed-genome split, so narrowPeak signalValue is raw per-library
fold-enrichment. Using it would confound occupancy loss with library depth.
"""
import csv
import os
import sys

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

FIG4 = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..")
OUT = sys.argv[2] if len(sys.argv) > 2 else FIG4

mpl.rcParams.update({
    "font.size": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
})

# colour bindings, consistent across panels
CLS_COL = {"MAINTAINED": "#3E7CB1", "LOST": "#C0453B", "GAINED": "#E1A73E"}
POS_COL = {"Promoter": "#2C6E9B", "Intron": "#7A9A5B", "Distal": "#B03A2E", "Other": "0.62"}
REG_COL = {"Pcdha": "#7A9A5B", "Pcdhb": "#B03A2E", "Pcdhg": "#5B7FA6"}


def read_csv_skip_comments(path):
    """Read a CSV whose leading '#' lines are provenance comments."""
    with open(path) as fh:
        lines = [ln for ln in fh.read().split("\n")
                 if ln and not ln.lstrip('"').startswith("#")]
    return list(csv.DictReader(lines))


# ---------------------------------------------------------------- panel A
def panel_A():
    """Per-peak signal change, grouped by MACS2 peak class.

    LOST is partly definitional (a component is WT-only largely because signal
    fell). MAINTAINED is the informative group: nothing in "called in both"
    constrains the signal change, yet its median is -0.61. GAINED changes only
    modestly (+0.11 median).
    """
    rows = read_csv_skip_comments(os.path.join(FIG4, "tables/Table_S7_CTCF_peak_signal.csv"))
    lf = {}
    for r in rows:
        lf.setdefault(r["peak_class"], []).append(float(r["log2FC_norm"]))
    lf = {k: np.asarray(v) for k, v in lf.items()}

    order = ["LOST", "MAINTAINED", "GAINED"]
    labels = [f"Lost\ncalled in WT only\n(n={len(lf['LOST']):,})",
              f"Maintained\ncalled in both\n(n={len(lf['MAINTAINED']):,})",
              f"Gained\ncalled in mutant only\n(n={len(lf['GAINED']):,})"]

    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    parts = ax.violinplot([lf[c] for c in order], positions=range(3),
                          widths=0.8, showextrema=False, showmedians=False)
    for body, cls in zip(parts["bodies"], order):
        body.set_facecolor(CLS_COL[cls]); body.set_alpha(0.72); body.set_edgecolor("none")
    for i, cls in enumerate(order):
        med = float(np.median(lf[cls]))
        ax.plot([i - 0.3, i + 0.3], [med, med], color="black", lw=1.5, zorder=6)
        ax.text(i, med - 0.40, f"{med:+.2f}", ha="center", va="top",
                fontsize=6.8, weight="bold")
    ax.axhline(0, color="0.25", lw=1.0, ls="--")
    ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=6.3)
    ax.set_ylabel("CTCF ChIP signal change in mutant\n"
                  "log$_2$(mutant / wild type), spike-in normalised", fontsize=7.2)
    ax.set_ylim(-4.8, 3.4)
    ax.text(2.44, 0.16, "no change", fontsize=5.9, color="0.4", va="bottom", ha="right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(os.path.join(OUT, "panel_A_global_occupancy_loss.pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- panel B
def panel_B():
    """U-motif content against occupancy-change decile.

    The x-axis is a RANK, deliberately. Plotting signed log2FC puts "more binding"
    on the right and runs the curve right-to-left; negating it to "occupancy lost"
    makes rightward mean less binding. Both force the reader to decode a direction.
    """
    d = read_csv_skip_comments(os.path.join(FIG4, "tables/Table_S8_umotif_by_loss_decile.csv"))
    rank = np.array([int(r["decile_rank"]) for r in d])
    pct = np.array([float(r["pct_proxy_positive"]) for r in d])
    lo = np.array([float(r["ci_lo"]) for r in d])
    hi = np.array([float(r["ci_hi"]) for r in d])
    n_tot = sum(int(r["n_components"]) for r in d)
    n_pos = sum(int(r["n_proxy_positive"]) for r in d)
    overall = 100.0 * n_pos / n_tot

    fig, ax = plt.subplots(figsize=(4.5, 2.9))
    ax.fill_between(rank, lo, hi, color="#B03A2E", alpha=0.18, lw=0)
    ax.plot(rank, pct, "-", color="#B03A2E", lw=1.4, zorder=3)
    ax.plot(rank, pct, "o", color="#B03A2E", ms=4.4, mec="white", mew=0.7, zorder=4)
    ax.axhline(overall, color="0.45", lw=0.9, ls=":", zorder=1)
    ax.text(10.35, overall + 0.5, f"all sites, {overall:.1f}%",
            fontsize=5.9, color="0.42", ha="right", va="bottom")
    ax.annotate(f"{pct[0]:.1f}%", xy=(1, pct[0]), xytext=(1.15, pct[0] + 1.5),
                fontsize=6.8, weight="bold", color="#B03A2E")
    ax.annotate(f"{pct[-1]:.1f}%", xy=(10, pct[-1]), xytext=(9.85, pct[-1] - 2.4),
                fontsize=6.8, weight="bold", color="#B03A2E", ha="right")
    ax.set_xticks(rank); ax.set_xticklabels([str(i) for i in rank], fontsize=6.4)
    ax.set_xlabel("CTCF union components ranked by occupancy change, decile\n"
                  "(1 = most binding lost   \u2192   10 = binding gained)", fontsize=7.2)
    ax.set_ylabel("Components with TGCAG-containing\nupstream-20-mer proxy (%)", fontsize=7.2)
    ax.set_ylim(0, 23); ax.set_xlim(0.4, 10.6)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(os.path.join(OUT, "panel_B_umotif_gradient.pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- panel C
def panel_C():
    """Genomic-position composition across all ten deciles.

    A composition (four categories summing to 100% within each decile) on a linear
    percentage axis. An earlier odds-ratio forest plot was abandoned after three
    failed attempts at a readable log axis: the chart type was the problem.
    """
    d = read_csv_skip_comments(os.path.join(FIG4, "tables/Table_S9_lost_site_location.csv"))
    rank = np.array([int(r["decile_rank"]) for r in d])
    series = {p: np.array([float(r["pct_" + p.lower()]) for r in d])
              for p in ["Promoter", "Intron", "Distal", "Other"]}

    fig, ax = plt.subplots(figsize=(4.7, 2.9))
    for pos in ["Promoter", "Distal", "Intron", "Other"]:
        ax.plot(rank, series[pos], "-o", color=POS_COL[pos], lw=1.4, ms=3.8,
                mec="white", mew=0.6, label=pos, zorder=3)
    for pos, yoff in [("Promoter", 2.6), ("Distal", -3.4), ("Intron", 2.4), ("Other", -3.0)]:
        v = series[pos]
        ax.text(10.25, v[-1] + yoff * 0.4, f"{v[-1]:.0f}%", fontsize=6.3,
                color=POS_COL[pos], va="center", ha="left", weight="bold")
        ax.text(0.72, v[0], f"{v[0]:.0f}%", fontsize=6.3,
                color=POS_COL[pos], va="center", ha="right", weight="bold")
    ax.set_xticks(rank); ax.set_xticklabels([str(i) for i in rank], fontsize=6.4)
    ax.set_xlabel("CTCF union components ranked by occupancy change, decile\n"
                  "(1 = most binding lost   \u2192   10 = binding gained)", fontsize=7.2)
    ax.set_ylabel("Share of components in each decile (%)", fontsize=7.2)
    ax.set_xlim(0.1, 11.0); ax.set_ylim(0, 66)
    ax.legend(frameon=False, fontsize=6.3, loc="upper center", bbox_to_anchor=(0.5, 1.16),
              ncol=4, handlelength=1.1, columnspacing=1.3)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(os.path.join(OUT, "panel_C_lost_site_location.pdf"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- panel D
def panel_D():
    """cPcdh locus: 4C contacts, CTCF occupancy and expression on one axis.

    Region assignment in the source CSV is by COORDINATE. The A4 table's
    nearest_pcdh_gene column reports "Pcdha4b" for 32 of its 41 peaks and is
    unusable. Subfamily intervals nest, so a disjoint assignment is used with
    Pcdhb and Pcdhg taking precedence.
    """
    d = read_csv_skip_comments(os.path.join(FIG4, "tables/Table_S10_cPcdh_locus.csv"))
    c4 = [r for r in d if r["track"] == "4C"]
    ctcf = [r for r in d if r["track"] == "CTCF"]
    expr = [r for r in d if r["track"] == "expression"]

    # Subfamily bands come from the ANNOTATION, parsed from the CSV header, not from
    # feature extents. Deriving them from the rows inflates Pcdha to 1037 kb (the first
    # and last feature do not sit at the interval edges) and overlaps the other bands.
    bounds = {}
    with open(os.path.join(FIG4, "tables/Table_S10_cPcdh_locus.csv")) as fh:
        for ln in fh:
            if not ln.lstrip('"').startswith("#"):
                break
            parts = ln.lstrip('"').lstrip("#").split()
            if len(parts) == 3 and parts[0] in ("Pcdha", "Pcdhb", "Pcdhg"):
                bounds[parts[0]] = [int(parts[1]), int(parts[2])]
    if set(bounds) != {"Pcdha", "Pcdhb", "Pcdhg"}:
        raise SystemExit("subfamily interval bounds not found in figure4_cpcdh_locus.csv header")
    b_s, b_e = bounds["Pcdhb"]
    x0, x1 = 36.90, 37.88

    fig = plt.figure(figsize=(6.0, 4.6))
    gs = fig.add_gridspec(4, 1, height_ratios=[0.26, 1.0, 1.0, 1.0], hspace=0.26)

    axr = fig.add_subplot(gs[0])
    for reg in ["Pcdha", "Pcdhb", "Pcdhg"]:
        s, e = bounds[reg]
        axr.axvspan(s / 1e6, e / 1e6, color=REG_COL[reg], alpha=0.55, lw=0)
        axr.text((s + e) / 2e6, 0.5, reg, ha="center", va="center", fontsize=7.0,
                 color="white", weight="bold", style="italic")
    axr.set_xlim(x0, x1); axr.set_ylim(0, 1); axr.axis("off")

    def shade(a):
        a.axvspan(b_s / 1e6, b_e / 1e6, color="#B03A2E", alpha=0.07, lw=0, zorder=0)

    # 4C contact change
    ax0 = fig.add_subplot(gs[1]); shade(ax0)
    xs = np.array([(int(r["chr18_start"]) + int(r["chr18_end"])) / 2e6 for r in c4])
    ds = np.array([float(r["value"]) for r in c4])
    o = np.argsort(xs); xs, ds = xs[o], ds[o]
    ax0.fill_between(xs, 0, ds, where=ds < 0, color="#C0453B", alpha=0.55, lw=0, interpolate=True)
    ax0.fill_between(xs, 0, ds, where=ds >= 0, color="#4E7A9B", alpha=0.55, lw=0, interpolate=True)
    ax0.axhline(0, color="0.3", lw=0.8)
    ax0.set_xlim(x0, x1)
    ax0.set_ylabel("4C contact change\nmut \u2212 WT", fontsize=6.6)
    ax0.set_xticklabels([]); ax0.spines[["top", "right"]].set_visible(False)
    ax0.text(x1 - 0.01, ax0.get_ylim()[1] * 0.72, "lost in mutant",
             fontsize=5.6, color="#C0453B", ha="right")

    # CTCF occupancy: lost peaks visually dominant
    ax1 = fig.add_subplot(gs[2]); shade(ax1)
    for r in ctcf:
        mid = (int(r["chr18_start"]) + int(r["chr18_end"])) / 2e6
        val = float(r["value"]); cls = r["peak_class"]
        col = CLS_COL["LOST"] if cls == "LOST" else ("#E1A73E" if cls == "GAINED" else "#7E8B99")
        is_lost = cls == "LOST"
        ax1.vlines(mid, 0, val, color=col, lw=1.4 if is_lost else 0.9,
                   alpha=1.0 if is_lost else 0.55, zorder=4 if is_lost else 2)
        ax1.plot(mid, val, "o", color=col, ms=4.0 if is_lost else 2.8, mec="white",
                 mew=0.5, alpha=1.0 if is_lost else 0.55, zorder=5 if is_lost else 3)
    ax1.axhline(0, color="0.3", lw=0.8)
    ax1.set_xlim(x0, x1); ax1.set_ylim(-4.6, 2.4)
    # pin ticks: matplotlib's automatic choice varies with version and would make
    # the regenerated panel differ cosmetically from the installed one
    ax1.set_yticks([-4, -2, 0, 2])
    ax1.set_ylabel("CTCF signal change\nlog$_2$(mut / WT)", fontsize=6.6)
    ax1.set_xticklabels([]); ax1.spines[["top", "right"]].set_visible(False)
    ax1.legend(handles=[Line2D([], [], color=CLS_COL["LOST"], marker="o", ls="none", ms=4.0, label="lost"),
                        Line2D([], [], color="#7E8B99", marker="o", ls="none", ms=2.8, label="maintained"),
                        Line2D([], [], color="#E1A73E", marker="o", ls="none", ms=2.8, label="gained")],
               frameon=False, fontsize=5.8, loc="lower left", ncol=3,
               handletextpad=0.3, columnspacing=0.9, bbox_to_anchor=(0.0, -0.04))

    # expression
    ax2 = fig.add_subplot(gs[3]); shade(ax2)
    for r in expr:
        mid = (int(r["chr18_start"]) + int(r["chr18_end"])) / 2e6
        val = float(r["value"]); sig = str(r["significant"]).strip().lower() == "true"
        col = REG_COL[r["region"]]
        ax2.vlines(mid, 0, val, color=col, lw=1.0, alpha=0.9 if sig else 0.28, zorder=2)
        ax2.plot(mid, val, "o", color=col, ms=3.2, mec="white", mew=0.5,
                 alpha=1.0 if sig else 0.32, zorder=3)
    ax2.axhline(0, color="0.3", lw=0.8)
    ax2.set_xlim(x0, x1)
    ax2.set_yticks([-6, -4, -2, 0, 2]); ax2.set_ylim(-7.4, 2.2)
    ax2.set_ylabel("Expression change\nlog$_2$(mut / WT)", fontsize=6.6)
    ax2.set_xlabel("chr18 position (Mb)", fontsize=7.0)
    ax2.spines[["top", "right"]].set_visible(False)

    fig.savefig(os.path.join(OUT, "panel_D_cpcdh_locus.pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    panel_A(); panel_B(); panel_C(); panel_D()
    print("Figure 4 panels A-D written to", os.path.abspath(OUT))
