#!/usr/bin/env python3
"""
Supplementary: aggregate CTCF profile at each peak class, and the width of the peaks.

WHY THIS IS IN THE SUPPLEMENT AND NOT THE MAIN FIGURE. The earlier build paired these
metaprofiles with a cumulative distribution of per-site change. That CDF said the same thing
main-figure panel B now says, that lost and gained are quantitative shifts with nothing at
baseline, so it was duplication and has been dropped. What the metaprofile adds, and panel B
cannot, is peak SHAPE.

THE QUESTION IT ANSWERS. CTCF peaks are expected to be sharp, and the first version of this
figure made them look broad. That was a display artefact, not the data: the old tracks were
log2(CTCF/input) at 10 bp bins with --pseudocount 1, where most bins hold 0 to 2 reads, so
log2((x+1)/(y+1)) can only take values near 0, +/-1, +/-2, and filling a signed variable
rendered that quantisation as solid wedges. The peaks were never broad. Measured directly,
wild-type peaks have a median width of 331 bp and mutant 406 bp, with summits at the midpoint;
0.4% and 2.6% respectively exceed 1 kb. This panel shows the same thing as a profile: a narrow
summit on a flat flank, in every class and both genotypes.

WHAT THE PROFILES SHOW
    lost      3.09 -> 2.01     signal falls by about a third, the summit stays a summit
    shared    4.06 -> 4.11     unchanged, the internal control
    gained    1.66 -> 2.60     rises from a pre-existing peak, not from flat background

The gained panel is the one worth looking at twice. Its wild-type curve is already a peak,
1.66 over a flank near 0.05, which is the aggregate version of the point that gained sites are
pre-existing weak sites rather than new ones.

COLOUR. The mutant fill was red in the previous build, the same red as the lost-peak class
colour, so one colour carried two unrelated meanings. Genotypes are blue and purple here, as
in the main tracks panel, and outside the class palette.

ONE LIBRARY PER GENOTYPE. No per-site significance test is possible, and a uniform offset
between genotypes would not be interpretable. Computed with pyBigWig one interval at a time:
the earlier rtracklayer version returned windows in contig order rather than query order and
is quarantined.

Date: 2026-08-16
"""
import os, csv, collections
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel_io import save_twin

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
OUT = os.path.join(Z, "E18p5_clean/results/chipseq_2026/figures")
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

CLS = ["lost", "shared", "gained"]
CLS_COL = {"lost": "#BF443A", "shared": "#6E7B8B", "gained": "#E0A53D"}
GT_COL = {"wt": "#3B6FD4", "mut": "#7B4EA8"}          # never the class reds and oranges

prof = collections.defaultdict(list)
nsites = {}
for r in csv.DictReader(l for l in open(os.path.join(T, "F4G_metaprofile_v2.csv"))
                        if l.strip() and not l.startswith("#")):
    g = "mut" if r["genotype"] not in ("wt", "wild type") else "wt"
    prof[(g, r["cls"])].append((int(r["pos"]), float(r["log2"])))
    nsites[r["cls"]] = int(r["n_sites"])
for k in prof:
    prof[k] = np.array(sorted(prof[k]))

fig, axes = plt.subplots(1, 3, figsize=(7.8, 3.0), sharey=True)
fig.subplots_adjust(left=0.085, right=0.985, top=0.795, bottom=0.275, wspace=0.13)
rows = []
for ax, c in zip(axes, CLS):
    for g, lab in (("wt", "wild type"), ("mut", "mutant")):
        a = prof[(g, c)]
        ax.plot(a[:, 0] / 1000, a[:, 1], color=GT_COL[g], lw=1.7, zorder=3, label=lab)
        peak = a[:, 1].max()
        flank = np.median(np.concatenate([a[:12, 1], a[-12:, 1]]))
        rows.append(dict(cls=c, genotype=g, n_sites=nsites[c], peak=round(float(peak), 3),
                         flank=round(float(flank), 4),
                         ratio_peak_flank=round(float(peak / max(flank, 1e-6)), 1)))
    w, m = prof[("wt", c)][:, 1].max(), prof[("mut", c)][:, 1].max()
    ax.set_title(f"{c}    n = {nsites[c]:,}", fontsize=8.8, color=CLS_COL[c],
                 fontweight="bold", pad=6)
    ax.text(0.03, 0.95, f"peak {w:.2f} → {m:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.0, color="0.3")
    ax.axvline(0, color="0.8", lw=0.6, ls=(0, (3, 3)), zorder=1)
    ax.set_xlim(-2, 2)
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.set_xticklabels(["−2", "−1", "0", "1", "2 kb"])
    ax.tick_params(length=2)
axes[0].set_ylabel("CTCF signal (log2 over input)", fontsize=8.4)
axes[0].legend(frameon=False, fontsize=7.4, loc="center left", bbox_to_anchor=(0.02, 0.66),
               handlelength=1.2)
axes[1].set_xlabel("distance from peak centre", fontsize=8.6)
sup = fig.suptitle("CTCF peaks are narrow in both genotypes and every class",
             fontsize=10.4, x=0.085, ha="left", y=0.965)
cap = fig.text(0.085, 0.020,
         "Measured directly, peaks have a median width of 331 bp in wild type and 406 bp in "
         "the mutant, with summits at the midpoint; 0.4% and 2.6% exceed 1 kb.\nAn earlier "
         "figure made them look broad because it plotted log2 ratios at 10 bp bins with a "
         "pseudocount of 1, where most bins hold 0 to 2 reads.\nNote the gained panel: its "
         "wild-type curve is already a peak at 1.66 over a flank near 0.05, so these are "
         "pre-existing weak sites that strengthen.",
         ha="left", va="bottom", fontsize=6.4, color="0.42")
save_twin(fig, OUT, "S4_metaprofile_by_class", captions=[cap, sup],
          bare_size=(7.8, 2.75), bare_adjust=dict(top=0.865, bottom=0.155))

with open(os.path.join(T, "S4_metaprofile_by_class.csv"), "w", newline="") as f:
    f.write("# Supplementary: aggregate CTCF profile per peak class, summit and flank.\n")
    f.write("# Source F4G_metaprofile_v2.csv (pyBigWig, order-safe; the rtracklayer version\n")
    f.write("#   returned windows in contig order and is quarantined).\n")
    f.write("# The cumulative-distribution half of the earlier panel was dropped: it\n")
    f.write("#   duplicated main-figure panel B. What survives is peak SHAPE.\n")
    f.write("# Peaks are narrow: median width 331 bp wild type, 406 bp mutant, summits at\n")
    f.write("#   the midpoint, 0.4% and 2.6% above 1 kb.\n")
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
for r in rows:
    print(f"[154] {r['cls']:<7}{r['genotype']:<5}n={r['n_sites']:>6,}  peak {r['peak']:.2f}  "
          f"flank {r['flank']:.3f}  peak/flank {r['ratio_peak_flank']:.0f}x")
print("[154] wrote S4_metaprofile_by_class.pdf and S4_metaprofile_by_class.csv")
