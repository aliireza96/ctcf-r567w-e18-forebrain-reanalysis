#!/usr/bin/env python3
"""
Figure 4 panel E, final. (Panel D was built here too; see the note below.)

PANEL D. Script 144 drew this as a two-point slope, which was ugly and threw away the
genomic-position breakdown that makes panel C readable. It is redrawn here in panel C's
visual grammar: the same four position categories, the same colour ramp, the same horizontal
layout, so the two panels read as a pair. Panel C says where each class of site sits; panel D
says what happened to the CTCF signal at each kind of location. Promoter-proximal sites gain,
everything else is flat.

This is the identified measurement, not the peak-call one. The outcome is the continuous
change in signal at every site, with no peak classes involved, so nothing is conditioned on
the outcome. Validation is in the supplement, not on the panel: matched on core motif score
the promoter difference is +0.546 in 5 of 5 strata (script 141), and promoters carrying no
CTCF peak sit at 0.08 CPM and do not rise, so this is not ChIP background at open chromatin.

PANEL E. Restored to the script 134 design, which paired a schematic of the motif with the
bars. Two boxes on a 5' to 3' line: the upstream motif and the CTCF core, drawn in neutral
greys deliberately outside the class palette so the schematic says what the motif IS while
the bars say who carries it.

NAMING. "Element" is dropped throughout. It said nothing. The sequence is the UPSTREAM MOTIF,
TGCAG, sitting 7 bp upstream of the CTCF core; Do et al. call the same sequence the U motif
and tie its dependence to ZF9-11, where R567W sits.

The core box carries the consensus of the matrix the peaks were called with,
CTCF_MOUSE.H11MO.0.A, restricted to its informative span. Positions 1-2 and 18-20 sit below
0.3 bits; positions 3-17 are the 15 bp core and are what is printed. Taking all 20 columns
would put five near-random letters in the box and imply the core is longer than it is.

Date: 2026-08-15
"""
import os, csv
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel_io import save_twin
from matplotlib.patches import FancyBboxPatch

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
OUT = os.path.join(Z, "E18p5_clean/results/chipseq_2026/figures")
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

CLS = ["lost", "shared", "gained"]
COL = {"lost": "#BF443A", "shared": "#6E7B8B", "gained": "#E0A53D"}
POS = ["Promoter", "Exon", "Intron", "Distal intergenic"]
PCOL = {"Promoter": "#2C4A63", "Exon": "#5F87AB",
        "Intron": "#A3BDD1", "Distal intergenic": "#DFE7ED"}
REG, MOTIF = 6, "TGCAG"
CORE = "GCCACCAGGGGGCGC"
MONO = "DejaVu Sans Mono"


def rd(name):
    return csv.DictReader(l for l in open(os.path.join(T, name))
                          if l.strip() and not l.startswith("#"))


def wilson(k, n, z=1.96):
    p, d = k / n, 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


# PANEL D IS NO LONGER BUILT HERE. Its position-category bars became the supplementary
# panel (script 146) once panel D was rebuilt on real Hi-C boundaries (script 147).
# The code was removed rather than left dormant so that re-running this script cannot
# resurrect a retired panel into the live figures directory.

# ================================================================ panel E
rows = list(rd("F4C_umotif_panel_assignments.csv"))
by = {c: [r for r in rows if r["cls"] == c] for c in CLS}

fig = plt.figure(figsize=(3.9, 4.0))
gs = fig.add_gridspec(2, 1, height_ratios=[1, 3.8], hspace=0.30,
                      left=0.195, right=0.975, top=0.965, bottom=0.135)

axs = fig.add_subplot(gs[0])
axs.set_xlim(0, 10); axs.set_ylim(0, 1); axs.axis("off")
axs.plot([0.2, 9.8], [0.40, 0.40], color="0.55", lw=1.0, zorder=1)
axs.add_patch(FancyBboxPatch((1.05, 0.18), 2.15, 0.44,
                             boxstyle="round,pad=0.02,rounding_size=0.08",
                             facecolor="#3F3F46", edgecolor="none", zorder=2))
axs.text(2.125, 0.405, MOTIF, ha="center", va="center", fontsize=7.6, color="white",
         fontweight="bold", zorder=3, family=MONO)
axs.text(2.125, 0.755, "upstream motif", ha="center", va="center", fontsize=7.6, color="0.30")
axs.add_patch(FancyBboxPatch((3.95, 0.18), 5.35, 0.44,
                             boxstyle="round,pad=0.02,rounding_size=0.08",
                             facecolor="#D9D9DE", edgecolor="none", zorder=2))
axs.text(6.625, 0.405, CORE, ha="center", va="center", fontsize=7.6, color="#3F3F46",
         fontweight="bold", zorder=3, family=MONO)
axs.text(6.625, 0.755, "CTCF core motif", ha="center", va="center", fontsize=7.6, color="0.30")
axs.annotate("", xy=(3.93, 0.40), xytext=(3.22, 0.40),
             arrowprops=dict(arrowstyle="-", color="0.55", lw=1.0))
axs.text(3.575, 0.135, "7 bp", ha="center", va="center", fontsize=6.4, color="0.55")
axs.text(0.2, 0.40, "5′ ", ha="right", va="center", fontsize=7.4, color="0.55")
axs.text(9.8, 0.40, " 3′", ha="left", va="center", fontsize=7.4, color="0.55")

ax = fig.add_subplot(gs[1])
rowsE = []
for i, c in enumerate(CLS):
    n = len(by[c])
    k = sum(r["upstream"].upper()[REG:REG + 5] == MOTIF for r in by[c])
    pct = 100 * k / n
    lo, hi = wilson(k, n)
    ax.bar(i, pct, width=0.62, color=COL[c], zorder=3, linewidth=0)
    ax.errorbar(i, pct, yerr=[[pct - lo], [hi - pct]], color="0.25", lw=0.9,
                capsize=3, zorder=4)
    ax.text(i, hi + 0.55, f"{pct:.1f}", ha="center", va="bottom", fontsize=8.4,
            color=COL[c], fontweight="bold")
    rowsE.append(dict(cls=c, n=n, n_carrying=k, pct=round(pct, 2),
                      wilson_low=round(lo, 2), wilson_high=round(hi, 2)))
    print(f"[145] E  {c:<7} n={n:>6,}  {pct:.2f}%")
ax.set_xticks(range(3))
ax.set_xticklabels(CLS)
for t, c in zip(ax.get_xticklabels(), CLS):
    t.set_color(COL[c])
ax.set_ylabel("sites with upstream TGCAG (%)", fontsize=9)
ax.set_ylim(0, 17)
ax.set_xlim(-0.62, 2.62)
ax.tick_params(length=2.2)
# this panel carries no explanatory footnote, so the twin is identical. It is still
# emitted so assembly can pick up "<stem>_bare" for every panel without special cases.
save_twin(fig, OUT, "F4D_upstream_motif", captions=[])
with open(os.path.join(T, "F4D_upstream_motif_source.csv"), "w", newline="") as f:
    f.write("# Panel E: carriage of the UPSTREAM MOTIF (TGCAG, 7 bp upstream of the CTCF\n")
    f.write("#   core) by peak class. Do et al. call the same sequence the U motif.\n")
    f.write("# Validation in the supplement: the enrichment holds in every wild-type signal\n")
    f.write("#   bin, MH OR 2.16 matched on wild-type and 1.70 matched on mutant, so it is\n")
    f.write("#   not a restatement of lost sites being weaker than shared sites.\n")
    w_ = csv.DictWriter(f, fieldnames=list(rowsE[0].keys()))
    w_.writeheader(); w_.writerows(rowsE)
print("[145] wrote F4D_upstream_motif.pdf")
