#!/usr/bin/env python3
"""
Figure 4 panels C and D, split into one message each, plus the supplementary position panel.

THE SPLIT. Script 139 drew position and anchor status as one two-part panel, and script 145
then added a separate panel showing the change in CTCF signal by position. That put the same
association in the figure twice: "promoters gain signal" and "gained calls are
promoter-enriched" are one fact viewed from two sides. The layout now separates the three
messages and demotes the duplicated one.

  panel C  genomic position by class            composition, main figure
  panel D  CTCF anchor status by class          a separate axis, main figure
  supp     change in CTCF signal by position    the identified measure, supplement

WHICH NUMBERS CAN CARRY A CLAIM. Peak classes are thresholds on a noisy difference, so a
class comparison conditions on a collider and is also confounded by site strength. Each
property was therefore tested twice, matched on wild-type signal and matched on mutant
signal; a real effect points the same way both times (script 140).

  promoter, lost vs shared      0.264 and 0.239    robust
  CTCF anchor, gained vs shared 0.525 and 0.132    robust
  promoter, gained vs shared    0.512 and 1.921    NOT IDENTIFIED
  CTCF anchor, lost vs shared   0.979 and 1.589    dissolves, no depletion at matched signal

So panel C's lost-side promoter depletion is solid and its gained-side promoter enrichment is
not; panel D's gained-site anchor deficit is solid and its lost-site deficit is not. Both
captions say so, and the supplementary panel carries the measure that does not depend on peak
calls at all.

WHAT PANEL D'S VARIABLE ACTUALLY IS, AND WHY IT WAS RENAMED. It was called "CTCF anchor",
which the provenance does not support.

  source      bonev_neural_CTCF_boundaries_merged.csv
  built from  Bonev et al. 2017, GSE96107, CTCF ChIP-seq narrowPeak from NPC (neural
              progenitor) AND CN (cortical neuron), merged ACROSS the two cell types
  contents    31,828 merged regions, median width 6 bp (59% are <=10 bp), 0.61% of the
              genome, median 2 contributing peaks per region, max 25
  used in     Chapter 1 as a PROXY for neural TAD boundaries, and its own integration plan
              describes it as "the CTCF-ChIP proxy boundaries" to be replaced by real Hi-C

So it is a two-cell-type merged CTCF peak set, deliberately assembled as a boundary proxy.
It is not a Hi-C boundary call, not an insulation measure, and not a loop-anchor annotation:
mouse has a few thousand TADs, not thirty thousand. For THIS figure the defensible reading is
the direct one, "is this also a CTCF site in independent neural ChIP data", and the panel is
labelled that way. Any architectural reading would have to come from the real 4DN Hi-C
boundary calls, which are on disk at public_3Dgenome_4DN/bonev_2017/4dn/ (NPC boundaries,
7,510 calls at 5 kb) and are NOT what this panel uses.

TWO LIMITS TO STATE IN THE TEXT. The Bonev cells are ES-derived neural progenitors and
cortical neurons grown in vitro, not E18.5 forebrain, so this is a system-matched-ish but not
tissue-matched comparison. And the exact merge operation is not documented in the repository:
only consumers of the file are tracked, not its build script, so the inputs and the file's
properties are known but the merge step is inferred from its contents.

It is also not a restatement of genomic position. Detection rates run 46.8% at promoters,
69.7% exon, 74.7% intron, 82.8% distal, and within distal intergenic sites alone the class
difference is undiminished: gained 37.5% against shared 92.1%.

CAVEAT FOR THE TEXT. Gained sites are weak in our own wild type (median 1.53 against 3.68
for shared), and weak sites are less likely to be called in any dataset, so this partly
overlaps with "gained sites are weak". It is not purely circular, since the deficit survives
matching on our own wild-type signal (MH OR 0.525, and 0.132 matched on mutant), but the two
facts are related and must not be presented as independent evidence.

Date: 2026-08-15
"""
import os, csv
import numpy as np
from scipy.stats import fisher_exact
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel_io import save_twin
from matplotlib.patches import Patch

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


def rd(name):
    return csv.DictReader(l for l in open(os.path.join(T, name))
                          if l.strip() and not l.startswith("#"))


def wilson(k, n, z=1.96):
    p, d = k / n, 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return 100 * (c - h), 100 * (c + h)


rows = list(rd("F4_peak_classes_with_position.csv"))
site = {c: [r for r in rows if r["cls"] == c] for c in CLS}
N = {c: len(v) for c, v in site.items()}
ypos = {c: i for i, c in enumerate(reversed(CLS))}       # lost at the top

# ================================================================ panel C: position
fig, ax = plt.subplots(figsize=(4.6, 2.5))
fig.subplots_adjust(left=0.215, right=0.975, top=0.79, bottom=0.31)
for c in CLS:
    left = 0.0
    for p in POS:
        v = 100 * sum(r["position"] == p for r in site[c]) / N[c]
        ax.barh(ypos[c], v, left=left, height=0.62, color=PCOL[p],
                edgecolor="white", linewidth=0.7, zorder=3)
        if v >= 6:
            ax.text(left + v / 2, ypos[c], f"{v:.0f}", ha="center", va="center",
                    fontsize=7.0, zorder=4,
                    color="white" if p in ("Promoter", "Exon") else "#2A2A2E")
        left += v
ax.set_xlim(0, 100)
ax.set_xlabel("% of sites", fontsize=8.6)
ax.set_yticks([ypos[c] for c in CLS])
ax.set_yticklabels([f"{c}\nn = {N[c]:,}" for c in CLS], fontsize=7.8)
for t, c in zip(ax.get_yticklabels(), CLS):
    t.set_color(COL[c])
ax.tick_params(length=2)
ax.spines["left"].set_visible(False)
ax.legend(handles=[Patch(facecolor=PCOL[p], edgecolor="white", label=p) for p in POS],
          frameon=False, fontsize=7.2, ncol=4, loc="upper center",
          bbox_to_anchor=(0.5, 1.26), handlelength=1.0, columnspacing=1.0,
          handletextpad=0.45)
capC = fig.text(0.5, 0.015, "Lost sites are distal and intronic; gained sites are promoter-proximal. "
         "Composition only:\nthe change in signal at each kind of location is in the "
         "supplement.", ha="center", va="bottom", fontsize=6.7, color="0.42")
save_twin(fig, OUT, "F4C_position", captions=[capC],
          bare_size=(4.6, 2.1), bare_adjust=dict(top=0.79, bottom=0.20))

# PANEL D IS NO LONGER BUILT HERE. It used the bonev_neural_CTCF_boundaries_merged
# proxy, which is a merged CTCF PEAK set rather than a boundary call. Script 147
# rebuilt it on the real 4DN Hi-C boundaries and it then moved to the supplement.
# The code is removed rather than left dormant so a re-run cannot resurrect it.

# ================================================================ supplement: signal change
S = {}
for r in rd("F4G_per_site_signal.csv"):
    try:
        w, m = float(r["wt_centre"]), float(r["mut_centre"])
    except (TypeError, ValueError):
        continue
    if not (np.isnan(w) or np.isnan(m)):
        S[(r["chr"], r["start"], r["end"])] = (w, m)
posn = {(r["chr"], r["start"], r["end"]): r["position"] for r in rows}
K = [k for k in S if k in posn]
fig, ax = plt.subplots(figsize=(4.2, 2.6))
fig.subplots_adjust(left=0.315, right=0.965, top=0.93, bottom=0.30)
rowsS = []
for i, p in enumerate(POS):
    v = np.array([S[k][1] - S[k][0] for k in K if posn[k] == p])
    se = 1.96 * v.std(ddof=1) / np.sqrt(len(v))
    y = len(POS) - 1 - i
    ax.barh(y, v.mean(), height=0.6, color=PCOL[p], zorder=3,
            edgecolor="#9AA3AC" if p == "Distal intergenic" else "none", linewidth=0.5)
    ax.plot([v.mean() - se, v.mean() + se], [y, y], color="#2A2A2E", lw=1.0, zorder=4,
            solid_capstyle="butt")
    ax.text(0.755, y, f"{v.mean():+.2f}", va="center", ha="right", fontsize=7.4, zorder=4)
    rowsS.append(dict(position=p, n=len(v), mean_change=round(float(v.mean()), 4),
                      ci_low=round(float(v.mean() - se), 4),
                      ci_high=round(float(v.mean() + se), 4)))
    print(f"[146] supp  {p:<18} n={len(v):>6,}  change {v.mean():+.3f}")
ax.axvline(0, color="#2A2A2E", lw=0.7, zorder=2)
ax.set_yticks(range(len(POS)))
ax.set_yticklabels(list(reversed(POS)), fontsize=8)
ax.set_xlabel("change in CTCF signal (log2, mutant minus wild type)", fontsize=8.2)
ax.set_xlim(-0.20, 0.76)
ax.set_xticks([-0.2, 0, 0.2, 0.4, 0.6])
ax.tick_params(length=2)
ax.spines["left"].set_visible(False)
capS = fig.text(0.5, 0.015, "The gain tracks proximity to the transcription start site. No peak "
         "classes are used,\nso this does not depend on a calling threshold. Bars are 95% "
         "confidence intervals.",
         ha="center", va="bottom", fontsize=6.7, color="0.42")
save_twin(fig, OUT, "S4_signal_change_by_position", captions=[capS],
          bare_size=(4.2, 2.25), bare_adjust=dict(top=0.95, bottom=0.235))
with open(os.path.join(T, "S4_signal_change_by_position.csv"), "w", newline="") as f:
    f.write("# Supplementary: change in CTCF signal by genomic position category.\n")
    f.write("# Continuous outcome at every site, no peak classes, so nothing is conditioned\n")
    f.write("#   on the outcome. This is the identified version of the promoter result that\n")
    f.write("#   main-figure panel C shows as composition.\n")
    f.write("# Further validation: matched on core motif score the promoter difference is\n")
    f.write("#   +0.546 in 5/5 strata, and promoters with no CTCF peak sit at 0.08 CPM and\n")
    f.write("#   do not rise, so this is not ChIP background at open chromatin (script 141).\n")
    w_ = csv.DictWriter(f, fieldnames=list(rowsS[0].keys()))
    w_.writeheader(); w_.writerows(rowsS)

for r in rowsD:
    print(f"[146] D  {r['cls']:<7} n={r['n']:>6,}  anchor {r['pct']:>5.1f}%  "
          f"OR vs shared {r['odds_ratio_vs_shared']:.3f}")
print("[146] wrote F4C_position.pdf, F4D_independent_ctcf.pdf, S4_signal_change_by_position.pdf")
print("[146] SUPERSEDED: F4C_annotation.pdf (two-part), F4D_promoter_gain.pdf (moved to supp)")
