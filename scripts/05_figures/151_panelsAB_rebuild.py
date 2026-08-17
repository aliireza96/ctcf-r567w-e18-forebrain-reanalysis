#!/usr/bin/env python3
"""
Figure 4 panels A and B, rebuilt. These were the only main panels untouched in the rebuild.

WHAT WAS WRONG WITH EACH.

PANEL A drew shared in BLUE while panels C, D and E all use grey for shared, so the figure's
palette broke at its first panel. It also stated only the three class counts, which misses
the fact that makes "redistribution" the right word rather than "loss": the total number of
CTCF peaks barely moves, 33,192 in wild type against 33,614 in the mutant, +1.3%. The mutant
does not have less CTCF binding; it has CTCF binding in different places.

  THE BOOKKEEPING THAT HAS TO BE STATED. Shared sites are counted in WILD-TYPE intervals.
  Adding shared + gained gives 33,950 against 33,614 mutant peaks called, a gap of 336,
  because one mutant peak can overlap several wild-type intervals. Counted the other way
  round shared is ~26,405 rather than 26,741. Both numbers are correct; they answer slightly
  different questions. A reader who adds the panel's numbers and compares them to the peak
  files will hit this, so the panel says which convention it uses.

PANEL B still drew a dotted line at 0.5 log2, labelled as a threshold. That is the cut from
the four-way subclass split ("de novo" / "strengthened", "fully lost" / "weakened") which was
tested and REJECTED on three independent grounds: Hartigan dip p = 1.0000 and Silverman
p = 0.635 on wild-type signal among gained sites, so the distribution is unimodal; continuous
signal explains 1.8 to 2.9 times more deviance than the step; and roughly 403 of the 732
"de novo" sites are attributable to regression to the mean. Leaving the line in implies a
threshold the figure no longer uses, and invites exactly the reading that was retired.

WHAT PANEL B IS FOR, and it is the load-bearing panel of the figure. Every downstream claim
depends on lost and gained being QUANTITATIVE shifts rather than presence and absence: lost
sites keep a median 1.76 in the mutant against a flank baseline near 0.05, and gained sites
already carry 1.53 in wild type. Shared sites do not move, 3.68 to 3.75, which is the internal
control showing the two libraries are comparable.

GENOTYPE ENCODING DIFFERS FROM PANEL E, DELIBERATELY. Here class is the primary grouping and
genotype is secondary, so genotype is light against dark within each class colour. In panel E
genotype is primary and class is an annotation, so genotype takes the blue/purple fills. Both
legends state their own convention.

Date: 2026-08-16
"""
import os, csv
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel_io import save_twin
from matplotlib.patches import Patch

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
D = os.path.join(Z, "E18p5_clean/data/ctcf_reanalysis_2026")
OUT = os.path.join(Z, "E18p5_clean/results/chipseq_2026/figures")
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

CLS = ["lost", "shared", "gained"]
COL = {"lost": "#BF443A", "shared": "#6E7B8B", "gained": "#E0A53D"}   # matches C, D, E

rows = [r for r in csv.DictReader(
    l for l in open(os.path.join(T, "F4_peak_classes_with_position.csv"))
    if l.strip() and not l.startswith("#"))]
N = {c: sum(1 for r in rows if r["cls"] == c) for c in CLS}
n_wt = sum(1 for _ in open(os.path.join(D, "wt_PRIMARY_peaks.narrowPeak")))
n_mut = sum(1 for _ in open(os.path.join(D, "homo_PRIMARY_peaks.narrowPeak")))
print(f"[151] classes: " + "  ".join(f"{c} {N[c]:,}" for c in CLS) + f"   union {len(rows):,}")
print(f"[151] peaks called: wild type {n_wt:,}, mutant {n_mut:,}  "
      f"({100*(n_mut-n_wt)/n_wt:+.1f}%)")

# ================================================================ panel A
fig, (axt, axb) = plt.subplots(
    2, 1, figsize=(5.6, 3.3), gridspec_kw=dict(height_ratios=[1, 1.5], hspace=0.60))
fig.subplots_adjust(left=0.150, right=0.975, top=0.92, bottom=0.40)

# top: total peaks called per genotype. The point is that it barely moves.
for i, (lab, n, col) in enumerate((("wild type", n_wt, "#3B6FD4"),
                                   ("mutant", n_mut, "#7B4EA8"))):
    axt.barh(1 - i, n, height=0.55, color=col, zorder=3)
    axt.text(n + 400, 1 - i, f"{n:,}", va="center", fontsize=7.4, color=col)
axt.set_yticks([1, 0]); axt.set_yticklabels(["wild type", "mutant"], fontsize=7.6)
axt.set_xlim(0, 44000); axt.set_xticks([])
axt.set_ylim(-0.55, 1.55)
axt.tick_params(length=0)
axt.spines["bottom"].set_visible(False); axt.spines["left"].set_visible(False)
ttl_a = axt.set_title("peaks called", fontsize=8.2, pad=3, loc="left", color="0.35")

# bottom: the union, partitioned
left = 0
for c in CLS:
    axb.barh(0, N[c], left=left, height=0.62, color=COL[c], zorder=3,
             edgecolor="white", linewidth=0.8)
    axb.text(left + N[c] / 2, 0, f"{N[c]:,}", ha="center", va="center", fontsize=7.6,
             color="white", zorder=4)
    left += N[c]
axb.set_ylim(-0.6, 0.6); axb.set_yticks([])
axb.set_xlim(0, 44000)
axb.set_xticks([0, 10000, 20000, 30000, 40000])
axb.set_xticklabels(["0", "10,000", "20,000", "30,000", "40,000"])
axb.set_xlabel("CTCF sites", fontsize=8.6)
axb.tick_params(length=2)
axb.spines["left"].set_visible(False)
ttl_b = axb.set_title("union of both peak sets", fontsize=8.2, pad=3, loc="left", color="0.35")
axb.legend(handles=[Patch(facecolor=COL[c], label=c) for c in CLS],
           frameon=False, fontsize=7.4, ncol=3, loc="upper center",
           bbox_to_anchor=(0.5, -0.62), handlelength=1.1, columnspacing=1.6)
capA = fig.text(0.5, 0.012,
         f"Total CTCF binding is not reduced: {n_wt:,} peaks in wild type against {n_mut:,} in "
         f"the mutant, {100*(n_mut-n_wt)/n_wt:+.1f}%.\nWhat changes is where. Shared sites are "
         "counted in wild-type intervals; counted in mutant\nintervals they number ~26,405, "
         "since one mutant peak can overlap several wild-type ones.",
         ha="center", va="bottom", fontsize=6.5, color="0.42")
save_twin(fig, OUT, "F4A_site_classes", captions=[capA, ttl_a, ttl_b],
          bare_size=(5.6, 2.5), bare_adjust=dict(top=0.90, bottom=0.28))

# ================================================================ panel B
S = {}
for r in csv.DictReader(l for l in open(os.path.join(T, "F4G_per_site_signal.csv"))
                        if l.strip() and not l.startswith("#")):
    try:
        w, m = float(r["wt_centre"]), float(r["mut_centre"])
    except (TypeError, ValueError):
        continue
    if not (np.isnan(w) or np.isnan(m)):
        S.setdefault(r["cls"], []).append((w, m))

fig, ax = plt.subplots(figsize=(4.6, 4.2))
fig.subplots_adjust(left=0.150, right=0.965, top=0.92, bottom=0.245)
pos, ticks, labs, out = 0, [], [], []
for c in CLS:
    arr = np.array(S[c])
    for gi, g in enumerate(("wild type", "mutant")):
        v = arr[:, gi]
        p = ax.violinplot([v], positions=[pos], widths=0.72, showextrema=False,
                          showmedians=True)
        for b in p["bodies"]:
            b.set_facecolor(COL[c]); b.set_alpha(0.42 if gi == 0 else 0.92)
            b.set_edgecolor("none")
        p["cmedians"].set_color("white"); p["cmedians"].set_linewidth(1.3)
        ax.text(pos, np.median(v) + 0.16, f"{np.median(v):.2f}", ha="center", va="bottom",
                fontsize=6.6, color="0.2")
        out.append(dict(cls=c, genotype=g, n=len(v), median=round(float(np.median(v)), 3),
                        q25=round(float(np.percentile(v, 25)), 3),
                        q75=round(float(np.percentile(v, 75)), 3)))
        ticks.append(pos); labs.append(g.replace(" ", "\n"))
        pos += 1
    pos += 0.6
# the flank baseline is the only horizontal reference this panel needs. The 0.5 line that
# used to sit here was the four-way subclass threshold, and that split was rejected.
ax.axhline(0.05, color="0.55", lw=0.8, ls=(0, (4, 3)), zorder=1)
ax.text(0.985, 0.10, "flank baseline", transform=ax.get_yaxis_transform(), va="bottom",
        ha="right", fontsize=6.2, color="0.5")
ax.set_xticks(ticks); ax.set_xticklabels(labs, fontsize=6.8)
ax.set_ylabel("CTCF signal (log2 over input)", fontsize=8.4)
ax.set_ylim(-0.5, 7.4)
ax.tick_params(length=2)
for c, x in (("lost", 0.5), ("shared", 3.1), ("gained", 5.7)):
    ax.text(x, 7.25, c, ha="center", va="top", fontsize=8.4, color=COL[c], fontweight="bold")
capB = fig.text(0.5, 0.010,
         "Nothing starts or ends at baseline: lost sites keep a median 1.76 in\nthe mutant, "
         "and gained sites already carry 1.53 in wild type. Shared\nsites do not move, 3.68 "
         "to 3.75, so the two libraries are comparable.",
         ha="center", va="bottom", fontsize=6.5, color="0.42")
save_twin(fig, OUT, "F4B_signal_distributions", captions=[capB],
          bare_size=(4.6, 3.7), bare_adjust=dict(top=0.92, bottom=0.135))

with open(os.path.join(T, "F4AB_source.csv"), "w", newline="") as f:
    f.write("# Panels A and B source.\n")
    f.write(f"# Peaks called: wild type {n_wt:,}, mutant {n_mut:,} ({100*(n_mut-n_wt)/n_wt:+.1f}%).\n")
    f.write(f"# Union partition: lost {N['lost']:,}, shared {N['shared']:,}, "
            f"gained {N['gained']:,}, total {len(rows):,}.\n")
    f.write("# Shared is counted in WILD-TYPE intervals. shared + gained = 33,950 against\n")
    f.write("#   33,614 mutant peaks called: 336 mutant peaks each overlap several wild-type\n")
    f.write("#   intervals. Counted in mutant intervals shared is ~26,405. Both are correct.\n")
    f.write("# Panel B carries NO 0.5 threshold line: the four-way subclass split it belonged\n")
    f.write("#   to was rejected (dip p = 1.0000, Silverman p = 0.635, and ~403 of 732 'de\n")
    f.write("#   novo' sites attributable to regression to the mean).\n")
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)
for r in out:
    print(f"[151] B  {r['cls']:<7}{r['genotype']:<11}n={r['n']:>6,}  median {r['median']:.2f}")
print("[151] wrote F4A_site_classes.pdf, F4B_signal_distributions.pdf, F4AB_source.csv")
