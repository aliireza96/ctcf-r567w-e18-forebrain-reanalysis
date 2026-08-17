#!/usr/bin/env python3
"""
Evf2 (Dlx6os1) transcript by cell type: the expression layer of the Dlx locus result.

WHY IT IS ITS OWN PANEL. It shares no axis with the locus tracks, which run in genomic
coordinates, and it is a different kind of evidence: single-nucleus differential expression
rather than ChIP coverage. Stacking it under the tracks made one panel carry two axes and two
data types.

WHICH CLASS SCHEME, AND WHY IT MATTERS. An earlier build used the older 13-type fine
annotation, whose only progenitor class is "Cycling Radial Glia (Dorsal/Gliogenic)". That
label is wrong for a ventral gene: subclustering the Cycling RG class shows 444 of its 1,147
cells are a Ventral IPC (LGE-leaning) state at 85% Dlx1-positive, so the +1.07 attributed to
a dorsal/gliogenic class was measured across a population that is 39% ventral.

This panel uses the UPDATED class scheme that Figures 2 and 3 use, in which Ventral IPC is its
own class. That is both the correct population and consistent with the rest of the chapter.

THE NUMBERS MOVE WHEN THE SCHEME CHANGES, AND NOT IN OUR FAVOUR. On the old scheme, striatal
D1 SPNs reached padj 7.9e-3. The updated scheme splits D1 into precursor-like and
differentiating, which dilutes it to padj 0.161. Only MGE-derived interneurons survives BH
correction here. Ventral IPC carries the largest effect (+0.838) at raw p = 4.0e-5 but does
not survive correction either.

WHAT THE PANEL CAN AND CANNOT SAY. Evf2 is tested in exactly eight classes, all ventral or
GABAergic, and in none of the dorsal or glial ones, which is the sanity check that it is being
measured where it is expressed. Every one of the eight moves in the same direction, up in the
mutant. One survives multiple-testing correction. That is a consistent direction with one
significant class, not eight significant classes, and the panel marks which is which.

AGGREGATION TRADE-OFF, RECORDED BECAUSE IT WILL COME UP. A broader scheme also exists in
which the two D1 classes are one "SPNs" class. There, SPNs reach +0.404 at padj 0.042, so
merging recovers a second significant class. It is not used here because that scheme has no
Ventral IPC: its only progenitor class is the mixed Cycling RG, 39% of which is Ventral IPC
by subclustering, which is the labelling problem this panel exists to avoid. Choosing the
correct population costs one significant class, and that is the right trade.

NOT ON THE PANEL, BUT FOR THE RESULTS TEXT: Dlx5 and Dlx6 themselves are unchanged in every
class tested. So this is not "Evf2 rises and drives Dlx5/Dlx6". And CTCF at a promoter can
activate or repress, so gained binding alongside a raised transcript is co-occurrence.

Date: 2026-08-16
"""
import os, csv, glob, re
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel_io import save_twin

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
DE = os.path.join(Z, "E18p5_clean/results/figure3_rebuild_122/de")
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
OUT = os.path.join(Z, "E18p5_clean/results/chipseq_2026/figures")
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

GENE = "Dlx6os1"
NICE = {"Ventral_IPC": "Ventral IPC",
        "MGE_IN": "MGE interneurons",
        "SPN_D1_differentiating_": "SPN D1, differentiating",
        "GABAergic_IN_origin_unassigned": "GABAergic IN, unassigned",
        "LGE_derived_GABAergic": "LGE-derived GABAergic",
        "CGE_IN": "CGE interneurons",
        "SPN_D1_precursor_like_": "SPN D1, precursor-like",
        "SPN_D2_differentiating_": "SPN D2, differentiating"}

rows = []
for f in sorted(glob.glob(os.path.join(DE, "122_DE_*_mut_vs_wt.csv"))):
    cls = re.sub(r"^122_DE_|_mut_vs_wt\.csv$", "", os.path.basename(f))
    for r in csv.DictReader(open(f)):
        if r["gene"] != GENE:
            continue
        rows.append(dict(cls=cls, label=NICE.get(cls, cls.replace("_", " ")),
                         log2FC=float(r["avg_log2FC"]), p=float(r["p_val"]),
                         padj=float(r["p_val_adj"]),
                         pct_mut=float(r["pct.1"]), pct_wt=float(r["pct.2"])))
        break
rows.sort(key=lambda r: r["log2FC"])
print(f"[153] {GENE} tested in {len(rows)} classes, "
      f"{sum(r['padj'] < 0.05 for r in rows)} significant after BH")

fig, ax = plt.subplots(figsize=(5.2, 3.5))
fig.subplots_adjust(left=0.315, right=0.965, top=0.885, bottom=0.290)
for i, r in enumerate(rows):
    sig = r["padj"] < 0.05
    ax.barh(i, r["log2FC"], height=0.62, zorder=3,
            color="#E0A53D" if sig else "#DEDEE3",
            edgecolor="#B8912F" if sig else "#C4C4CB", linewidth=0.6)
    txt = f"{r['log2FC']:+.2f}" + (f"   padj {r['padj']:.0e}" if sig else "")
    ax.text(r["log2FC"] + 0.02, i, txt, va="center", fontsize=6.8,
            color="#2A2A2E" if sig else "0.5")
ax.axvline(0, color="#2A2A2E", lw=0.7, zorder=2)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r["label"] for r in rows], fontsize=7.2)
for t, r in zip(ax.get_yticklabels(), rows):
    if r["padj"] < 0.05:
        t.set_color("#8A6B22")
ax.set_xlim(-0.06, 1.32)
ax.set_xlabel("$\\it{Evf2}$ expression change, log2 (mutant vs wild type)", fontsize=8.4)
ax.tick_params(length=2)
ax.spines["left"].set_visible(False)
ttl = ax.set_title("$\\it{Evf2}$ is raised in ventral populations", fontsize=9.6, pad=7,
                   loc="left")
cap = fig.text(0.030, 0.018,
         "Tested in the eight ventral and GABAergic classes where it is expressed, and in "
         "none of the dorsal\nor glial ones. All eight move the same way; one survives "
         "Benjamini-Hochberg correction (filled).\nUpdated class scheme, as in Figures 2 "
         "and 3.",
         ha="left", va="bottom", fontsize=6.5, color="0.42")
save_twin(fig, OUT, "F4H_evf2_expression", captions=[cap, ttl],
          bare_size=(5.2, 2.9), bare_adjust=dict(top=0.965, bottom=0.175))

with open(os.path.join(T, "F4H_evf2_expression_source.csv"), "w", newline="") as f:
    f.write("# Evf2 (Dlx6os1) mutant vs wild type, UPDATED class scheme (Figures 2 and 3),\n")
    f.write("#   from results/figure3_rebuild_122/de/. Tested in 8 ventral/GABAergic classes.\n")
    f.write("# The older 13-type scheme had no Ventral IPC class: its only progenitor class,\n")
    f.write("#   'Cycling Radial Glia (Dorsal/Gliogenic)', is 39% Ventral IPC by subclustering\n")
    f.write("#   (444 of 1,147 cells, 85% Dlx1+), so that label was wrong for a ventral gene.\n")
    f.write("# Scheme change costs significance: D1 SPNs were padj 7.9e-3 when D1 was one\n")
    f.write("#   class; split into precursor-like and differentiating it is padj 0.161.\n")
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
for r in reversed(rows):
    print(f"[153] {r['label']:<28}{r['log2FC']:+7.3f}  p {r['p']:.1e}  padj {r['padj']:.2e}"
          + ("  *" if r["padj"] < 0.05 else ""))
print("[153] wrote F4H_evf2_expression.pdf and F4H_evf2_expression_source.csv")
