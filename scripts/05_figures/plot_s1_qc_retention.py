#!/usr/bin/env python3
"""Plot the retained fraction of input barcodes for Supplementary Figure S1F."""
import csv
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt


fig_dir = sys.argv[1] if len(sys.argv) > 1 else ".."
out_dir = sys.argv[2] if len(sys.argv) > 2 else fig_dir
table_path = os.path.join(fig_dir, "tables", "S1_qc_retention.csv")

with open(table_path, newline="") as handle:
    rows = list(csv.DictReader(handle))

labels = [row["genotype"] for row in rows]
percent = [float(row["retained_percent"]) for row in rows]
counts = [(int(row["retained_nuclei"]), int(row["input_barcodes"])) for row in rows]

mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
})

fig, ax = plt.subplots(figsize=(3.35, 2.25))
colours = ["#3A6EA5", "#C0453B"]
bars = ax.bar(labels, percent, color=colours, width=0.62, edgecolor="none")
ax.set_ylim(0, 100)
ax.set_ylabel("Input barcodes retained (%)", fontsize=7.2)
ax.tick_params(axis="x", labelsize=7)
ax.tick_params(axis="y", labelsize=6.5)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="y", color="0.9", linewidth=0.5, zorder=0)
ax.set_axisbelow(True)

for bar, value, (retained, total) in zip(bars, percent, counts):
    ax.text(bar.get_x() + bar.get_width() / 2, value + 2.4,
            f"{value:.1f}%", ha="center", va="bottom", fontsize=7.2)
    ax.text(bar.get_x() + bar.get_width() / 2, value / 2,
            f"{retained:,} / {total:,}", ha="center", va="center",
            fontsize=6.4, color="white", fontweight="bold")

ax.text(0.5, -0.24, "Identical fixed QC thresholds", transform=ax.transAxes,
        ha="center", va="top", fontsize=6.4, color="0.35")
fig.savefig(os.path.join(out_dir, "panel_F_qc_retention.pdf"), bbox_inches="tight")
plt.close(fig)
