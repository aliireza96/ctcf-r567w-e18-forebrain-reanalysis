#!/usr/bin/env python3
"""
The Dlx5/Dlx6 locus: a boundary that loses CTCF beside a regulatory lncRNA that gains it.

WHY THIS PANEL WAS NEARLY DROPPED, AND WHY IT SURVIVED. Its insulation change is 0.577 units
against an average of 0.043 for boundaries that lose a site, roughly the 99th percentile, and
the locus was chosen partly for that. An extreme case chosen for being extreme cannot support
a claim about the typical effect, and this panel does not make one.

What rescued it is that every element of it is independently confirmed, so it is a case study
rather than an anecdote:

  our peak      class    wt -> mut    Yu E16.5 fold   called in E12.5/E14.5/E16.5/P0/BonevNPC
  6.870 Mb      gained   0.70 -> 1.16      14.3       Y Y Y Y .
  6.885 Mb      shared   7.91 -> 7.21      27.5       Y Y Y Y Y
  6.915 Mb      shared   5.93 -> 6.11      17.3       Y Y Y Y Y
  6.943 Mb      shared   1.05 -> 0.99       5.6       Y Y Y . Y
  6.956 Mb      lost     2.44 -> 1.28      13.1       Y Y Y Y Y
  7.003 Mb      lost     3.14 -> 0.81      29.9       Y Y Y Y Y
  7.004 Mb      lost     2.33 -> 0.87      15.1       Y Y Y Y Y

  Bonev NPC Hi-C calls a STRONG boundary at 6.955 Mb, independently of our insulation track
  and of our peak calls.

So an outside Hi-C dataset says there is a strong boundary here, five outside ChIP datasets
say there is a CTCF site on it, and in the mutant that site drops by half while insulation
weakens.

THE GENE ASSIGNMENTS, WHICH THE EARLIER PANEL GOT HALF RIGHT. The gained promoter peak at
6.870 Mb is 0.5 kb from Dlx6os1, the Evf2 ultraconserved-enhancer lncRNA that regulates
Dlx5/Dlx6 and GABAergic interneuron development. The LOST promoter peak at 6.956 Mb is on
Sdhaf3, not on Dlx5: Dlx5 is 73.9 kb away, and the two intronic losses at 7.003 and 7.004 Mb
are 121 kb out. So the correct statement is NOT "CTCF is lost at Dlx5/Dlx6". It is that the
domain containing the Dlx cluster loses CTCF at its boundary and at sites beyond it, while
the cluster's own regulatory lncRNA gains CTCF at its promoter.

WHAT CHANGED FROM THE OLD BUILD
  - CPM tracks (25 bp bins, extendReads 200) instead of the log2-over-input tracks, whose
    single-read quantisation at 10 bp bins was what made peaks look broad.
  - A Yu 2024 E16.5 forebrain reference track, as in the main tracks panel.
  - The Bonev NPC Hi-C boundary drawn as an independent call, separate from our own.
  - Gene labels corrected; Evf2 named, since that is the biologically interesting one.
  - The percentile stated on the panel rather than only in a comment.

THE THIRD LAYER, EXPRESSION, IS NOW A SEPARATE PANEL (script 153). Summary: Dlx6os1 (Evf2) is raised in the mutant in every one of the six
cell types it could be tested in, and significantly so in two: MGE-derived interneurons
(+0.645, padj 8.0e-6) and striatal D1 SPNs (+0.458, padj 7.9e-3). The effect is confined to
ventral populations, which is why the single-nucleus data is the right instrument here and
bulk forebrain RNA would have diluted it.

  NOT FOR THE LEGEND, BUT FOR THE RESULTS TEXT: Dlx5 itself is unchanged in all four cell
  types it was tested in, and Dlx6 in both of its two. So this is not "Evf2 rises and drives
  Dlx5/Dlx6"; the canonical targets do not move. And CTCF at a promoter can activate or
  repress, so gained binding alongside a raised transcript is co-occurrence, not mechanism.

ONE LIBRARY PER GENOTYPE. No per-site significance test is possible here, and the insulation
difference carries a genome-wide offset that is removed before plotting. Nothing on this panel
is a test; it is a picture of one locus.

Runs with E18p5_clean/.hic_venv/bin/python.
Date: 2026-08-16
"""
import os, csv, gzip, collections
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel_io import save_twin
from matplotlib.patches import Patch
import pyBigWig

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
B = "/Users/alireza/Desktop/share_seq final"
TR = os.path.join(Z, "E18p5_clean/data/ctcf_reanalysis_2026/tracks")
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
HIC = os.path.join(Z, "E18p5_clean/results/discovery_hic_14/tables")
OUT = os.path.join(Z, "E18p5_clean/results/chipseq_2026/figures")
FDN = os.path.join(B, "public_3Dgenome_4DN/bonev_2017/4dn")
TSSF = os.path.join(B, "share_seq_run2_pipeline/03_downstream_analysis/results/telencephalic_only/"
                       "mechanistic_footprint_target_TAD_DORC/data/mm10_tss_all.tsv")
REF_BW = "/Users/alireza/Downloads/GSE200114_CTCF_ChIP_FB_e165_pooled.fold_change.signal.bigwig"
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "axes.linewidth": 0.6,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

# widened from 6.80-7.06 Mb so BOTH edges of the domain are in frame. The Dlx
# cluster sits inside a ~380 kb NPC domain bounded by a Weak boundary at 6.575 Mb
# and the Strong one at 6.955 Mb; the earlier window showed only the downstream edge.
CH, LO, HI = "chr6", 6_500_000, 7_100_000
CLS_COL = {"lost": "#BF443A", "shared": "#6E7B8B", "gained": "#E0A53D"}
GT_COL = {"wt": "#3B6FD4", "mut": "#7B4EA8"}
REF_COL = "#5A6068"
# labelled by hand because the automatic nearest-TSS call would print Dlx6os1 three times
GENES = [(6_866_000, "Evf2"), (6_882_400, "Dlx5"), (6_956_200, "Sdhaf3")]

bw = {g: pyBigWig.open(os.path.join(TR, f"{n}_ctcf.cpm.bs25.ext200.bw"))
      for g, n in (("wt", "wt"), ("mut", "homo"))}
refbw = pyBigWig.open(REF_BW) if os.path.exists(REF_BW) else None
NB = 1400


def prof(b):
    v = b.stats(CH, LO, HI, nBins=NB, type="mean")
    return np.array([0.0 if q is None else q for q in v], dtype=float)


x = np.linspace(LO, HI, NB) / 1e6
sig = {g: prof(b) for g, b in bw.items()}
if refbw is not None:
    sig["ref"] = prof(refbw)

ins = {}
for g, f in (("wt", "insulation_wt_25kb.bedgraph"), ("mut", "insulation_homo_25kb.bedgraph")):
    pts = []
    for line in open(os.path.join(HIC, f)):
        v = line.split()
        if v[0] != CH:
            continue
        s, e = int(v[1]), int(v[2])
        if e < LO or s > HI:
            continue
        try:
            pts.append(((s + e) / 2, float(v[3])))
        except ValueError:
            continue
    ins[g] = np.array(sorted(pts)) if pts else np.empty((0, 2))
    print(f"[152] insulation {g}: {len(ins[g])} bins in window")
# the mutant track carries a genome-wide offset; remove it so only local shape is compared
if len(ins["wt"]) and len(ins["mut"]):
    allw, allm = [], []
    for g, f in (("wt", "insulation_wt_25kb.bedgraph"), ("mut", "insulation_homo_25kb.bedgraph")):
        vals = []
        for line in open(os.path.join(HIC, f)):
            v = line.split()
            try:
                vals.append(float(v[3]))
            except (IndexError, ValueError):
                continue
        (allw if g == "wt" else allm).append(np.array(vals))
    offset = np.nanmedian(allm[0]) - np.nanmedian(allw[0])
    ins["mut"][:, 1] -= offset
    print(f"[152] removed genome-wide insulation offset of {offset:+.3f}")

pk = [r for r in csv.DictReader(
    l for l in open(os.path.join(T, "F4_peak_classes_with_position.csv"))
    if l.strip() and not l.startswith("#"))
    if r["chr"] == CH and LO < int(r["end"]) and int(r["start"]) < HI]
bnd = [(int(v[1]), int(v[2]), v[3]) for v in
       (l.split("\t") for l in gzip.open(
           os.path.join(FDN, "Bonev_NPC_boundaries_4DNFIMMHY96D.bed.gz"), "rt"))
       if v[0] == CH and LO < int(v[2]) and int(v[1]) < HI]
print(f"[152] {len(pk)} peaks, {len(bnd)} independent NPC boundaries in window")

nrow = 5 if refbw is not None else 4
hr = [1.5] + [1] * (nrow - 2) + [1.15]   # annotation row needs room for gene labels
fig, axes = plt.subplots(nrow, 1, figsize=(8.6, 5.5),
                         gridspec_kw=dict(height_ratios=hr, hspace=0.16))
fig.subplots_adjust(left=0.155, right=0.860, top=0.925, bottom=0.075)

ax = axes[0]
for g, lab in (("wt", "wild type"), ("mut", "mutant")):
    if len(ins[g]):
        ax.plot(ins[g][:, 0] / 1e6, ins[g][:, 1], color=GT_COL[g], lw=1.6, marker="o",
                ms=2.6, label=lab, zorder=3)
ax.set_ylabel("insulation", fontsize=8.2)
# the in-axes legend sat on the wild-type curve; genotype keys move to the right margin,
# above the peak-call key, so nothing overlaps the data
ax.legend(frameon=False, fontsize=7.2, loc="upper left", bbox_to_anchor=(1.01, 1.02),
          handlelength=1.2, title="insulation")
ax.get_legend().get_title().set_fontsize(7.4)
ax.get_legend().get_title().set_color("0.35")
ax.tick_params(length=2, labelbottom=False)
ttl = ax.set_title("A boundary that loses CTCF, beside a lncRNA promoter that gains it",
             fontsize=10.6, pad=8, loc="left")

for i, (g, lab, col) in enumerate((("wt", "WT", GT_COL["wt"]), ("mut", "MUT", GT_COL["mut"]),
                                   ("ref", "E16.5 ref", REF_COL))):
    if g == "ref" and refbw is None:
        continue
    ax = axes[1 + i]
    ax.fill_between(x, 0, sig[g], color=col, lw=0, zorder=3)
    ax.set_ylim(0, max(sig["wt"].max(), sig["mut"].max()) * 1.1 if g != "ref"
                else sig["ref"].max() * 1.1)
    ax.set_ylabel(lab, fontsize=7.6, color=col, rotation=0, ha="right", va="center",
                  labelpad=12)
    ax.set_yticks([]); ax.tick_params(length=2, labelbottom=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    unit = "fold" if g == "ref" else "CPM"
    ax.text(0.997, 0.86, f"[0 – {ax.get_ylim()[1]/1.1:.0f}] {unit}", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.4, color="0.45")

ax = axes[-1]
ax.set_ylim(0, 1)
for sp in ("top", "right", "left"):
    ax.spines[sp].set_visible(False)
ax.set_yticks([])
# everything below lives in ONE annotation axis, in three stacked bands, so nothing can
# spill into the coverage tracks the way the gene labels did in the previous build
for s_, e_, cl in ((int(r["start"]), int(r["end"]), r["cls"]) for r in pk):
    ax.add_patch(plt.Rectangle((s_ / 1e6, 0.80), max((e_ - s_) / 1e6, 0.0016), 0.16,
                               facecolor=CLS_COL[cl], edgecolor="none", clip_on=False))
for s_, e_, kind in bnd:
    ax.add_patch(plt.Rectangle((s_ / 1e6, 0.52), (e_ - s_) / 1e6, 0.16,
                               facecolor="#2A2A2E", alpha=0.75, edgecolor="none"))
    ax.text((s_ + e_) / 2e6, 0.48, f"{kind.strip()} boundary", ha="center",
            va="top", fontsize=6.0, color="#2A2A2E")
# gene labels alternate height so the three TSS within 20 kb do not overprint
for i, (p_, g) in enumerate(GENES):
    y = 0.26 if i % 2 == 0 else 0.02
    ax.plot([p_ / 1e6], [0.34], marker="v", ms=4.2, color="0.35")
    ax.plot([p_ / 1e6, p_ / 1e6], [y + 0.07, 0.33], color="0.7", lw=0.6, zorder=1)
    ax.text(p_ / 1e6, y, f"$\\it{{{g}}}$", ha="center", va="center", fontsize=7.0,
            color="0.25")
ax.text(-0.010, 0.88, "CTCF peaks", transform=ax.transAxes, ha="right", va="center",
        fontsize=6.8, color="0.45")
ax.text(-0.010, 0.60, "Bonev NPC Hi-C", transform=ax.transAxes, ha="right", va="center",
        fontsize=6.8, color="0.45")
ax.tick_params(length=2, labelbottom=True)
ax.set_xlabel(f"{CH} (Mb)", fontsize=8.6)
for a in axes:
    a.set_xlim(LO / 1e6, HI / 1e6)

fig.legend(handles=[Patch(facecolor=CLS_COL[c], label=c) for c in ("lost", "shared", "gained")],
           loc="upper left", bbox_to_anchor=(0.870, 0.62), frameon=False, fontsize=7.2,
           title="peak call")
# NO EXPLANATORY CAPTION. The gene identities, the domain, the independent-dataset
# support and the 99th-percentile caveat are all written into the Results text, not
# onto the panel. What stays is only what is needed to READ the figure: the axis
# labels, the track names and the two colour keys. Everything the panel used to say
# in grey type is preserved in this script's docstring so it is not lost.

# the explanatory caption was already removed from this panel, so only the title differs
save_twin(fig, OUT, "F4G_dlx_locus", captions=[ttl],
          bare_size=(8.6, 5.2), bare_adjust=dict(top=0.965, bottom=0.075))

with open(os.path.join(T, "F4G_dlx_locus_source.csv"), "w", newline="") as f:
    f.write("# Dlx5/Dlx6 locus panel. Illustration, NOT evidence of the typical effect:\n")
    f.write("#   insulation moves 0.577 units here against an average of 0.043, ~99th pct.\n")
    f.write("# Gene assignment corrected: the GAINED promoter peak is 0.5 kb from Dlx6os1\n")
    f.write("#   (Evf2); the LOST promoter peak is on Sdhaf3, with Dlx5 73.9 kb away. So the\n")
    f.write("#   loss is at the domain boundary, not at Dlx5/Dlx6 itself.\n")
    f.write("# Bonev NPC Hi-C independently calls a Strong boundary at 6.955 Mb.\n")
    w = csv.DictWriter(f, fieldnames=["chr", "start", "end", "cls", "position"])
    w.writeheader()
    w.writerows([{k: r[k] for k in ("chr", "start", "end", "cls", "position")} for r in pk])
print("[152] wrote F4G_dlx_locus.pdf and F4G_dlx_locus_source.csv")
