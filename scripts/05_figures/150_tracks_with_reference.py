#!/usr/bin/env python3
"""
Tracks panel plus an independent forebrain CTCF reference. Two variants, to compare.

THE PROBLEM THIS ADDRESSES. Our changed peaks are the weak tail of the CTCF landscape: median
peak-genotype CPM 1.40 lost and 1.22 gained against 3.49 shared, and none of the sampled
changed sites reached 10 CPM. A reader looking at a 1.5 CPM trace is entitled to ask whether
the peak is real at all. An independent wild-type dataset answers that, and the answer is
reassuring for the losses and interesting for the gains:

  % of our peaks also called in each reference    lost   shared   gained
    Yu forebrain E12.5                           74.9%    92.5%    40.3%
    Yu forebrain E13.5                           77.2%    94.8%    46.7%
    Yu forebrain E14.5                           92.0%    99.2%    71.5%
    Yu forebrain E15.5                           86.3%    98.3%    69.6%
    Yu forebrain E16.5                           81.0%    96.7%    51.5%
    Yu forebrain P0                              56.5%    88.0%    22.9%
    Bonev NPC (ES-derived, in vitro)             61.8%    77.0%    19.7%
    Bonev CN  (ES-derived, in vitro)             52.5%    77.8%     8.7%

Shared > lost > gained in every one of the eight, including P0 which brackets our E18.5 from
the other side. So the lost sites are bona fide forebrain CTCF sites that lose signal, and the
gained sites are substantially non-canonical.

WHY YU AND NOT BONEV. Yu et al. 2024 (GSE200114) is mouse FOREBRAIN at E12.5, E13.5, E14.5,
E15.5, E16.5 and P0: same tissue, six stages bracketing our E18.5. Bonev is ES-derived
in-vitro NPC and cortical neurons. Bonev is kept in the table as the weaker comparison, not
as the reference.

READ THE COLUMNS, NOT THE ROWS. The reference sets differ in size, 43,838 to 59,365 peaks, so
E14.5 scoring highest is partly sequencing depth and peak-calling sensitivity. Comparisons
across our three classes WITHIN one reference are the valid ones, and those are consistent
everywhere.

WHICH LOCI SURVIVED, AND WHY THE OTHERS WENT. The reference decided one of these outright.

  kept    Hopx    lost,   gliogenic / radial glia   called at all 6 Yu stages, ref fold 18.4
  kept    Foxp2   lost,   dorsal deep layer         all 6 stages, ref fold 26.5, best support
  kept    Pou3f3  gained, dorsal upper layer        the ONLY true promoter example, +0.1 kb
  kept    Sp9     gained, LGE-derived               all 6 stages, ref fold 24.3
  kept    Foxp1   gained, SPN                       5 of 6 stages, ref fold 6.3

  DROPPED Sox6    lost,   MGE interneuron   NOT CALLED at any of the six Yu forebrain stages,
          and its reference track is flat noise at 3.3 fold with nothing under our site. Our
          own call rests on 1.51 CPM. It was the only qualifying MGE candidate, so the panel
          now has no MGE example rather than an unsupported one.
  DROPPED Sox2    gained, cycling progenitor   well supported, but the site is 21.2 kb from
          the TSS, so the attribution to the gene is much weaker than the others.
  DROPPED Nr2f1   gained, CGE interneuron      supported at 5 of 6 stages but 9.9 kb out, and
          the chapter's snRNA-seq scope is cortical and striatal rather than interneuron.

The surviving five cover dorsal (Hopx, Foxp2, Pou3f3) and ventral striatal (Sp9, Foxp1),
which is the scope of the accompanying single-nucleus analysis, and split 2 lost / 3 gained.

THE TWO VARIANTS BUILT HERE
  A  "ticks": a row of six marks per locus, filled where that stage called a peak. Costs
     nothing, uses only the peak calls already on disk, and shows the developmental profile.
  B  "track": a third coverage row from the Yu E16.5 pooled fold-change bigWig, the closest
     stage to E18.5 that ships a signal track (P0 ships peaks only). More vivid, and shows
     peak shape, but it is ONE stage and loses the developmental axis.

UNITS DO NOT MATCH AND MUST NOT SHARE AN AXIS. Our tracks are CPM; the Yu bigWig is pooled
fold change over control, genome max 110. The reference row therefore carries its own
autoscale and its own unit label. It authenticates the SITE. It says nothing about our
genotype comparison: different study, different animals, E16.5 not E18.5, wild type only.

Runs with E18p5_clean/.hic_venv/bin/python.
Date: 2026-08-16
"""
import os, csv, collections, gzip
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pyBigWig

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
B = "/Users/alireza/Desktop/share_seq final"
TR = os.path.join(Z, "E18p5_clean/data/ctcf_reanalysis_2026/tracks")
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
OUT = os.path.join(Z, "E18p5_clean/results/chipseq_2026/figures")
YU = os.path.join(B, "public_3Dgenome_4DN/yu_2024_forebrain")
REF_BW = "/Users/alireza/Downloads/GSE200114_CTCF_ChIP_FB_e165_pooled.fold_change.signal.bigwig"
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "axes.linewidth": 0.6})

FLANK, BIN, PAD = 2000, 25, 1.18
CLS_COL = {"lost": "#BF443A", "shared": "#6E7B8B", "gained": "#E0A53D"}
GT_COL = {"wt": "#3B6FD4", "mut": "#7B4EA8"}
REF_COL = "#5A6068"
STAGES = [("E12.5", "GSE200114_CTCF_ChIP_FB_e125.replicated_peaks.narrowPeak.gz"),
          ("E13.5", "GSE200114_CTCF_ChIP_FB_e135.replicated_peaks.narrowPeak.gz"),
          ("E14.5", "GSE200114_CTCF_ChIP_FB_e145.replicated_peaks.narrowPeak.gz"),
          ("E15.5", "GSE200114_CTCF_ChIP_FB_e155.replicated_peaks.narrowPeak.gz"),
          ("E16.5", "GSE200114_CTCF_ChIP_FB_e165.replicated_peaks.narrowPeak.gz"),
          ("P0",    "GSE200114_P0_FB_CTCT_ChIP_pseudo_replicated_peaks.narrowPeak.gz")]

# loci come from the table script 148 writes, so the two scripts cannot drift apart
KEEP = ["REF", "Hopx", "Foxp2", "Pou3f3", "Sp9", "Foxp1"]
all_loci = [r for r in csv.DictReader(
    l for l in open(os.path.join(T, "F4E_tracks_source.csv"))
    if l.strip() and not l.startswith("#"))]
loci = sorted([r for r in all_loci if r["gene"] in KEEP], key=lambda r: KEEP.index(r["gene"]))
dropped = [r["gene"] for r in all_loci if r["gene"] not in KEEP]
print(f"[150] kept {len(loci)}: {', '.join(r['gene'] for r in loci)}")
print(f"[150] dropped {len(dropped)}: {', '.join(dropped)}")

ref_pk = {}
for name, fn in STAGES:
    d = collections.defaultdict(list)
    for line in gzip.open(os.path.join(YU, fn), "rt"):
        v = line.split("\t")
        d[v[0]].append((int(v[1]), int(v[2])))
    for c in d:
        d[c].sort()
    ref_pk[name] = d
    print(f"[150] {name}: {sum(len(v) for v in d.values()):,} peaks")

pk = collections.defaultdict(list)
for r in csv.DictReader(l for l in open(os.path.join(T, "F4_peak_classes_with_position.csv"))
                        if l.strip() and not l.startswith("#")):
    pk[r["chr"]].append((int(r["start"]), int(r["end"]), r["cls"]))
for c in pk:
    pk[c].sort()
bw = {g: pyBigWig.open(os.path.join(TR, f"{n}_ctcf.cpm.bs25.ext200.bw"))
      for g, n in (("wt", "wt"), ("mut", "homo"))}
refbw = pyBigWig.open(REF_BW) if os.path.exists(REF_BW) else None
if refbw is None:
    print(f"[150] NOTE reference bigWig not found, variant B skipped:\n      {REF_BW}")


def called(stage, ch, s, e):
    for a, b in ref_pk[stage].get(ch, ()):
        if b >= s and a <= e:
            return True
        if a > e:
            break
    return False


nb = 2 * FLANK // BIN
x = (np.arange(nb) * BIN) - FLANK + BIN / 2


def build(mode):
    """mode 'ticks' = reference presence row; mode 'track' = reference coverage row."""
    ntrack = 3 if mode == "track" else 2
    hr = [1] * ntrack + [0.17] + ([0.30] if mode == "ticks" else [])
    fig = plt.figure(figsize=(9.8, 5.5 if mode == "track" else 5.3))
    outer = fig.add_gridspec(2, 3, hspace=0.68, wspace=0.20,
                             left=0.075, right=0.988, top=0.805, bottom=0.185)
    for k, r in enumerate(loci):
        gene, ch, ctr = r["gene"], r["chr"], int(r["centre"])
        cls, is_ref = r["cls"], r["gene"] == "REF"
        inner = outer[k // 3, k % 3].subgridspec(len(hr), 1, height_ratios=hr, hspace=0.10)
        sig = {}
        for g, b in bw.items():
            v = b.stats(ch, ctr - FLANK, ctr + FLANK, nBins=nb, type="mean")
            sig[g] = np.array([0.0 if q is None else q for q in v], dtype=float)
        ymax = max(sig["wt"].max(), sig["mut"].max()) * PAD or 1.0
        order = [("wt", ymax), ("mut", ymax)]
        if mode == "track":
            v = refbw.stats(ch, ctr - FLANK, ctr + FLANK, nBins=nb, type="mean")
            sig["ref"] = np.array([0.0 if q is None else q for q in v], dtype=float)
            order.append(("ref", (sig["ref"].max() * PAD) or 1.0))
        for j, (g, top) in enumerate(order):
            ax = fig.add_subplot(inner[j])
            col = REF_COL if g == "ref" else GT_COL[g]
            ax.fill_between(x, 0, sig[g], color=col, lw=0, zorder=3)
            ax.set_xlim(-FLANK, FLANK); ax.set_ylim(0, top)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ("top", "right", "bottom"):
                ax.spines[sp].set_visible(False)
            lab = {"wt": "WT", "mut": "MUT", "ref": "E16.5"}[g]
            ax.text(-0.045, 0.5, lab, transform=ax.transAxes, ha="right", va="center",
                    fontsize=6.4, color=col)
            if g == "ref":
                ax.text(0.985, 0.88, f"[0 – {top/PAD:.0f}] fold", transform=ax.transAxes,
                        ha="right", va="top", fontsize=6.6, color="0.45")
            if j == 0:
                ax.text(0.985, 0.88, f"[0 – {top/PAD:.1f}] CPM", transform=ax.transAxes,
                        ha="right", va="top", fontsize=7.2,
                        color="#2A2A2E" if is_ref else "0.35",
                        fontweight="bold" if is_ref else "normal")
                ax.set_title(f"$\\it{{{gene}}}$" if not is_ref else "stable site",
                             fontsize=9.4, pad=7)
                if not is_ref:
                    w, m = float(r["wt_cpm"]), float(r["mut_cpm"])
                    fold = (m / max(w, 0.05)) if cls == "gained" else (w / max(m, 0.05))
                    ax.text(0.015, 0.86, f"{cls}  {fold:.1f}×", transform=ax.transAxes,
                            ha="left", va="top", fontsize=7.0, color=CLS_COL[cls],
                            fontweight="bold")
                    ax.text(0.015, 0.40, f"{int(r['dist_to_tss'])/1000:+.1f} kb from TSS",
                            transform=ax.transAxes, ha="left", va="top", fontsize=5.9,
                            color="0.45")
        axp = fig.add_subplot(inner[ntrack])
        axp.set_xlim(-FLANK, FLANK); axp.set_ylim(0, 1); axp.axis("off")
        for s, e, cl in pk.get(ch, []):
            if e < ctr - FLANK or s > ctr + FLANK:
                continue
            axp.add_patch(plt.Rectangle((s - ctr, 0.25), max(e - s, 60), 0.5,
                                        facecolor=CLS_COL[cl], edgecolor="none"))
        if mode == "ticks":
            axt = fig.add_subplot(inner[ntrack + 1])
            axt.set_xlim(-0.6, len(STAGES) - 0.4); axt.set_ylim(0, 1); axt.axis("off")
            s0, e0 = ctr - 150, ctr + 150
            for i, (name, _f) in enumerate(STAGES):
                on = called(name, ch, s0, e0)
                axt.add_patch(plt.Rectangle((i - 0.34, 0.30), 0.68, 0.46,
                                            facecolor=REF_COL if on else "none",
                                            edgecolor=REF_COL, linewidth=0.7))
                if k == 0:
                    axt.text(i, -0.30, name.replace("E", "").replace(".5", "5"),
                             ha="center", va="top", fontsize=5.2, color="0.45")
            axt.text(-0.75, 0.53, "ref", ha="right", va="center", fontsize=5.8,
                     color="0.45", transform=axt.transData)
        else:
            axp.text(0, -0.55, "±2 kb", transform=axp.transAxes, ha="left", va="top",
                     fontsize=6.2, color="0.5")

    fig.suptitle("CTCF coverage at changed sites near identity genes",
                 fontsize=11.5, x=0.075, ha="left", y=0.988)
    extra = ("Bottom row of each panel: whether independent Yu 2024 forebrain CTCF ChIP "
             "called a peak here, at E12.5 to P0. Filled = called."
             if mode == "ticks" else
             "Third track is Yu 2024 forebrain CTCF ChIP at E16.5, pooled fold change over "
             "control. Separate unit, separate scale: it authenticates the site, not the "
             "genotype difference.")
    sub_txt = fig.text(0.075, 0.912,
             "CPM coverage, 25 bp bins, reads extended to 200 bp. IGV group autoscale, shared "
             "between genotypes within a locus. First cell is a stable shared site, for "
             f"scale.\n{extra}",
             ha="left", va="bottom", fontsize=6.9, color="0.42")
    foot_txt = fig.text(0.075, 0.020,
             "Independent replication genome-wide: at E14.5, 92.0% of lost and 99.2% of "
             "shared sites are called by Yu, against 71.5% of gained; at P0 the ordering holds "
             "at 56.5, 88.0 and 22.9%. The losses are\nbona fide forebrain CTCF sites and the "
             "gains are substantially non-canonical. Every locus here is called at 5 or 6 of "
             "the 6 Yu stages. Sox6 was dropped for failing at all six, Sox2 and Nr2f1 for "
             "weak\ngene attribution at 21.2 and 9.9 kb. Reference sets differ in size, so "
             "compare classes within a stage, not stages with each other.",
             ha="left", va="bottom", fontsize=6.4, color="0.42")
    leg = fig.legend(handles=[Patch(facecolor=CLS_COL[c], label=c)
                              for c in ("lost", "shared", "gained")],
                     loc="lower center", bbox_to_anchor=(0.5, 0.062), ncol=3, frameon=False,
                     fontsize=7.4, title="peak call, bar under each track")
    leg.get_title().set_fontsize(6.8); leg.get_title().set_color("0.45")
    stem = "F4E_tracks" if mode == "track" else "F4E_tracks_ticks_variant"
    f = os.path.join(OUT, stem)
    fig.savefig(f + ".pdf"); fig.savefig(f + ".png", dpi=210)
    print(f"[150] wrote {stem}.pdf")

    if mode == "track":
        # caption-free twin for assembly: same figure object with the text artists hidden and
        # the canvas shrunk, so the two can never drift apart. The peak-call key stays, since
        # it is needed to read the panel.
        for artist in (fig._suptitle, sub_txt, foot_txt):
            artist.set_visible(False)
        fig.set_size_inches(9.8, 4.5)
        fig.subplots_adjust(top=0.945, bottom=0.115)
        leg.set_bbox_to_anchor((0.5, 0.010))
        fig.savefig(f + "_bare.pdf"); fig.savefig(f + "_bare.png", dpi=210)
        print(f"[150] wrote {stem}_bare.pdf")
    plt.close(fig)


# the coverage variant was chosen over the tick row: it shows peak SHAPE, so the
# reference peak can be seen sitting under ours, and it exposed Sox6 as unsupported
# in a way six empty boxes did not. The tick builder is kept and still runs.
if refbw is not None:
    build("track")
# the tick-row variant lost the comparison and is not rebuilt; see the docstring.
# Re-enable by calling build("ticks") if the reference bigWig is ever unavailable.

rows = []
for r in loci:
    ch, ctr = r["chr"], int(r["centre"])
    d = dict(gene=r["gene"], cls=r["cls"], chr=ch, centre=ctr)
    for name, _f in STAGES:
        d[f"yu_{name}"] = called(name, ch, ctr - 150, ctr + 150)
    if refbw is not None:
        v = refbw.stats(ch, ctr - 150, ctr + 150, type="max")
        d["yu_E16.5_fold"] = round(float(v[0]), 2) if v and v[0] is not None else ""
    rows.append(d)
with open(os.path.join(T, "F4E_tracks_reference_support.csv"), "w", newline="") as h:
    h.write("# Independent support for each locus in the tracks panel.\n")
    h.write("# Yu et al. 2024 GSE200114, mouse FOREBRAIN CTCF ChIP, E12.5 to P0. Same tissue,\n")
    h.write("#   six stages bracketing our E18.5. Called = a Yu peak overlaps +/-150 bp of\n")
    h.write("#   our site centre. yu_E16.5_fold is the pooled fold-change bigWig max.\n")
    h.write("# Genome-wide, at E14.5: 92.0% of lost, 99.2% of shared, 71.5% of gained are\n")
    h.write("#   called by Yu; at P0, 56.5 / 88.0 / 22.9%. Ordering holds in all 6 stages and\n")
    h.write("#   in both Bonev sets. Reference sets differ in size (43,838 to 59,365 peaks),\n")
    h.write("#   so compare classes within a stage, never stages with each other.\n")
    w = csv.DictWriter(h, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("[150] wrote F4E_tracks_reference_support.csv")
for r in rows:
    print(f"[150] {r['gene']:<8}{r['cls']:<8}" +
          "".join("Y" if r[f"yu_{n}"] else "." for n, _ in STAGES) +
          (f"   fold {r.get('yu_E16.5_fold','')}" if refbw else ""))
