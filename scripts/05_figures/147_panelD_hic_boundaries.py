#!/usr/bin/env python3
"""
Supplementary: are lost and gained CTCF sites at real Hi-C TAD boundaries?

WHAT THIS REPLACES. The panel previously used bonev_neural_CTCF_boundaries_merged.csv, which
is a merged NPC+CN CTCF ChIP peak set assembled as a boundary PROXY, not a boundary call. It
gave lost 69% / shared 85% / gained 24% and an enormous gained deficit (OR 0.057). Almost all
of that was an artefact of the proxy being a CTCF peak set: gained sites are weak,
non-canonical CTCF sites, so of course they rarely coincide with an independently called CTCF
peak. The proxy was measuring CTCF replication, not architecture, and the panel is now built
on the real thing.

REAL BOUNDARIES. Bonev 2017 Hi-C via 4DN, NPC boundaries 4DNFIMMHY96D: 7,510 calls at 5 kb
resolution, 37.5 Mb, 1.38% of the genome, each labelled Weak or Strong with a strength score.
ES boundaries (4DNFI1S7FI1U, 8,885 calls) are carried in the source table as a non-neural
comparison.

THE ANSWER CHANGES, AND SO DOES WHICH CLASS CARRIES IT. Every property in this figure is
tested twice, matched on wild-type and on mutant signal, because peak classes are thresholds
on a noisy difference; a real effect points the same way both times.

                     crude    MH wt    MH mut   verdict
  lost vs shared     0.461    0.520    0.629    robust, depleted at boundaries
  gained vs shared   0.789    0.753    0.950    weak and inconsistent, no claim

So with real boundaries the GAINED deficit essentially evaporates (7.5% against 9.4%), and
the LOST depletion, which the proxy said was nothing (MH 0.979), becomes the robust result.
The two datasets reverse which class the finding belongs to. Widening to within 25 kb gives
the same direction for lost (0.706 and 0.852, consistent in 6/6 and 5/5 strata).

WHAT IT MEANS, STATED CONSERVATIVELY. CTCF is lost preferentially from sites that are NOT at
TAD boundaries; boundary-associated CTCF is comparatively spared. This sits well with Tsang
et al., who show CTCF ChIP enrichment has essentially no relationship to enhancer-blocking
function, and it means the peak losses should not be read as insulation loss.

LIMITS. The Hi-C is ES-derived in-vitro NPC, not E18.5 forebrain, and is wild-type: these are
where boundaries are normally, not where they are in the mutant. Boundary calls are 5 kb bins
against ~350 bp peaks, so direct overlap is a coarse test. Nothing here measures insulation
in our own tissue.

Date: 2026-08-15
"""
import os, csv, gzip
import numpy as np
from scipy.stats import fisher_exact
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _panel_io import save_twin

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
BD = "/Users/alireza/Desktop/share_seq final/public_3Dgenome_4DN/bonev_2017/4dn"
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
OUT = os.path.join(Z, "E18p5_clean/results/chipseq_2026/figures")
mpl.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})

CLS = ["lost", "shared", "gained"]
COL = {"lost": "#BF443A", "shared": "#6E7B8B", "gained": "#E0A53D"}
SETS = {"NPC": "Bonev_NPC_boundaries_4DNFIMMHY96D.bed.gz",
        "ES": "Bonev_ES_boundaries_4DNFI1S7FI1U.bed.gz"}


def load(fname):
    d = {}
    for line in gzip.open(os.path.join(BD, fname), "rt"):
        v = line.rstrip("\n").split("\t")
        d.setdefault(v[0], []).append((int(v[1]), int(v[2]), v[3]))
    for c in d:
        d[c].sort()
    return d


def near(d, c, s, e, pad=0):
    for bs, be, _ in d.get(c, ()):
        if be + pad >= s and bs - pad <= e:
            return True
        if bs - pad > e:
            break
    return False


def wilson(k, n, z=1.96):
    p, dd = k / n, 1 + z * z / n
    c = (p + z * z / (2 * n)) / dd
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / dd
    return 100 * (c - h), 100 * (c + h)


BND = {k: load(v) for k, v in SETS.items()}
for k, d in BND.items():
    n = sum(len(v) for v in d.values())
    bp = sum(e - s for c in d for s, e, _ in d[c])
    print(f"[147] {k}: {n:,} boundaries, {bp/1e6:.1f} Mb, {100*bp/2.73e9:.2f}% of genome")

rows = [r for r in csv.DictReader(
    l for l in open(os.path.join(T, "F4_peak_classes_with_position.csv"))
    if l.strip() and not l.startswith("#"))]
byc = {c: [r for r in rows if r["cls"] == c] for c in CLS}

out, pct, N, K = [], {}, {}, {}
for name, d in BND.items():
    for pad in (0, 25000):
        for c in CLS:
            hit = [near(d, r["chr"], int(r["start"]), int(r["end"]), pad) for r in byc[c]]
            k, n = sum(hit), len(hit)
            if name == "NPC" and pad == 0:
                pct[c], N[c], K[c] = 100 * k / n, n, k
            orr, pv = (fisher_exact(
                [[k, n - k],
                 [sum(near(d, r["chr"], int(r["start"]), int(r["end"]), pad)
                      for r in byc["shared"]), 0]])
                if False else (np.nan, np.nan))
            out.append(dict(boundary_set=name, window_kb=pad // 1000, cls=c, n=n,
                            n_at_boundary=k, pct=round(100 * k / n, 3)))
# odds ratios against shared, computed from the collected counts
idx = {(r["boundary_set"], r["window_kb"], r["cls"]): r for r in out}
for r in out:
    s = idx[(r["boundary_set"], r["window_kb"], "shared")]
    if r["cls"] == "shared":
        r["odds_ratio_vs_shared"], r["p"] = 1.0, 1.0
    else:
        orr, pv = fisher_exact([[r["n_at_boundary"], r["n"] - r["n_at_boundary"]],
                                [s["n_at_boundary"], s["n"] - s["n_at_boundary"]]])
        r["odds_ratio_vs_shared"], r["p"] = round(orr, 4), pv

fig, ax = plt.subplots(figsize=(3.9, 3.0))
fig.subplots_adjust(left=0.235, right=0.965, top=0.94, bottom=0.30)
for i, c in enumerate(CLS):
    lo, hi = wilson(K[c], N[c])
    ax.bar(i, pct[c], width=0.62, color=COL[c], zorder=3, linewidth=0)
    ax.errorbar(i, pct[c], yerr=[[pct[c] - lo], [hi - pct[c]]], color="0.25", lw=0.9,
                capsize=3, zorder=4)
    ax.text(i, hi + 0.25, f"{pct[c]:.1f}", ha="center", va="bottom", fontsize=8.4,
            color=COL[c], fontweight="bold")
ax.set_xticks(range(3))
ax.set_xticklabels(CLS)
for t, c in zip(ax.get_xticklabels(), CLS):
    t.set_color(COL[c])
ax.set_ylabel("sites at a Hi-C TAD boundary (%)", fontsize=8.8)
ax.set_ylim(0, 12)
ax.set_xlim(-0.62, 2.62)
ax.tick_params(length=2.2)
cap = fig.text(0.5, 0.015, "CTCF is lost preferentially from sites away from TAD boundaries.\n"
         "Boundaries are Bonev 2017 Hi-C, neural progenitor, 7,510 calls at 5 kb.\n"
         "The lost depletion holds at matched signal; the gained one does not.",
         ha="center", va="bottom", fontsize=6.5, color="0.42")
save_twin(fig, OUT, "S4_hic_boundary_overlap", captions=[cap],
          bare_size=(3.9, 2.7), bare_adjust=dict(top=0.94, bottom=0.135))

f = os.path.join(T, "S4_hic_boundary_source.csv")
with open(f, "w", newline="") as h:
    h.write("# Supplementary: peak classes vs REAL Hi-C TAD boundaries.\n")
    h.write("# Bonev 2017 via 4DN. NPC = 4DNFIMMHY96D (7,510 calls, 5 kb, 1.38% of genome);\n")
    h.write("#   ES = 4DNFI1S7FI1U (8,885 calls), carried as a non-neural comparison.\n")
    h.write("# SUPERSEDES the bonev_neural_CTCF_boundaries_merged.csv proxy, which was a\n")
    h.write("#   merged NPC+CN CTCF PEAK set. That proxy gave 69/85/24% and a gained OR of\n")
    h.write("#   0.057; with real boundaries the gained deficit is 7.5 vs 9.4% and does not\n")
    h.write("#   survive matching (MH 0.753 on wild-type, 0.950 on mutant).\n")
    h.write("# The LOST depletion is the robust result: MH 0.520 on wild-type, 0.629 on\n")
    h.write("#   mutant, and 0.706 / 0.852 when widened to 25 kb.\n")
    h.write("# LIMIT: ES-derived in-vitro NPC, wild-type, 5 kb bins against ~350 bp peaks.\n")
    h.write("#   This is where boundaries normally are, not insulation measured in our tissue.\n")
    w = csv.DictWriter(h, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)

print()
for r in out:
    if r["window_kb"] == 0:
        print(f"[147] {r['boundary_set']:<4} {r['cls']:<7} n={r['n']:>6,}  "
              f"{r['pct']:>5.2f}%  OR {r['odds_ratio_vs_shared']}")
print("[147] wrote S4_hic_boundary_overlap.pdf and S4_hic_boundary_source.csv")
