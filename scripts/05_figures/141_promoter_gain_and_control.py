#!/usr/bin/env python3
"""
Is the CTCF gain at promoters real, or is it ChIP background at open chromatin?

THE PROBLEM THIS SETTLES. Gained peak calls are promoter-enriched (32.2% against 17.9% for
shared). That comparison is confounded twice over. Peak classes are thresholds on a noisy
difference, so conditioning on a class conditions on a collider; and the classes sit in very
different places on the signal axis, where site properties vary steeply for reasons that have
nothing to do with genotype. Matching on wild-type signal and matching on mutant signal gave
opposite answers (OR 0.51 against 1.92), because both matchings are biased by regression to
the mean in opposite directions: wild-type signal is measured with noise, so within a bin the
weak sites that landed there by positive noise fall back on remeasurement.

TEST A, THE IDENTIFIED VERSION. Drop peak classes entirely and use the continuous outcome,
mutant minus wild-type signal, across all sites. Match on CORE MOTIF SCORE, which is read off
the genome and carries no ChIP noise, so it cannot drive regression to the mean. There is no
selection on the outcome, so there is no collider. This asks the question the peak-call
framing was trying to ask, and answers it cleanly.

TEST B, THE CONTROL THAT DECIDES WHETHER TEST A MEANS ANYTHING. Promoters are open chromatin,
and open chromatin yields background in any ChIP regardless of the antibody. If the mutant
library carries more background, promoters would appear to gain CTCF with no CTCF biology
involved. So the same change is measured at promoters that contain NO called CTCF peak in
either genotype. If those rise too, the effect is background. If they are flat while
CTCF-containing promoters rise, the effect is CTCF.

Non-promoter CTCF peaks are measured alongside as a third anchor.

MEASUREMENT. CPM bigWigs, 25 bp bins, reads extended to 200 bp, built on depth-matched BAMs.
Read one interval at a time through pyBigWig: rtracklayer's import(which=) returns intervals
in the bigWig's contig order rather than query order, which silently corrupted an earlier
version of this figure. Never batch the query.

Runs with E18p5_clean/.hic_venv/bin/python.
Date: 2026-08-15
"""
import os, csv, random
import numpy as np
from scipy.stats import mannwhitneyu
import pyBigWig

B = "/Users/alireza/Desktop/share_seq final"
Z = os.path.join(B, "Zhang et al 2024 Nat Comm sc-RNA-seq Analysis")
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
D = os.path.join(Z, "E18p5_clean/data/ctcf_reanalysis_2026")
TRK = os.path.join(D, "tracks")
TSSF = os.path.join(B, "share_seq_run2_pipeline/03_downstream_analysis/results/telencephalic_only/"
                       "mechanistic_footprint_target_TAD_DORC/data/mm10_tss_all.tsv")
PSEUDO = 0.05          # CPM units, well below the promoter background level
NCTRL = 6000
SEED = 42
EDG = [7, 12, 14, 16, 18, 28]


def rd(name):
    return csv.DictReader(l for l in open(os.path.join(T, name))
                          if l.strip() and not l.startswith("#"))


# ---------------------------------------------------------------- TEST A: identified version
S, C, pos = {}, {}, {}
for r in rd("F4G_per_site_signal.csv"):
    try:
        w, m = float(r["wt_centre"]), float(r["mut_centre"])
    except (TypeError, ValueError):
        continue
    if not (np.isnan(w) or np.isnan(m)):
        S[(r["chr"], r["start"], r["end"])] = (w, m, r["cls"])
for r in rd("F4C_umotif_panel_assignments.csv"):
    try:
        C[(r["chr"], r["start"], r["end"])] = float(r["core_score"])
    except (TypeError, ValueError):
        pass
for r in rd("F4_peak_classes_with_position.csv"):
    pos[(r["chr"], r["start"], r["end"])] = r
K = [k for k in S if k in C and k in pos]
print(f"[141] TEST A on {len(K):,} sites with signal and a scored core motif\n")
print("TEST A  change in CTCF signal at promoters vs elsewhere, matched on core motif score")
print("        no peak classes used, so no selection on the outcome")
print(f"{'core score':>12}{'n prom':>9}{'d prom':>9}{'n other':>10}{'d other':>9}{'diff':>8}{'p':>10}")
rowsA, allp, allo = [], [], []
for i in range(len(EDG) - 1):
    lo, hi = EDG[i], EDG[i + 1]
    sub = [k for k in K if lo < C[k] <= hi]
    p = np.array([S[k][1] - S[k][0] for k in sub if pos[k]["position"] == "Promoter"])
    o = np.array([S[k][1] - S[k][0] for k in sub if pos[k]["position"] != "Promoter"])
    if len(p) < 25 or len(o) < 25:
        continue
    _, pv = mannwhitneyu(p, o)
    allp.append(p); allo.append(o)
    print(f"{f'{lo}-{hi}':>12}{len(p):>9,}{p.mean():>9.3f}{len(o):>10,}{o.mean():>9.3f}"
          f"{p.mean()-o.mean():>8.3f}{pv:>10.1e}")
    rowsA.append(dict(core_low=lo, core_high=hi, n_promoter=len(p),
                      mean_change_promoter=round(float(p.mean()), 4), n_other=len(o),
                      mean_change_other=round(float(o.mean()), 4),
                      difference=round(float(p.mean() - o.mean()), 4), p=pv))
P, O = np.concatenate(allp), np.concatenate(allo)
print(f"\n  pooled: promoters {P.mean():+.3f}   other {O.mean():+.3f}   "
      f"difference {P.mean()-O.mean():+.3f}, same direction in "
      f"{sum(a.mean()>b.mean() for a,b in zip(allp,allo))}/{len(allp)} strata")

# ---------------------------------------------------------------- TEST B: background control
peaks = []
for g in ("wt", "homo"):
    f = os.path.join(D, f"{g}_PRIMARY_peaks.narrowPeak")
    for line in open(f):
        v = line.split("\t")
        if len(v) >= 3:
            peaks.append((v[0], int(v[1]), int(v[2])))
bychr = {}
for c, s, e in peaks:
    bychr.setdefault(c, []).append((s, e))
for c in bychr:
    bychr[c].sort()
print(f"\n[141] union of called peaks, both genotypes: {len(peaks):,}")


def hits_peak(c, s, e, pad=5000):
    """Any called peak within pad of the window. Linear scan per chromosome is fast enough."""
    for ps, pe in bychr.get(c, ()):
        if pe + pad >= s and ps - pad <= e:
            return True
        if ps - pad > e:
            break
    return False


prom = []
for line in open(TSSF):
    v = line.rstrip("\n").split("\t")
    if len(v) < 4:
        continue
    c, t, strand = v[0], int(v[1]), v[2]
    if not c.startswith("chr") or "_" in c:
        continue
    s, e = (t - 200, t + 2000) if strand == "-" else (t - 2000, t + 200)
    prom.append((c, max(1, s), e))
prom = sorted(set(prom))
random.Random(SEED).shuffle(prom)
print(f"[141] promoter windows from the project TSS table: {len(prom):,}")

free = []
for c, s, e in prom:
    if not hits_peak(c, s, e):
        free.append((c, s, e))
    if len(free) >= NCTRL:
        break
print(f"[141] CTCF-free promoter windows sampled: {len(free):,} "
      f"(no called peak in either genotype within 5 kb)")

bw = {g: pyBigWig.open(os.path.join(TRK, f"{n}_ctcf.cpm.bs25.ext200.bw"))
      for g, n in (("wt", "wt"), ("mut", "homo"))}
chroms = set(bw["wt"].chroms())


def cpm(g, c, s, e):
    if c not in chroms:
        return np.nan
    e = min(e, bw[g].chroms()[c])
    if e <= s:
        return np.nan
    v = bw[g].stats(c, s, e, type="mean")          # ONE interval at a time, never batched
    return np.nan if v is None or v[0] is None else float(v[0])


def measure(regions, label):
    out = []
    for c, s, e in regions:
        a, b = cpm("wt", c, s, e), cpm("mut", c, s, e)
        if not (np.isnan(a) or np.isnan(b)):
            out.append((a, b, np.log2((b + PSEUDO) / (a + PSEUDO))))
    arr = np.array(out)
    print(f"  {label:<34} n={len(arr):>6,}  wt CPM {arr[:,0].mean():6.3f}  "
          f"mut CPM {arr[:,1].mean():6.3f}  log2 change {np.median(arr[:,2]):+.3f}")
    return arr


print("\nTEST B  same measurement on CPM tracks, three region sets")
sets = {}
sets["CTCF-free promoters (control)"] = measure(free, "CTCF-free promoters (control)")
pk_prom = [(k[0], int(k[1]), int(k[2])) for k in K if pos[k]["position"] == "Promoter"]
pk_other = [(k[0], int(k[1]), int(k[2])) for k in K if pos[k]["position"] != "Promoter"]
random.Random(SEED).shuffle(pk_other)
sets["CTCF peaks at promoters"] = measure(pk_prom, "CTCF peaks at promoters")
sets["CTCF peaks elsewhere"] = measure(pk_other[:NCTRL], "CTCF peaks elsewhere")

ctrl = sets["CTCF-free promoters (control)"][:, 2]
real = sets["CTCF peaks at promoters"][:, 2]
_, pv = mannwhitneyu(real, ctrl)
print(f"\n  CTCF promoters vs CTCF-free promoters: median {np.median(real):+.3f} against "
      f"{np.median(ctrl):+.3f}, p = {pv:.1e}")
verdict = ("BACKGROUND: control rises too" if np.median(ctrl) > 0.5 * np.median(real)
           else "REAL: the rise needs a CTCF peak, the control is flat")
print(f"  VERDICT: {verdict}")

with open(os.path.join(T, "F4_promoter_gain_identified.csv"), "w", newline="") as f:
    f.write("# TEST A: change in CTCF signal at promoters vs elsewhere, matched on core motif\n")
    f.write("#   score. No peak classes are used, so nothing is conditioned on the outcome.\n")
    f.write("# Motif score is read off the genome and carries no ChIP noise, so it cannot\n")
    f.write("#   drive the regression to the mean that made the signal-matched versions\n")
    f.write("#   disagree (OR 0.51 matched on wild-type against 1.92 matched on mutant).\n")
    w = csv.DictWriter(f, fieldnames=list(rowsA[0].keys()))
    w.writeheader(); w.writerows(rowsA)

with open(os.path.join(T, "F4_promoter_gain_background_control.csv"), "w", newline="") as f:
    f.write("# TEST B: is the promoter gain CTCF, or ChIP background at open chromatin?\n")
    f.write("# CPM bigWigs, 25 bp bins, extendReads 200, depth-matched BAMs, read one\n")
    f.write("#   interval at a time through pyBigWig (never batched: rtracklayer's\n")
    f.write("#   import(which=) returns contig order, not query order).\n")
    f.write("# Control = project promoter windows with NO called peak in either genotype\n")
    f.write(f"#   within 5 kb. log2 change uses a pseudocount of {PSEUDO} CPM.\n")
    w = csv.writer(f)
    w.writerow(["region_set", "n", "mean_cpm_wt", "mean_cpm_mut",
                "median_log2_change", "q25_log2_change", "q75_log2_change"])
    for name, arr in sets.items():
        w.writerow([name, len(arr), round(float(arr[:, 0].mean()), 4),
                    round(float(arr[:, 1].mean()), 4),
                    round(float(np.median(arr[:, 2])), 4),
                    round(float(np.percentile(arr[:, 2], 25)), 4),
                    round(float(np.percentile(arr[:, 2], 75)), 4)])
print("\n[141] wrote F4_promoter_gain_identified.csv and "
      "F4_promoter_gain_background_control.csv")
