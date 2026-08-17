#!/usr/bin/env python3
"""
Which Figure 4 site properties are genotype effects, and which are CTCF signal strength?

WHY THIS EXISTS. Panel C reports position, anchor status and chromatin marks per peak class.
The gradient panel then showed that all three properties vary steeply with CTCF signal
strength in SHARED sites, which by definition did not change class. Since the classes sit in
very different places on the signal axis (median wild-type log2: gained 1.53, lost 2.77,
shared 3.68), any crude class comparison is confounded by strength. This script decides, for
each property, whether anything is left after conditioning on it.

METHOD. Mantel-Haenszel odds ratios across signal strata, each class against SHARED, never
lost against gained. Strata need >= 25 sites on both sides or they are dropped.

THE COLLIDER PROBLEM, AND THE TEST FOR IT. A site is called gained when mutant minus wild-type
exceeds a threshold. Conditioning on wild-type signal AND on class therefore conditions
implicitly on mutant signal, which is a collider. A stratified estimate can be manufactured
entirely by that selection. The test is to repeat the whole analysis matching on MUTANT signal
instead: a real property effect should point the same way on either axis, while a
selection artefact can invert.

That test is not decoration. It overturned the headline. Matched on wild-type, gained sites
look promoter-DEPLETED (MH 0.51, all 5 strata); matched on mutant, the same sites look
promoter-ENRICHED (MH 1.92, all 6 strata); crude says enriched (32.5% vs 18.0%). Three
different answers from the same data, so the promoter status of gained sites is NOT
IDENTIFIED by this design and cannot be claimed in either direction.

The mechanical reason is overlap. Gained sites have median wild-type signal 1.53 against 3.68
for shared, so the wild-type-matched comparison leans on strata where shared sites are scarce
(68 in the lowest bin). Lost sites overlap shared well and their answers are stable.

VERDICT ASSIGNED BY THIS SCRIPT
  robust        same direction on both matching axes, and crude agrees or is more extreme
  not identified  direction flips between the two matching axes
  dissolves     neither matched estimate supports the crude effect

Date: 2026-08-15
"""
import os, csv
import numpy as np
from scipy.stats import fisher_exact

Z = "/Users/alireza/Desktop/share_seq final/Zhang et al 2024 Nat Comm sc-RNA-seq Analysis"
T = os.path.join(Z, "E18p5_clean/results/chipseq_2026/tables")
EDG = [-2, 1.5, 2.0, 2.5, 3.0, 3.5, 8]
MINN = 25


def rd(name):
    return csv.DictReader(l for l in open(os.path.join(T, name))
                          if l.strip() and not l.startswith("#"))


S = {}
for r in rd("F4G_per_site_signal.csv"):
    try:
        w, m = float(r["wt_centre"]), float(r["mut_centre"])
    except (TypeError, ValueError):
        continue
    if not (np.isnan(w) or np.isnan(m)):
        S[(r["chr"], r["start"], r["end"])] = (w, m, r["cls"])
pos = {(r["chr"], r["start"], r["end"]): r for r in rd("F4_peak_classes_with_position.csv")}
el = {(r["chr"], r["start"], r["end"]): r["upstream"].upper()[6:11] == "TGCAG"
      for r in rd("F4C_umotif_panel_assignments.csv")}

PROPS = [
    ("upstream element", [k for k in S if k in el], lambda k: el[k]),
    ("promoter",  [k for k in S if k in pos], lambda k: pos[k]["position"] == "Promoter"),
    ("H3K4me3",   [k for k in S if k in pos], lambda k: pos[k]["H3K4me3"] == "TRUE"),
    ("H3K27ac",   [k for k in S if k in pos], lambda k: pos[k]["H3K27ac"] == "TRUE"),
    ("CTCF anchor", [k for k in S if k in pos], lambda k: pos[k]["ctcf_anchor"] == "TRUE"),
    ("distal intergenic", [k for k in S if k in pos],
     lambda k: pos[k]["position"] == "Distal intergenic"),
]


def mh(keys, fn, cls, axis):
    D = [(S[k][axis], S[k][2], fn(k)) for k in keys]
    g = [d for d in D if d[1] == cls]
    s = [d for d in D if d[1] == "shared"]
    num = den = 0.0
    det = []
    for i in range(len(EDG) - 1):
        lo, hi = EDG[i], EDG[i + 1]
        cs = [d for d in g if lo < d[0] <= hi]
        ss = [d for d in s if lo < d[0] <= hi]
        if len(cs) < MINN or len(ss) < MINN:
            continue
        a = sum(d[2] for d in cs); b = sum(d[2] for d in ss)
        c_, dd = len(cs) - a, len(ss) - b
        N = len(cs) + len(ss)
        num += a * dd / N; den += c_ * b / N
        o, _ = fisher_exact([[a, c_], [b, dd]])
        det.append(o)
    return (num / den if den else float("nan")), det


out = []
for name, keys, fn in PROPS:
    for cls in ("lost", "gained"):
        g = [k for k in keys if S[k][2] == cls]
        s = [k for k in keys if S[k][2] == "shared"]
        a = sum(fn(k) for k in g); b = sum(fn(k) for k in s)
        crude, pv = fisher_exact([[a, len(g) - a], [b, len(s) - b]])
        mw, dw = mh(keys, fn, cls, 0)
        mm, dm = mh(keys, fn, cls, 1)
        # direction agreement between the two matching axes decides the verdict
        sw, sm = np.sign(np.log(mw)), np.sign(np.log(mm))
        if sw != sm:
            verdict = "NOT IDENTIFIED"
        elif min(abs(np.log(mw)), abs(np.log(mm))) < np.log(1.25):
            verdict = "dissolves"
        else:
            verdict = "robust"
        out.append(dict(
            property=name, comparison=f"{cls} vs shared",
            pct_class=round(100 * a / len(g), 2), pct_shared=round(100 * b / len(s), 2),
            crude_OR=round(crude, 3), crude_p=pv,
            MH_OR_matched_wt=round(mw, 3), n_strata_wt=len(dw),
            MH_OR_matched_mut=round(mm, 3), n_strata_mut=len(dm),
            verdict=verdict,
            strata_OR_wt=" ".join(f"{d:.2f}" for d in dw),
            strata_OR_mut=" ".join(f"{d:.2f}" for d in dm)))

f = os.path.join(T, "F4_property_audit_stratified.csv")
with open(f, "w", newline="") as h:
    h.write("# Are the Figure 4 class differences genotype effects or CTCF signal strength?\n")
    h.write("# Mantel-Haenszel OR across signal strata, every class against SHARED.\n")
    h.write("# Run TWICE: matched on wild-type signal and matched on mutant signal.\n")
    h.write("# 'gained' is defined by mut minus wt, so matching on either axis conditions on a\n")
    h.write("#   collider. A real effect points the same way on both; an artefact can invert.\n")
    h.write("# NOT IDENTIFIED = the two axes disagree in sign. Do not claim a direction.\n")
    w = csv.DictWriter(h, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)

print(f"[140] wrote {os.path.basename(f)}\n")
print(f"{'property':<19}{'comparison':<18}{'crude':>8}{'MH wt':>8}{'MH mut':>8}   verdict")
print("-" * 82)
for v in ("robust", "NOT IDENTIFIED", "dissolves"):
    for r in out:
        if r["verdict"] == v:
            print(f"{r['property']:<19}{r['comparison']:<18}{r['crude_OR']:>8.3f}"
                  f"{r['MH_OR_matched_wt']:>8.3f}{r['MH_OR_matched_mut']:>8.3f}   {r['verdict']}")
