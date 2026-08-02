#!/usr/bin/env python3
"""
Upstream (U) motif frequency across deciles of CTCF occupancy change.

Reproduces the analysis behind Figure 4B and Supplementary Table S8, and the
per-peak signal table (Supplementary Table S7).

Method, following Zhang et al. (2024):

  1. Peaks are classified maintained / lost / gained by reciprocal overlap
     between the two genotypes (see 02_chipseq_peak_classification.R).
  2. Normalised ChIP signal is read from the deposited bigWig tracks over each
     peak interval, and the log2 ratio between genotypes gives the occupancy
     change. The deposited tracks are already normalised; no further scaling is
     applied.
  3. The CTCF core motif is located within +/-250 bp of each peak centre using
     the HOCOMOCO mouse PFM CTCF_MOUSE.H11MO.0.A at a match threshold of
     p < 1e-4, and the 20 bp immediately upstream of that core motif, on the
     motif's own strand, is extracted.
  4. Peaks are ranked by occupancy change and split into deciles. U-motif
     frequency is computed per decile with Wilson 95% intervals.

The position anchoring in step 3 is essential. Scanning a fixed consensus
sequence across whole peaks does not recover the effect, because the U motif is
defined by its position relative to the core motif rather than by a sequence
occurring anywhere in the peak.

Requires: pyBigWig, pyfaidx, numpy, scipy, and an mm10 genome FASTA.

Usage:
    python 11_chipseq_umotif_by_occupancy.py \
        --peaks    results/chipseq/tables/chipseq_02_all_peaks_classified.bed \
        --wt-bw    data/chipseq/Brain_CTCF_wt.mm10.bw \
        --mut-bw   data/chipseq/Brain_CTCF_homo.mm10.bw \
        --genome   data/genome/mm10.fa \
        --pfm      CTCF_MOUSE.H11MO.0.A.pcm \
        --outdir   results/chipseq/tables
"""
import argparse
import csv
import os

import numpy as np
from scipy.stats import norm

FLANK = 250          # search window either side of the peak centre
UPSTREAM = 20        # length of the upstream window, in bp
P_THRESHOLD = 1e-4   # motif match threshold
N_DECILES = 10
U_CONSENSUS = ("TGCAG", "CTGCAG")   # substrings defining a U-motif-positive 20-mer
COMPLEMENT = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peaks", required=True, help="classified peak BED (chrom start end name . . class)")
    p.add_argument("--wt-bw", required=True, help="wild-type normalised bigWig")
    p.add_argument("--mut-bw", required=True, help="mutant normalised bigWig")
    p.add_argument("--genome", required=True, help="mm10 genome FASTA")
    p.add_argument("--pfm", required=True, help="HOCOMOCO CTCF position count matrix")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_pfm(path):
    """Read a HOCOMOCO .pcm file and return a log-odds PWM against uniform background."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            rows.append([float(x) for x in line.split()])
    counts = np.asarray(rows)              # positions x 4, in ACGT order
    freqs = (counts + 0.25) / (counts.sum(axis=1, keepdims=True) + 1.0)
    return np.log2(freqs / 0.25)


def score_threshold(pwm, p_value):
    """Score cutoff corresponding to p_value under a uniform-background null."""
    mu = pwm.mean(axis=1).sum()
    sd = np.sqrt((pwm.var(axis=1)).sum())
    return mu + norm.isf(p_value) * sd


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def scan(seq, pwm, cutoff):
    """Best PWM hit on either strand. Returns (start, strand, score) or None."""
    w = pwm.shape[0]
    idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    best = None
    for strand, s in (("+", seq), ("-", revcomp(seq))):
        for i in range(len(s) - w + 1):
            sub = s[i:i + w].upper()
            if any(c not in idx for c in sub):
                continue
            sc = sum(pwm[j, idx[c]] for j, c in enumerate(sub))
            if sc >= cutoff and (best is None or sc > best[2]):
                best = (i, strand, sc)
    return best


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (100 * (centre - half), 100 * (centre + half))


def load_peaks(path):
    peaks = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3 or f[0].startswith("track"):
                continue
            peaks.append(dict(chrom=f[0], start=int(f[1]), end=int(f[2]),
                              name=f[3] if len(f) > 3 else f"{f[0]}:{f[1]}",
                              signal=float(f[4]) if len(f) > 4 and f[4] not in (".", "") else float("nan"),
                              peak_class=f[-1] if len(f) > 6 else ""))
    return peaks


def main():
    args = parse_args()
    import pyBigWig
    from pyfaidx import Fasta

    os.makedirs(args.outdir, exist_ok=True)
    peaks = load_peaks(args.peaks)
    pwm = read_pfm(args.pfm)
    cutoff = score_threshold(pwm, P_THRESHOLD)

    bw_wt = pyBigWig.open(args.wt_bw)
    bw_mut = pyBigWig.open(args.mut_bw)
    genome = Fasta(args.genome, as_raw=True, sequence_always_upper=True)

    # --- per-peak normalised signal and occupancy change ---------------------
    for pk in peaks:
        wt = bw_wt.stats(pk["chrom"], pk["start"], pk["end"], type="mean")[0] or 0.0
        mut = bw_mut.stats(pk["chrom"], pk["start"], pk["end"], type="mean")[0] or 0.0
        pk["wt_norm"], pk["mut_norm"] = wt, mut
        pk["log2fc"] = np.log2((mut + 0.01) / (wt + 0.01))
    bw_wt.close()
    bw_mut.close()

    with open(os.path.join(args.outdir, "Table_S7_CTCF_peak_signal.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# Per-peak normalised CTCF signal. log2FC_norm = log2(mutant / wild type)"])
        w.writerow(["# over the peak interval, from the deposited normalised bigWig tracks."])
        w.writerow([])
        w.writerow(["chrom", "start", "end", "peak_name", "peak_class",
                    "macs2_signalValue", "wt_norm_signal", "mut_norm_signal", "log2FC_norm"])
        for pk in peaks:
            w.writerow([pk["chrom"], pk["start"], pk["end"], pk["name"], pk["peak_class"],
                        round(pk["signal"], 4) if pk["signal"] == pk["signal"] else "",
                        round(pk["wt_norm"], 4), round(pk["mut_norm"], 4), round(pk["log2fc"], 4)])

    # --- core motif and upstream 20-mer -------------------------------------
    for pk in peaks:
        centre = (pk["start"] + pk["end"]) // 2
        lo, hi = max(0, centre - FLANK), centre + FLANK
        try:
            window = str(genome[pk["chrom"]][lo:hi])
        except (KeyError, ValueError):
            pk["upstream"] = None
            continue
        hit = scan(window, pwm, cutoff)
        if hit is None:
            pk["upstream"] = None
            continue
        i, strand, _ = hit
        if strand == "+":
            up = window[max(0, i - UPSTREAM):i]
        else:
            rc = revcomp(window)
            up = rc[max(0, i - UPSTREAM):i]
        pk["upstream"] = up if len(up) == UPSTREAM else None

    scored = [pk for pk in peaks if pk["upstream"]]
    for pk in scored:
        pk["u_motif"] = any(c in pk["upstream"] for c in U_CONSENSUS)

    # --- deciles of occupancy change ----------------------------------------
    # rank 1 = most binding lost, rank 10 = most gained
    scored.sort(key=lambda x: x["log2fc"])
    edges = np.linspace(0, len(scored), N_DECILES + 1).astype(int)
    rows = []
    for d in range(N_DECILES):
        chunk = scored[edges[d]:edges[d + 1]]
        k = sum(1 for pk in chunk if pk["u_motif"])
        lo, hi = wilson(k, len(chunk))
        rows.append(dict(decile_rank=d + 1,
                         median_log2FC_norm=round(float(np.median([pk["log2fc"] for pk in chunk])), 4),
                         n_peaks=len(chunk), n_umotif_pos=k,
                         pct_umotif=round(100 * k / len(chunk), 2),
                         ci_lo=round(lo, 2), ci_hi=round(hi, 2)))

    with open(os.path.join(args.outdir, "Table_S8_umotif_by_loss_decile.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["# Panel B source. CTCF sites ranked by normalised occupancy change and split into deciles."])
        w.writerow(["# decile_rank 1 = most binding lost, 10 = binding gained. CI = Wilson 95%."])
        w.writerow([])
        w.writerow(list(rows[0].keys()))
        for r in rows:
            w.writerow(list(r.values()))

    print(f"{len(peaks)} peaks, {len(scored)} with a core-motif hit and usable upstream 20-mer")
    print(f"U-motif frequency: {rows[0]['pct_umotif']}% (most lost) to {rows[-1]['pct_umotif']}% (most gained)")


if __name__ == "__main__":
    main()
