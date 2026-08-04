#!/usr/bin/env python3
"""
Build the quantitative Figure 4A-C tables on a true CTCF peak union.

The input BED must contain the disjoint union components written by
02_chipseq_peak_classification.R:

  chrom, start, end, union_id, score, strand, peak_class,
  n_wt_unique_intervals, n_mut_unique_intervals,
  n_wt_peak_records, n_mut_peak_records

Exact-coordinate MACS2 multi-summit records are collapsed before the union is
built, and all overlapping WT and mutant intervals are reduced to one genomic
component. Consequently, every base belongs to at most one row and Wilson
intervals do not treat duplicated MACS2 records as independent sites.

The sequence analysis is a position-anchored proxy motivated by Zhang et al.
(2024), not a reproduction of their FIMO-plus-Hamming-clustering U1 class:

  1. Read spike-in-normalised WT and mutant bigWig signal over each complete
     union component and calculate log2((mutant + 0.01)/(WT + 0.01)).
  2. Locate the best HOCOMOCO CTCF core-motif match within +/-250 bp of the
     component centre. The custom cutoff is a normal approximation to a
     nominal 1e-4 uniform-background tail probability; it is not a calibrated
     FIMO motif P value.
  3. Extract the strand-oriented upstream 20-mer and label it proxy-positive
     when it contains TGCAG.
  4. Rank callable components by occupancy change and split them into deciles.

Outputs in --outdir:
  Table_S7_CTCF_peak_signal.csv
  Table_S8_umotif_by_loss_decile.csv
  Table_S9_lost_site_location.csv          (when --annotations is supplied)
  figure4_normalised_class_summary.csv
  figure4_peak_umotif_status.csv

Requires numpy, scipy and pyBigWig; pyfaidx is additionally required when a
genome FASTA is supplied. For an already exported 500-bp sequence table, use
--sequence-tsv instead of --genome.
"""

import argparse
import csv
import math
import os

import numpy as np
from scipy.stats import norm


FLANK = 250
UPSTREAM = 20
P_THRESHOLD = 1e-4
N_DECILES = 10
U_SUBSTRING = "TGCAG"
VALID_CLASSES = {"LOST", "MAINTAINED", "GAINED"}
COMPLEMENT = str.maketrans("ACGTacgtN", "TGCAtgcaN")
BASE_INDEX = np.full(256, 4, dtype=np.int8)
for _base, _index in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    BASE_INDEX[ord(_base)] = _index
    BASE_INDEX[ord(_base.lower())] = _index


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--peaks", required=True, help="classified non-overlapping union-component BED")
    p.add_argument("--wt-bw", required=True, help="wild-type spike-in-normalised bigWig")
    p.add_argument("--mut-bw", required=True, help="mutant spike-in-normalised bigWig")
    seq = p.add_mutually_exclusive_group(required=True)
    seq.add_argument("--genome", help="mm10 genome FASTA (requires pyfaidx)")
    seq.add_argument(
        "--sequence-tsv",
        help="pre-exported 500-bp component-centred sequences with columns union_id/id and sequence",
    )
    p.add_argument("--pfm", required=True, help="HOCOMOCO CTCF position count matrix")
    p.add_argument(
        "--annotations",
        help="optional ChIPseeker CSV from 03_chipseq_annotation.R; enables Table S9",
    )
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_pfm(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            rows.append([float(x) for x in line.split()])
    counts = np.asarray(rows, dtype=float)
    if counts.ndim != 2 or counts.shape[1] != 4:
        raise ValueError(f"Expected a positions-by-4 A/C/G/T matrix in {path}")
    freqs = (counts + 0.25) / (counts.sum(axis=1, keepdims=True) + 1.0)
    return np.log2(freqs / 0.25)


def score_threshold(pwm, p_value):
    mu = pwm.mean(axis=1).sum()
    sd = np.sqrt(pwm.var(axis=1).sum())
    return float(mu + norm.isf(p_value) * sd)


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def scan(seq, pwm, cutoff):
    """Return the best (start in scanned-strand window, strand, score), or None."""
    width = pwm.shape[0]

    def strand_best(scanned):
        encoded = BASE_INDEX[np.frombuffer(scanned.encode("ascii"), dtype=np.uint8)]
        windows = np.lib.stride_tricks.sliding_window_view(encoded, width)
        valid = np.all(windows < 4, axis=1)
        safe = np.where(windows < 4, windows, 0)
        scores = pwm[np.arange(width)[None, :], safe].sum(axis=1)
        scores[~valid] = -np.inf
        index = int(np.argmax(scores))
        return index, float(scores[index])

    plus_i, plus_score = strand_best(seq)
    minus_i, minus_score = strand_best(revcomp(seq))
    # Match the original deterministic tie behaviour: forward strand wins ties.
    if plus_score >= minus_score:
        best = (plus_i, "+", plus_score)
    else:
        best = (minus_i, "-", minus_score)
    return best if best[2] >= cutoff else None


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return 100 * (centre - half), 100 * (centre + half)


def load_peaks(path):
    peaks = []
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                raise ValueError(f"{path}:{line_no}: expected at least seven BED columns")
            peak_class = fields[6] if fields[6] in VALID_CLASSES else fields[-1]
            if peak_class not in VALID_CLASSES:
                raise ValueError(f"{path}:{line_no}: unrecognised peak class {peak_class!r}")
            counts = [None, None, None, None]
            if len(fields) >= 11 and fields[6] in VALID_CLASSES:
                counts = [int(x) for x in fields[7:11]]
            peaks.append(
                {
                    "chrom": fields[0],
                    "start": int(fields[1]),
                    "end": int(fields[2]),
                    "name": fields[3],
                    "peak_class": peak_class,
                    "n_wt_unique_intervals": counts[0],
                    "n_mut_unique_intervals": counts[1],
                    "n_wt_peak_records": counts[2],
                    "n_mut_peak_records": counts[3],
                }
            )
    if not peaks:
        raise ValueError(f"No peaks read from {path}")
    for left, right in zip(peaks, peaks[1:]):
        if left["chrom"] == right["chrom"] and right["start"] < left["end"]:
            raise ValueError(
                "Input intervals overlap; run 02_chipseq_peak_classification.R to build the reduced union"
            )
    if len({peak["name"] for peak in peaks}) != len(peaks):
        raise ValueError("Union component identifiers are not unique")
    return peaks


def load_sequence_table(path, peaks):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None or "sequence" not in reader.fieldnames:
            raise ValueError(f"{path} must contain a sequence column")
        key_col = "union_id" if "union_id" in reader.fieldnames else "id"
        rows = list(reader)
    by_id = {str(row[key_col]): row["sequence"].upper() for row in rows}
    sequences = []
    for index, peak in enumerate(peaks, 1):
        sequence = by_id.get(peak["name"])
        if sequence is None:
            sequence = by_id.get(str(index))
        if sequence is None and len(rows) == len(peaks):
            sequence = rows[index - 1]["sequence"].upper()
        sequences.append(sequence)
    return sequences


def fetch_sequences_from_fasta(path, peaks):
    try:
        from pyfaidx import Fasta
    except ImportError as exc:
        raise SystemExit("--genome requires pyfaidx; install it or use --sequence-tsv") from exc
    genome = Fasta(path, as_raw=True, sequence_always_upper=True)
    sequences = []
    for peak in peaks:
        centre = (peak["start"] + peak["end"]) // 2
        lo, hi = max(0, centre - FLANK), centre + FLANK
        try:
            sequence = str(genome[peak["chrom"]][lo:hi]).upper()
        except (KeyError, ValueError):
            sequence = None
        sequences.append(sequence)
    genome.close()
    return sequences


def write_commented_csv(path, comments, fieldnames, rows):
    with open(path, "w", newline="") as fh:
        for comment in comments:
            fh.write(f"# {comment}\n")
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_annotations(path):
    with open(path, newline="") as fh:
        reader = csv.DictReader(row for row in fh if not row.startswith("#"))
        rows = list(reader)
    if not rows or "union_id" not in rows[0]:
        raise ValueError("Annotation table must contain union_id")
    out = {}
    for row in rows:
        feature = row.get("figure4_feature", "")
        if feature not in {"Promoter", "Intron", "Distal", "Other"}:
            annotation = row.get("annotation", "").lower()
            if "promoter" in annotation:
                feature = "Promoter"
            elif "intron" in annotation:
                feature = "Intron"
            elif "intergenic" in annotation:
                feature = "Distal"
            else:
                feature = "Other"
        out[row["union_id"]] = feature
    return out


def main():
    args = parse_args()
    import pyBigWig

    os.makedirs(args.outdir, exist_ok=True)
    peaks = load_peaks(args.peaks)
    if len(peaks) != 45021:
        raise ValueError(f"Expected 45,021 union components, found {len(peaks):,}")

    class_counts = {label: sum(p["peak_class"] == label for p in peaks) for label in VALID_CLASSES}
    expected_counts = {"MAINTAINED": 28964, "LOST": 6666, "GAINED": 9391}
    if class_counts != expected_counts:
        raise ValueError(f"Unexpected union class counts: {class_counts}")

    bw_wt = pyBigWig.open(args.wt_bw)
    bw_mut = pyBigWig.open(args.mut_bw)
    for peak in peaks:
        wt = bw_wt.stats(peak["chrom"], peak["start"], peak["end"], type="mean")[0]
        mut = bw_mut.stats(peak["chrom"], peak["start"], peak["end"], type="mean")[0]
        wt = 0.0 if wt is None or not math.isfinite(wt) else float(wt)
        mut = 0.0 if mut is None or not math.isfinite(mut) else float(mut)
        peak["wt_norm_signal"] = wt
        peak["mut_norm_signal"] = mut
        peak["log2FC_norm"] = float(np.log2((mut + 0.01) / (wt + 0.01)))
    bw_wt.close()
    bw_mut.close()

    table_s7_fields = [
        "chrom",
        "start",
        "end",
        "union_id",
        "peak_class",
        "n_wt_unique_intervals",
        "n_mut_unique_intervals",
        "n_wt_peak_records",
        "n_mut_peak_records",
        "wt_norm_signal",
        "mut_norm_signal",
        "log2FC_norm",
    ]
    table_s7_rows = []
    for peak in peaks:
        row = dict(peak)
        row["union_id"] = peak["name"]
        for field in ("wt_norm_signal", "mut_norm_signal", "log2FC_norm"):
            row[field] = f"{row[field]:.8f}"
        table_s7_rows.append(row)
    write_commented_csv(
        os.path.join(args.outdir, "Table_S7_CTCF_peak_signal.csv"),
        [
            "One row per reduced non-overlapping CTCF union component (mm10, BED-style 0-based start).",
            "Exact-coordinate MACS2 multi-summit records were collapsed before WT and mutant intervals were reduced.",
            "Signal is the mean of the deposited spike-in-normalised bigWig over the complete union component.",
            "log2FC_norm = log2((mutant + 0.01)/(wild type + 0.01)).",
        ],
        table_s7_fields,
        table_s7_rows,
    )

    pwm = read_pfm(args.pfm)
    cutoff = score_threshold(pwm, P_THRESHOLD)
    sequences = (
        load_sequence_table(args.sequence_tsv, peaks)
        if args.sequence_tsv
        else fetch_sequences_from_fasta(args.genome, peaks)
    )
    for peak, window in zip(peaks, sequences):
        peak.update(
            core_score=None,
            core_strand=None,
            motif_start_in_window=None,
            upstream20=None,
            u_proxy=False,
            callable=False,
            decile_rank=None,
        )
        if window is None or len(window) != 2 * FLANK:
            continue
        hit = scan(window, pwm, cutoff)
        if hit is None:
            continue
        start_in_scanned, strand, score = hit
        scanned = window if strand == "+" else revcomp(window)
        upstream = scanned[max(0, start_in_scanned - UPSTREAM) : start_in_scanned]
        peak["core_score"] = score
        peak["core_strand"] = strand
        peak["motif_start_in_window"] = start_in_scanned
        if len(upstream) != UPSTREAM:
            continue
        peak["upstream20"] = upstream
        peak["u_proxy"] = U_SUBSTRING in upstream
        peak["callable"] = True

    callable_peaks = [peak for peak in peaks if peak["callable"]]
    callable_peaks.sort(key=lambda peak: (peak["log2FC_norm"], peak["name"]))
    edges = np.linspace(0, len(callable_peaks), N_DECILES + 1).astype(int)
    decile_rows = []
    decile_chunks = []
    for index in range(N_DECILES):
        chunk = callable_peaks[edges[index] : edges[index + 1]]
        decile_chunks.append(chunk)
        for peak in chunk:
            peak["decile_rank"] = index + 1
        positives = sum(peak["u_proxy"] for peak in chunk)
        ci_lo, ci_hi = wilson(positives, len(chunk))
        decile_rows.append(
            {
                "decile_rank": index + 1,
                "log2FC_min": f"{min(p['log2FC_norm'] for p in chunk):.6f}",
                "log2FC_max": f"{max(p['log2FC_norm'] for p in chunk):.6f}",
                "median_log2FC_norm": f"{np.median([p['log2FC_norm'] for p in chunk]):.6f}",
                "n_components": len(chunk),
                "n_proxy_positive": positives,
                "pct_proxy_positive": f"{100 * positives / len(chunk):.4f}",
                "ci_lo": f"{ci_lo:.4f}",
                "ci_hi": f"{ci_hi:.4f}",
            }
        )

    write_commented_csv(
        os.path.join(args.outdir, "Table_S8_umotif_by_loss_decile.csv"),
        [
            "TGCAG-containing, core-motif-anchored upstream-20-mer proxy across occupancy-change deciles.",
            "Decile 1 = most binding lost; decile 10 = most binding gained. CI = Wilson 95%.",
            "This custom proxy is related to but not identical with Zhang et al.'s FIMO-plus-clustering U1 class.",
        ],
        list(decile_rows[0]),
        decile_rows,
    )

    class_summary_rows = []
    for label in ("MAINTAINED", "LOST", "GAINED"):
        subset = [peak for peak in peaks if peak["peak_class"] == label]
        med_wt = float(np.median([peak["wt_norm_signal"] for peak in subset]))
        med_mut = float(np.median([peak["mut_norm_signal"] for peak in subset]))
        med_fc = float(np.median([peak["log2FC_norm"] for peak in subset]))
        class_summary_rows.append(
            {
                "peak_class": label,
                "n_components": len(subset),
                "median_wt_norm": f"{med_wt:.6f}",
                "median_mut_norm": f"{med_mut:.6f}",
                "median_log2FC": f"{med_fc:.6f}",
                "pct_components_lower_in_mutant": f"{100 * np.mean([p['log2FC_norm'] < 0 for p in subset]):.4f}",
            }
        )
    write_commented_csv(
        os.path.join(args.outdir, "figure4_normalised_class_summary.csv"),
        ["Per-class summary of the 45,021 reduced non-overlapping union components."],
        list(class_summary_rows[0]),
        class_summary_rows,
    )

    status_fields = [
        "union_id",
        "chrom",
        "start",
        "end",
        "peak_class",
        "wt_norm_signal",
        "mut_norm_signal",
        "log2FC_norm",
        "core_score",
        "core_strand",
        "motif_start_in_window",
        "upstream20",
        "u_proxy",
        "callable",
        "decile_rank",
    ]
    status_rows = []
    for peak in peaks:
        row = dict(peak)
        row["union_id"] = peak["name"]
        for field in ("wt_norm_signal", "mut_norm_signal", "log2FC_norm"):
            row[field] = f"{row[field]:.8f}"
        if peak["core_score"] is not None:
            row["core_score"] = f"{peak['core_score']:.8f}"
        status_rows.append(row)
    write_commented_csv(
        os.path.join(args.outdir, "figure4_peak_umotif_status.csv"),
        [
            "Per-component signal and position-anchored sequence-proxy status for all 45,021 union components.",
            "callable is true only when a core match and complete upstream 20-mer were obtained.",
        ],
        status_fields,
        status_rows,
    )

    if args.annotations:
        annotations = load_annotations(args.annotations)
        missing = [peak["name"] for peak in peaks if peak["name"] not in annotations]
        if missing:
            raise ValueError(f"Annotation table lacks {len(missing):,} union components; first is {missing[0]}")
        location_rows = []
        feature_order = ("Promoter", "Intron", "Distal", "Other")
        for rank, chunk in enumerate(decile_chunks, 1):
            counts = {feature: sum(annotations[p["name"]] == feature for p in chunk) for feature in feature_order}
            row = {
                "decile_rank": rank,
                "median_log2FC_norm": f"{np.median([p['log2FC_norm'] for p in chunk]):.6f}",
                "n_components": len(chunk),
            }
            for feature in feature_order:
                key = feature.lower()
                row[f"n_{key}"] = counts[feature]
                row[f"pct_{key}"] = f"{100 * counts[feature] / len(chunk):.4f}"
            location_rows.append(row)
        write_commented_csv(
            os.path.join(args.outdir, "Table_S9_lost_site_location.csv"),
            [
                "ChIPseeker genomic-position composition across all callable-component occupancy deciles.",
                "Promoter = +/-3 kb from a known-gene TSS; percentages sum to 100 within each decile.",
            ],
            list(location_rows[0]),
            location_rows,
        )

    print(
        f"{len(peaks):,} non-overlapping components; {len(callable_peaks):,} callable; "
        f"{sum(p['u_proxy'] for p in callable_peaks):,} TGCAG-upstream-proxy positive"
    )
    print(
        f"Proxy frequency: {decile_rows[0]['pct_proxy_positive']}% (most lost) to "
        f"{decile_rows[-1]['pct_proxy_positive']}% (most gained)"
    )


if __name__ == "__main__":
    main()
