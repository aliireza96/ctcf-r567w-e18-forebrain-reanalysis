#!/usr/bin/env python3
"""
16_hic_insulation_tad_boundaries.py
Compute genome-wide insulation scores and call TAD boundaries from the public
Zhang et al. 2024 brain Hi-C files using hicstraw (KR normalisation, 25 kb
resolution). Outputs:
  - insulation_wt_25kb.bedgraph    per-bin insulation score, WT
  - insulation_homo_25kb.bedgraph  per-bin insulation score, HOMO
  - insulation_diff_25kb.bedgraph  HOMO - WT difference
  - tad_boundaries_wt_25kb.bed     boundary positions in WT (local minima)
  - tad_boundaries_homo_25kb.bed   boundary positions in HOMO
  - tad_boundaries_wt_only.bed     boundaries present in WT, lost in HOMO
  - tad_boundary_insulation_summary.csv  per-boundary stats for ChIP overlap

Usage:
  python3 scripts/16_hic_insulation_tad_boundaries.py
"""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

try:
    import hicstraw
except ImportError:
    sys.exit("hicstraw not installed. Run: pip install hicstraw")

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent
ROOT = SCRIPT_DIR.parent
HIC_DIR = ROOT.parent / "HIC_GSE214687_RAW"
OUT_DIR = ROOT / "results" / "discovery_hic_14"
TABLES_DIR = OUT_DIR / "tables"
TABLES_DIR.mkdir(parents=True, exist_ok=True)

WT_HIC   = HIC_DIR / "GSM6614279_brain_wt.allValidPairs.hic"
HOMO_HIC = HIC_DIR / "GSM6614277_brain_homo.allValidPairs.hic"

RESOLUTION   = 25_000   # 25 kb
WINDOW_BINS  = 20       # insulation window = 20 bins = 500 kb
MIN_BINS     = 40       # skip chromosomes shorter than this

# mm10 chromosomes to analyse (skip chrM, unplaced scaffolds)
CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX"]

# ── Helpers ──────────────────────────────────────────────────────────────────
def load_matrix(hic_path: Path, chrom: str, resolution: int) -> np.ndarray | None:
    """Return KR-normalised contact matrix for one chromosome at given resolution."""
    try:
        hic = hicstraw.HiCFile(str(hic_path))
        mzd = hic.getMatrixZoomData(chrom, chrom, "observed", "KR", "BP", resolution)
        records = mzd.getRecords(0, int(2e9), 0, int(2e9))
    except Exception as exc:
        print(f"  [WARN] {chrom} failed: {exc}", file=sys.stderr)
        return None

    if not records:
        return None

    # Determine matrix size from chromosome length
    chrom_sizes = {c.name: c.length for c in hicstraw.HiCFile(str(hic_path)).getChromosomes()}
    if chrom not in chrom_sizes:
        return None
    n = math.ceil(chrom_sizes[chrom] / resolution) + 1
    mat = np.zeros((n, n), dtype=np.float32)
    for r in records:
        i = r.binX // resolution
        j = r.binY // resolution
        if i < n and j < n:
            mat[i, j] = r.counts
            mat[j, i] = r.counts
    # Replace NaN/inf from KR normalisation failures
    mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
    return mat


def insulation_score(mat: np.ndarray, window: int) -> np.ndarray:
    """
    Compute insulation score per bin as log2(mean contact in off-diagonal
    diamond window / genome mean). Bins near the edges get NaN.
    """
    n = mat.shape[0]
    scores = np.full(n, np.nan, dtype=np.float64)
    total_signal = mat[mat > 0].mean() if (mat > 0).any() else 1.0

    for i in range(window, n - window):
        block = mat[i - window:i, i:i + window]
        mean_val = block.mean()
        if mean_val > 0 and total_signal > 0:
            scores[i] = np.log2(mean_val / total_signal)
        else:
            scores[i] = np.nan
    return scores


def call_boundaries(scores: np.ndarray, min_delta: float = 0.1) -> list[int]:
    """
    Return bin indices of local minima in insulation score where the score
    drop (from flanking bins) exceeds min_delta. These are TAD boundaries.
    """
    n = len(scores)
    boundaries = []
    for i in range(2, n - 2):
        if np.isnan(scores[i]):
            continue
        left  = np.nanmean(scores[max(0, i-3):i])
        right = np.nanmean(scores[i+1:min(n, i+4)])
        if scores[i] < left and scores[i] < right:
            delta = min(left - scores[i], right - scores[i])
            if delta >= min_delta:
                boundaries.append(i)
    return boundaries


# ── Main loop ────────────────────────────────────────────────────────────────
print(f"Resolution: {RESOLUTION // 1000} kb | Window: {WINDOW_BINS * RESOLUTION // 1000} kb")
print(f"WT  .hic : {WT_HIC}")
print(f"HOMO .hic: {HOMO_HIC}")

wt_bdry_rows   = []   # for bedgraph and BED output
homo_bdry_rows = []
diff_rows      = []
wt_only_rows   = []
summary_rows   = []

wt_bed_rows   = []
homo_bed_rows = []

for chrom in CHROMS:
    print(f"\n  Processing {chrom}...")
    mat_wt   = load_matrix(WT_HIC,   chrom, RESOLUTION)
    mat_homo = load_matrix(HOMO_HIC, chrom, RESOLUTION)

    if mat_wt is None or mat_homo is None:
        print(f"    Skipping {chrom} (matrix not available)")
        continue

    n_wt   = mat_wt.shape[0]
    n_homo = mat_homo.shape[0]
    n_bins = min(n_wt, n_homo)
    if n_bins < MIN_BINS:
        print(f"    Skipping {chrom} (only {n_bins} bins)")
        continue

    mat_wt   = mat_wt[:n_bins, :n_bins]
    mat_homo = mat_homo[:n_bins, :n_bins]

    ins_wt   = insulation_score(mat_wt,   WINDOW_BINS)
    ins_homo = insulation_score(mat_homo, WINDOW_BINS)
    ins_diff = ins_homo - ins_wt

    # Write insulation bedgraph rows
    for i in range(n_bins):
        start = i * RESOLUTION
        end   = start + RESOLUTION
        if not np.isnan(ins_wt[i]):
            wt_bdry_rows.append((chrom, start, end, ins_wt[i]))
        if not np.isnan(ins_homo[i]):
            homo_bdry_rows.append((chrom, start, end, ins_homo[i]))
        if not (np.isnan(ins_wt[i]) or np.isnan(ins_homo[i])):
            diff_rows.append((chrom, start, end, ins_diff[i]))

    # Call boundaries
    bdry_wt   = call_boundaries(ins_wt)
    bdry_homo = call_boundaries(ins_homo)
    bdry_homo_set = {b for b in bdry_homo}  # for overlap check

    for b in bdry_wt:
        start = b * RESOLUTION
        end   = start + RESOLUTION
        score = ins_wt[b]
        wt_bed_rows.append((chrom, start, end, f"{chrom}:{start}", round(score, 4)))

        # Check if lost in HOMO: no boundary within ±3 bins
        wt_score = ins_wt[b]
        homo_score = ins_homo[b] if b < len(ins_homo) else np.nan
        delta = homo_score - wt_score if not np.isnan(homo_score) else np.nan

        is_lost = not any(abs(b - bh) <= 3 for bh in bdry_homo_set)
        if is_lost:
            wt_only_rows.append((chrom, start, end, f"{chrom}:{start}",
                                  round(wt_score, 4) if not np.isnan(wt_score) else "NA",
                                  round(delta, 4) if not np.isnan(delta) else "NA"))

        summary_rows.append({
            "chrom": chrom, "start": start, "end": end,
            "wt_insulation": round(float(wt_score), 4) if not np.isnan(wt_score) else None,
            "homo_insulation": round(float(homo_score), 4) if not np.isnan(homo_score) else None,
            "delta_homo_minus_wt": round(float(delta), 4) if not np.isnan(delta) else None,
            "boundary_lost_in_homo": is_lost
        })

    for b in bdry_homo:
        start = b * RESOLUTION
        end   = start + RESOLUTION
        score = ins_homo[b]
        homo_bed_rows.append((chrom, start, end, f"{chrom}:{start}", round(score, 4)))

    print(f"    WT boundaries: {len(bdry_wt)}  |  HOMO boundaries: {len(bdry_homo)}  |  "
          f"WT-only (lost): {sum(1 for r in wt_only_rows if r[0] == chrom)}")


# ── Write outputs ─────────────────────────────────────────────────────────────
def write_bedgraph(rows, path):
    with open(path, "w") as f:
        for chrom, start, end, score in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{score:.6f}\n")
    print(f"  Written: {path}")

def write_bed(rows, path):
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")
    print(f"  Written: {path}")

print("\nWriting output files...")
write_bedgraph(wt_bdry_rows,   TABLES_DIR / "insulation_wt_25kb.bedgraph")
write_bedgraph(homo_bdry_rows, TABLES_DIR / "insulation_homo_25kb.bedgraph")
write_bedgraph(diff_rows,      TABLES_DIR / "insulation_diff_homo_minus_wt_25kb.bedgraph")
write_bed(wt_bed_rows,         TABLES_DIR / "tad_boundaries_wt_25kb.bed")
write_bed(homo_bed_rows,       TABLES_DIR / "tad_boundaries_homo_25kb.bed")
write_bed(wt_only_rows,        TABLES_DIR / "tad_boundaries_wt_only_lost_in_homo.bed")

# Summary CSV (for R ChIP overlap analysis)
if summary_rows:
    keys = ["chrom","start","end","wt_insulation","homo_insulation",
            "delta_homo_minus_wt","boundary_lost_in_homo"]
    with open(TABLES_DIR / "tad_boundary_insulation_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"  Written: {TABLES_DIR / 'tad_boundary_insulation_summary.csv'}")

print(f"\nDone. Total WT boundaries: {len(wt_bed_rows)}")
print(f"      WT-only (lost in HOMO): {len(wt_only_rows)}")
print(f"      HOMO boundaries: {len(homo_bed_rows)}")
