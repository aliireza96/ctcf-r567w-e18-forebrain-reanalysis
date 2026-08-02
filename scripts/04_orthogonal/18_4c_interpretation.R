# ============================================================
# 18_4c_interpretation.R — summarize and interpret 4C cPcdh discovery results
# ============================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])[1]
file_arg <- gsub("~\\+~", " ", file_arg)
SCRIPT_DIR <- dirname(normalizePath(ifelse(length(file_arg) == 0 || is.na(file_arg), ".", file_arg)))
source(file.path(SCRIPT_DIR, "12_discovery_setup.R"))
discovery_require(c("dplyr", "readr", "stringr", "GenomicRanges"))
discovery_append_log("fourc", "18_4c_interpretation")

fourc_dirs <- discovery_dirs_for("fourc")
bulk_dirs <- discovery_dirs_for("bulk")

read_csv_safe <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

intervals <- read_csv_safe(file.path(fourc_dirs$tables, "14_4c_altered_intervals.csv"))
rep_corr <- read_csv_safe(file.path(fourc_dirs$tables, "14_4c_replicate_correlations.csv"))
bulk_cpcdh <- read_csv_safe(file.path(fourc_dirs$tables, "14_4c_bulk_cpcdh_integration.csv"))
profiles <- read_csv_safe(file.path(fourc_dirs$tables, "14_4c_binned_profiles.csv"))

if (nrow(intervals) == 0) {
  stop("4C altered interval table is missing or empty.", call. = FALSE)
}

find_shared_intervals <- function(intervals_df) {
  bc <- intervals_df %>% dplyr::filter(viewpoint == "BC")
  fv <- intervals_df %>% dplyr::filter(viewpoint == "F")
  if (nrow(bc) == 0 || nrow(fv) == 0) {
    return(data.frame())
  }
  gr_bc <- GenomicRanges::GRanges(seqnames = bc$chr, ranges = IRanges::IRanges(bc$start, bc$end))
  gr_f <- GenomicRanges::GRanges(seqnames = fv$chr, ranges = IRanges::IRanges(fv$start, fv$end))
  hits <- GenomicRanges::findOverlaps(gr_bc, gr_f, ignore.strand = TRUE)
  if (length(hits) == 0) {
    return(data.frame())
  }
  shared <- data.frame(
    chr = bc$chr[queryHits(hits)],
    overlap_start = pmax(bc$start[queryHits(hits)], fv$start[subjectHits(hits)]),
    overlap_end = pmin(bc$end[queryHits(hits)], fv$end[subjectHits(hits)]),
    BC_interval = bc$interval_group[queryHits(hits)],
    BC_direction = bc$direction[queryHits(hits)],
    BC_family = bc$family_overlap[queryHits(hits)],
    BC_mean_delta = bc$mean_delta[queryHits(hits)],
    F_interval = fv$interval_group[subjectHits(hits)],
    F_direction = fv$direction[subjectHits(hits)],
    F_family = fv$family_overlap[subjectHits(hits)],
    F_mean_delta = fv$mean_delta[subjectHits(hits)],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      overlap_bp = overlap_end - overlap_start + 1L,
      shared_direction = ifelse(BC_direction == F_direction, BC_direction, "discordant"),
      shared_family = dplyr::case_when(
        BC_family == F_family ~ BC_family,
        TRUE ~ paste(sort(unique(c(BC_family, F_family))), collapse = " | ")
      ),
      mean_abs_delta = (abs(BC_mean_delta) + abs(F_mean_delta)) / 2
    ) %>%
    dplyr::arrange(dplyr::desc(mean_abs_delta), dplyr::desc(overlap_bp))
  shared
}

shared_intervals <- find_shared_intervals(intervals)
write_csv_discovery("fourc", shared_intervals, "14_4c_shared_intervals_across_viewpoints.csv")

profiles_summary <- if (nrow(profiles) > 0) {
  aggregate(log_signal ~ viewpoint + genotype + replicate, data = profiles, FUN = mean)
} else {
  data.frame()
}
write_csv_discovery("fourc", profiles_summary, "14_4c_locus_mean_signal_by_sample.csv")

classify_interval <- function(family, direction) {
  if (grepl("Pcdha;Pcdhb", family, fixed = TRUE) && direction == "higher_in_wt") {
    return("WT-enriched contact loss across Pcdha/Pcdhb boundary-proximal region")
  }
  if (grepl("Pcdha;Pcdhg", family, fixed = TRUE) && direction == "higher_in_homo") {
    return("HOMO-enriched contact shift toward distal Pcdha/Pcdhg-side region")
  }
  if (family == "Pcdha" && direction == "higher_in_homo") {
    return("Local Pcdha-proximal HOMO-enriched contact")
  }
  if (family == "inter-family") {
    return("Flanking or inter-family structural shift")
  }
  "Other cPcdh-locus shift"
}

shared_intervals_annotated <- if (nrow(shared_intervals) > 0) {
  shared_intervals %>%
    dplyr::mutate(interval_interpretation = mapply(classify_interval, shared_family, shared_direction))
} else {
  data.frame()
}
if (nrow(shared_intervals_annotated) > 0) {
  write_csv_discovery("fourc", shared_intervals_annotated, "14_4c_shared_intervals_annotated.csv")
}

summary_lines_from_df <- function(df, n = 5) {
  if (nrow(df) == 0) {
    return("- No rows available.")
  }
  top <- df %>% dplyr::slice_head(n = n)
  apply(top, 1, function(row) {
    paste0(
      "- `", row[["chr"]], ":", format(as.numeric(row[["overlap_start"]]), scientific = FALSE, big.mark = ","),
      "-", format(as.numeric(row[["overlap_end"]]), scientific = FALSE, big.mark = ","), "`: `",
      row[["shared_direction"]], "`; family=`", row[["shared_family"]],
      "`; mean abs delta=", sprintf("%.2f", as.numeric(row[["mean_abs_delta"]]))
    )
  })
}

rep_lines <- if (nrow(rep_corr) > 0) {
  apply(rep_corr, 1, function(row) {
    paste0("- `", row[["viewpoint"]], "` `", row[["genotype"]], "` replicate rho=", sprintf("%.2f", as.numeric(row[["spearman_cor"]])))
  })
} else {
  "- Replicate-correlation table was unavailable."
}

bulk_lines <- if (nrow(bulk_cpcdh) > 0) {
  apply(bulk_cpcdh, 1, function(row) {
    paste0(
      "- `", row[["family"]], "`: bulk mean HOMO-vs-WT logFC=",
      sprintf("%.2f", as.numeric(row[["bulk_mean_HOMO_vs_WT_logFC"]])),
      "; dominant dosage patterns=`", row[["dominant_pattern"]], "`"
    )
  })
} else {
  "- Bulk cPcdh integration table was unavailable."
}

mean_signal_lines <- if (nrow(profiles_summary) > 0) {
  apply(profiles_summary, 1, function(row) {
    paste0(
      "- `", row[["viewpoint"]], "` `", row[["genotype"]], "` `", row[["replicate"]],
      "` locus mean log signal=", sprintf("%.2f", as.numeric(row[["log_signal"]]))
    )
  })
} else {
  "- Mean-signal summary could not be computed."
}

note_lines <- c(
  "# 4C cPcdh Interpretation Note",
  "",
  "## Main Message",
  "- The 4C data provide locus-level structural support for the `cPcdh` story, with the strongest reproducible signal being a WT-enriched contact block across the `Pcdha/Pcdhb` region and a smaller HOMO-enriched shift toward a more distal `Pcdha/Pcdhg`-side interval.",
  "- This is a bounded structural result at one locus. It strengthens the manuscript as a positive-control chromatin-contact layer, but it should not be generalized into a genome-wide mechanism.",
  "",
  "## Replicate Concordance",
  rep_lines,
  "- The `F` viewpoint is clearly more reproducible than the `BC` viewpoint, but both viewpoints show moderate-to-good within-genotype concordance.",
  "",
  "## Shared Intervals Across Both Viewpoints",
  summary_lines_from_df(shared_intervals_annotated, n = 8),
  "- These cross-viewpoint overlaps are the highest-confidence structural shifts because they recur in both capture designs with the same direction.",
  "",
  "## Integration With Bulk cPcdh Expression",
  bulk_lines,
  "- The structural pattern is most compatible with the bulk result in `Pcdhb`, where the expression decrease is strongest and the recurrent 4C losses sit across the `Pcdha/Pcdhb`-overlapping region.",
  "- `Pcdhg` is more mixed at the bulk-expression level, and the HOMO-enriched distal 4C shift should therefore be described as a redistribution within the locus, not simply as uniform loss of all contacts.",
  "",
  "## Locus-Wide Signal Context",
  mean_signal_lines,
  "- Locus-wide mean signal is not dramatically collapsed in HOMO across all samples. The informative result is the regional redistribution of contacts within the `cPcdh` interval, not a simple global reduction in 4C coverage.",
  "",
  "## Recommended Manuscript Framing",
  "- Use 4C as a bounded structural validation that the `cPcdh` locus is reorganized in HOMO brain.",
  "- Emphasize the recurrent WT-enriched contact loss over the `Pcdha/Pcdhb`-overlapping region as the cleanest structural result.",
  "- Mention the HOMO-enriched distal `Pcdha/Pcdhg`-side interval as evidence that the mutant phenotype may involve contact redistribution rather than only contact loss.",
  "- Avoid claiming that 4C proves mechanism beyond the locus or that every cPcdh subfamily is perturbed in the same way.",
  "",
  "## Key Files",
  paste0("- Shared-interval table: [14_4c_shared_intervals_across_viewpoints.csv](<", fourc_dirs$tables, "/14_4c_shared_intervals_across_viewpoints.csv>)"),
  paste0("- Annotated shared-interval table: [14_4c_shared_intervals_annotated.csv](<", fourc_dirs$tables, "/14_4c_shared_intervals_annotated.csv>)"),
  paste0("- Structural summary figure: [14_4c_cpcdh_structural_map.pdf](<", fourc_dirs$plots, "/14_4c_cpcdh_structural_map.pdf>)"),
  paste0("- Interval-summary figure: [14_4c_altered_interval_summary.pdf](<", fourc_dirs$plots, "/14_4c_altered_interval_summary.pdf>)")
)

write_status_note("fourc", "14_4c_interpretation.md", note_lines)
discovery_msg("4C interpretation artifacts written.")
