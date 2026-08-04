# ============================================================
# 02_chipseq_peak_classification.R — build a true CTCF peak union
# ============================================================

source(file.path("scripts", "chipseq", "00_chipseq_setup.R"))
chipseq_require(c("GenomicRanges", "rtracklayer", "ggplot2", "dplyr", "readr", "scales"))

peak_obj_path <- file.path(CHIP_DIRS$objects, "chipseq_01_peak_objects.rds")
if (!file.exists(peak_obj_path)) {
  stop("Run 01_chipseq_peak_overview.R first: peak objects not found.", call. = FALSE)
}
peak_obj <- readRDS(peak_obj_path)
wt_peaks <- peak_obj$wt
mut_peaks <- peak_obj$mut

# MACS2 can emit several summits (a/b/c suffixes) for one peak interval. Those
# records are not independent genomic sites. Collapse exact-coordinate records
# before constructing the union, then reduce all WT and mutant intervals into
# disjoint connected components. A component is MAINTAINED when it contains at
# least one interval from each genotype, LOST when it contains WT only, and
# GAINED when it contains mutant only.
collapse_exact_coordinates <- function(gr, source_label) {
  key <- paste0(as.character(GenomicRanges::seqnames(gr)), ":", start(gr), "-", end(gr))
  keep <- !duplicated(key)
  out <- gr[keep]
  counts <- table(key)
  out$n_source_records <- as.integer(counts[key[keep]])
  out$source_genotype <- source_label
  out
}

aggregate_component_values <- function(components, source_gr, values, fun, default = NA_real_) {
  hits <- GenomicRanges::findOverlaps(components, source_gr, minoverlap = 1L, ignore.strand = TRUE)
  out <- rep(default, length(components))
  if (length(hits) == 0L) return(out)
  grouped <- split(values[S4Vectors::subjectHits(hits)], S4Vectors::queryHits(hits))
  safe_apply <- function(x) {
    if (all(is.na(x))) return(default)
    fun(x[!is.na(x)])
  }
  vals <- vapply(grouped, safe_apply, numeric(1))
  out[as.integer(names(vals))] <- vals
  out
}

count_component_overlaps <- function(components, source_gr) {
  as.integer(GenomicRanges::countOverlaps(components, source_gr, minoverlap = 1L, ignore.strand = TRUE))
}

wt_unique <- collapse_exact_coordinates(wt_peaks, "WT")
mut_unique <- collapse_exact_coordinates(mut_peaks, "MUT")

union_components <- GenomicRanges::reduce(
  c(wt_unique, mut_unique),
  min.gapwidth = 1L,
  ignore.strand = TRUE
)
strand(union_components) <- "*"

union_components$n_wt_unique_intervals <- count_component_overlaps(union_components, wt_unique)
union_components$n_mut_unique_intervals <- count_component_overlaps(union_components, mut_unique)
union_components$n_wt_peak_records <- count_component_overlaps(union_components, wt_peaks)
union_components$n_mut_peak_records <- count_component_overlaps(union_components, mut_peaks)
union_components$class <- dplyr::case_when(
  union_components$n_wt_unique_intervals > 0L & union_components$n_mut_unique_intervals > 0L ~ "MAINTAINED",
  union_components$n_wt_unique_intervals > 0L ~ "LOST",
  TRUE ~ "GAINED"
)
union_components$union_id <- sprintf("CTCF_union_%05d", seq_along(union_components))
union_components$name <- union_components$union_id
union_components$score <- 0

# Retain explicit source-track metadata for compatibility with downstream
# exploratory scripts. The Figure 4 quantitative analysis uses the deposited
# spike-in-normalised bigWigs, not these MACS2 fields.
union_components$wt_max_signalValue <- aggregate_component_values(
  union_components, wt_peaks, wt_peaks$signalValue, max
)
union_components$mut_max_signalValue <- aggregate_component_values(
  union_components, mut_peaks, mut_peaks$signalValue, max
)
union_components$wt_max_qValue <- aggregate_component_values(
  union_components, wt_peaks, wt_peaks$qValue, max
)
union_components$mut_max_qValue <- aggregate_component_values(
  union_components, mut_peaks, mut_peaks$qValue, max
)
union_components$signalValue <- ifelse(
  union_components$n_wt_unique_intervals > 0L,
  union_components$wt_max_signalValue,
  union_components$mut_max_signalValue
)
union_components$qValue <- ifelse(
  union_components$n_wt_unique_intervals > 0L,
  union_components$wt_max_qValue,
  union_components$mut_max_qValue
)
union_components$pValue <- NA_real_
union_components$peak <- NA_integer_

lost_peaks <- union_components[union_components$class == "LOST"]
gained_peaks <- union_components[union_components$class == "GAINED"]
maintained_peaks <- union_components[union_components$class == "MAINTAINED"]

stopifnot(
  length(union_components) == 45021L,
  length(maintained_peaks) == 28964L,
  length(lost_peaks) == 6666L,
  length(gained_peaks) == 9391L,
  all(GenomicRanges::countOverlaps(union_components, union_components) == 1L)
)

chipseq_msg(
  "True union: ", length(union_components), " non-overlapping components (",
  length(maintained_peaks), " maintained, ", length(lost_peaks), " lost, ",
  length(gained_peaks), " gained)."
)

# Preserve the record-level overlap calls only as an audit layer. findOverlaps
# returns all pairs, but unique query/subject indices prevent pair expansion.
# These record counts must not be described as a genomic-site union.
chipseq_msg("Computing record-level and reciprocal-overlap audit summaries.")
ov_1bp_raw <- GenomicRanges::findOverlaps(wt_peaks, mut_peaks, minoverlap = 1L, ignore.strand = TRUE)
wt_record_class <- rep("LOST", length(wt_peaks))
wt_record_class[unique(S4Vectors::queryHits(ov_1bp_raw))] <- "MAINTAINED"
mut_record_class <- rep("GAINED", length(mut_peaks))
mut_record_class[unique(S4Vectors::subjectHits(ov_1bp_raw))] <- "MAINTAINED"
wt_peaks$record_class <- wt_record_class
mut_peaks$record_class <- mut_record_class

ov_1bp_unique <- GenomicRanges::findOverlaps(wt_unique, mut_unique, minoverlap = 1L, ignore.strand = TRUE)
ov_df <- find_best_overlap_pairs(wt_unique, mut_unique, ov_1bp_unique)
if (nrow(ov_df) > 0L) {
  ov_df$wt_width <- width(wt_unique)[ov_df$query]
  ov_df$mut_width <- width(mut_unique)[ov_df$subject]
  ov_df$frac_wt <- ov_df$overlap_width / ov_df$wt_width
  ov_df$frac_mut <- ov_df$overlap_width / ov_df$mut_width
  ov_df$reciprocal_50 <- ov_df$frac_wt >= 0.5 & ov_df$frac_mut >= 0.5
} else {
  ov_df <- data.frame()
}

wt_maintained_1bp <- unique(S4Vectors::queryHits(ov_1bp_unique))
mut_maintained_1bp <- unique(S4Vectors::subjectHits(ov_1bp_unique))
wt_maintained_50 <- unique(ov_df$query[ov_df$reciprocal_50])
mut_maintained_50 <- unique(ov_df$subject[ov_df$reciprocal_50])

sensitivity_summary <- dplyr::bind_rows(
  data.frame(
    overlap_rule = "1bp_any_overlap_unique_source_intervals",
    wt_lost = length(wt_unique) - length(wt_maintained_1bp),
    wt_maintained = length(wt_maintained_1bp),
    mut_gained = length(mut_unique) - length(mut_maintained_1bp),
    mut_maintained = length(mut_maintained_1bp)
  ),
  data.frame(
    overlap_rule = "50pct_reciprocal_overlap_unique_source_intervals",
    wt_lost = length(wt_unique) - length(wt_maintained_50),
    wt_maintained = length(wt_maintained_50),
    mut_gained = length(mut_unique) - length(mut_maintained_50),
    mut_maintained = length(mut_maintained_50)
  )
)
chipseq_write_csv(sensitivity_summary, "chipseq_02_peak_overlap_sensitivity_summary.csv")

summary_df <- data.frame(
  class = c("LOST", "MAINTAINED", "GAINED"),
  n_peaks = c(length(lost_peaks), length(maintained_peaks), length(gained_peaks))
) %>%
  dplyr::mutate(
    percent_of_union = n_peaks / sum(n_peaks) * 100,
    class = factor(class, levels = c("MAINTAINED", "LOST", "GAINED"))
  )
chipseq_write_csv(summary_df, "chipseq_02_peak_classification_summary.csv")

reconciliation_df <- data.frame(
  unit = c(
    "MACS2 peak records in asymmetric WT-plus-mutant-only catalogue",
    "exact-coordinate-deduplicated asymmetric catalogue",
    "reduced symmetric non-overlapping union components"
  ),
  total = c(
    sum(wt_record_class %in% c("LOST", "MAINTAINED")) + sum(mut_record_class == "GAINED"),
    length(unique(c(
      paste0(seqnames(wt_peaks), ":", start(wt_peaks), "-", end(wt_peaks)),
      paste0(seqnames(mut_peaks[mut_record_class == "GAINED"]), ":", start(mut_peaks[mut_record_class == "GAINED"]), "-", end(mut_peaks[mut_record_class == "GAINED"]))
    ))),
    length(union_components)
  ),
  maintained = c(sum(wt_record_class == "MAINTAINED"), length(wt_maintained_1bp), length(maintained_peaks)),
  lost = c(sum(wt_record_class == "LOST"), length(wt_unique) - length(wt_maintained_1bp), length(lost_peaks)),
  gained = c(sum(mut_record_class == "GAINED"), length(mut_unique) - length(mut_maintained_1bp), length(gained_peaks))
)
chipseq_write_csv(reconciliation_df, "chipseq_02_peak_unit_reconciliation.csv")

chipseq_save_rds(list(
  wt_peaks = wt_peaks,
  mut_peaks = mut_peaks,
  wt_unique_peaks = wt_unique,
  mut_unique_peaks = mut_unique,
  union_components = union_components,
  lost_peaks = lost_peaks,
  gained_peaks = gained_peaks,
  maintained_peaks = maintained_peaks,
  raw_overlaps = ov_df
), "chipseq_02_classified_peaks.rds")

export_union_component_bed <- function(gr, fname) {
  bed <- data.frame(
    chrom = as.character(seqnames(gr)),
    start = start(gr) - 1L,
    end = end(gr),
    union_id = gr$union_id,
    score = 0,
    strand = ".",
    peak_class = gr$class,
    n_wt_unique_intervals = gr$n_wt_unique_intervals,
    n_mut_unique_intervals = gr$n_mut_unique_intervals,
    n_wt_peak_records = gr$n_wt_peak_records,
    n_mut_peak_records = gr$n_mut_peak_records,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    bed,
    file = file.path(CHIP_DIRS$tables, fname),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  invisible(bed)
}

export_union_component_bed(union_components, "chipseq_02_all_peaks_classified.bed")
export_union_component_bed(union_components, "chipseq_02_union_components_classified.bed")
export_peak_bed(lost_peaks, "chipseq_CTCF_WT_only_LOST_peaks.bed")
export_peak_bed(gained_peaks, "chipseq_CTCF_MUT_only_GAINED_peaks.bed")
export_peak_bed(maintained_peaks, "chipseq_CTCF_MAINTAINED_peaks.bed")

p_bar <- ggplot(summary_df, aes(x = class, y = n_peaks, fill = class)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = scales::comma(n_peaks)), vjust = -0.25, fontface = "bold") +
  scale_fill_manual(values = c(
    "LOST" = chipseq_palette$lost,
    "MAINTAINED" = chipseq_palette$maintained,
    "GAINED" = chipseq_palette$gained
  )) +
  scale_y_continuous(labels = scales::comma_format(), expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "CTCF union-component classification",
    subtitle = "Exact duplicate records collapsed; overlapping WT and mutant intervals reduced to disjoint components",
    x = NULL,
    y = "Union-component count"
  )
chipseq_save_plot(p_bar, "chipseq_02_peak_classification_bar", subdir = "peak_overview", w = 7.2, h = 5.4)

p_stack <- ggplot(summary_df, aes(x = "Classified components", y = percent_of_union, fill = class)) +
  geom_col(width = 0.34, color = "white", linewidth = 0.55, position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent_of_union)),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    color = "white",
    fontface = "bold",
    size = 4.2
  ) +
  scale_fill_manual(values = c(
    "LOST" = chipseq_palette$lost,
    "MAINTAINED" = chipseq_palette$maintained,
    "GAINED" = chipseq_palette$gained
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "Peak-union component proportions",
    subtitle = "45,021 disjoint genomic components",
    x = NULL,
    y = "Percent of components"
  )
chipseq_save_plot(p_stack, "chipseq_02_peak_classification_stacked", subdir = "peak_overview", w = 4.4, h = 4.6)

if (nrow(ov_df) > 0L) {
  maintained_pairs <- ov_df %>%
    dplyr::mutate(
      wt_signal = wt_unique$signalValue[query],
      mut_signal = mut_unique$signalValue[subject],
      wt_q = wt_unique$qValue[query],
      mut_q = mut_unique$qValue[subject],
      wt_name = wt_unique$name[query],
      mut_name = mut_unique$name[subject]
    )
  chipseq_write_csv(maintained_pairs, "chipseq_02_maintained_peak_signal_pairs.csv")
}

chipseq_msg("Script 02 complete: true union components, audit summaries, BED exports and plots written.")
