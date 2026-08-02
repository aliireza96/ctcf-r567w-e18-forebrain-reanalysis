# ============================================================
# 02_chipseq_peak_classification.R — classify LOST / GAINED / MAINTAINED peaks
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

chipseq_msg("Classifying peaks using default 1 bp overlap criterion.")
ov_1bp <- GenomicRanges::findOverlaps(wt_peaks, mut_peaks, minoverlap = 1L, ignore.strand = TRUE)

wt_class <- rep("LOST", length(wt_peaks))
wt_class[unique(S4Vectors::queryHits(ov_1bp))] <- "MAINTAINED"
mut_class <- rep("GAINED", length(mut_peaks))
mut_class[unique(S4Vectors::subjectHits(ov_1bp))] <- "MAINTAINED"

wt_peaks$class <- wt_class
mut_peaks$class <- mut_class

lost_peaks <- wt_peaks[wt_peaks$class == "LOST"]
gained_peaks <- mut_peaks[mut_peaks$class == "GAINED"]
maintained_peaks <- wt_peaks[wt_peaks$class == "MAINTAINED"]
maintained_mut_peaks <- mut_peaks[mut_peaks$class == "MAINTAINED"]

chipseq_msg("Computing stricter reciprocal-overlap sensitivity summary.")
ov_df <- find_best_overlap_pairs(wt_peaks, mut_peaks, ov_1bp)
if (nrow(ov_df) > 0) {
  ov_df$wt_width <- width(wt_peaks)[ov_df$query]
  ov_df$mut_width <- width(mut_peaks)[ov_df$subject]
  ov_df$frac_wt <- ov_df$overlap_width / ov_df$wt_width
  ov_df$frac_mut <- ov_df$overlap_width / ov_df$mut_width
  ov_df$reciprocal_50 <- ov_df$frac_wt >= 0.5 & ov_df$frac_mut >= 0.5
} else {
  ov_df <- data.frame()
}

wt_maintained_50 <- unique(ov_df$query[ov_df$reciprocal_50])
mut_maintained_50 <- unique(ov_df$subject[ov_df$reciprocal_50])

sensitivity_summary <- dplyr::bind_rows(
  data.frame(
    overlap_rule = "1bp_any_overlap",
    wt_lost = length(lost_peaks),
    wt_maintained = length(maintained_peaks),
    mut_gained = length(gained_peaks),
    mut_maintained = length(maintained_mut_peaks)
  ),
  data.frame(
    overlap_rule = "50pct_reciprocal_overlap",
    wt_lost = length(wt_peaks) - length(wt_maintained_50),
    wt_maintained = length(wt_maintained_50),
    mut_gained = length(mut_peaks) - length(mut_maintained_50),
    mut_maintained = length(mut_maintained_50)
  )
)
chipseq_write_csv(sensitivity_summary, "chipseq_02_peak_overlap_sensitivity_summary.csv")

summary_df <- data.frame(
  class = c("LOST", "MAINTAINED", "GAINED"),
  n_peaks = c(length(lost_peaks), length(maintained_peaks), length(gained_peaks))
) %>%
  dplyr::mutate(
    percent_of_wt_mut_union = n_peaks / sum(n_peaks) * 100,
    class = factor(class, levels = c("MAINTAINED", "LOST", "GAINED"))
  )
chipseq_write_csv(summary_df, "chipseq_02_peak_classification_summary.csv")

chipseq_save_rds(list(
  wt_peaks = wt_peaks,
  mut_peaks = mut_peaks,
  lost_peaks = lost_peaks,
  gained_peaks = gained_peaks,
  maintained_peaks = maintained_peaks,
  maintained_mut_peaks = maintained_mut_peaks,
  overlaps = ov_df
), "chipseq_02_classified_peaks.rds")

export_peak_bed(lost_peaks, "chipseq_CTCF_WT_only_LOST_peaks.bed")
export_peak_bed(gained_peaks, "chipseq_CTCF_MUT_only_GAINED_peaks.bed")
export_peak_bed(maintained_peaks, "chipseq_CTCF_MAINTAINED_peaks.bed")
export_peak_bed(wt_peaks, "chipseq_CTCF_all_WT_peaks_classified.bed")

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
    title = "CTCF peak classification by overlap",
    subtitle = "Default classification uses any WT/MUT overlap (>=1 bp); sensitivity table reports a stricter 50% reciprocal overlap alternative",
    x = NULL,
    y = "Peak count"
  )
chipseq_save_plot(p_bar, "chipseq_02_peak_classification_bar", subdir = "peak_overview", w = 7.2, h = 5.4)

p_stack <- ggplot(summary_df, aes(x = "Classified peaks", y = percent_of_wt_mut_union, fill = class)) +
  geom_col(width = 0.34, color = "white", linewidth = 0.55, position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = sprintf("%.1f%%", percent_of_wt_mut_union)),
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
    plot.title = element_text(face = "bold", color = chipseq_palette$navy, size = 15),
    plot.subtitle = element_text(color = "#4B5563", size = 10.5),
    axis.text.y = element_text(color = chipseq_palette$navy),
    axis.title.y = element_text(color = chipseq_palette$navy),
    axis.line = element_line(color = chipseq_palette$navy, linewidth = 0.45),
    axis.ticks = element_line(color = chipseq_palette$navy, linewidth = 0.35),
    legend.title = element_blank(),
    legend.text = element_text(color = chipseq_palette$navy, size = 10.5),
    legend.key.size = grid::unit(0.45, "cm"),
    legend.position = "right"
  ) +
  labs(
    title = "Peak class proportions",
    subtitle = ">=1 bp WT/MUT overlap rule",
    x = NULL,
    y = "Percent of classified peaks"
  )
chipseq_save_plot(p_stack, "chipseq_02_peak_classification_stacked", subdir = "peak_overview", w = 4.4, h = 4.6)

maintained_pairs <- ov_df
if (nrow(maintained_pairs) > 0) {
  maintained_pairs <- maintained_pairs %>%
    dplyr::mutate(
      wt_signal = wt_peaks$signalValue[query],
      mut_signal = mut_peaks$signalValue[subject],
      wt_q = wt_peaks$qValue[query],
      mut_q = mut_peaks$qValue[subject],
      wt_name = wt_peaks$name[query],
      mut_name = mut_peaks$name[subject]
    )
  chipseq_write_csv(maintained_pairs, "chipseq_02_maintained_peak_signal_pairs.csv")

  cor_val <- suppressWarnings(stats::cor(maintained_pairs$wt_signal, maintained_pairs$mut_signal, use = "pairwise.complete.obs"))
  p_scatter <- ggplot(maintained_pairs, aes(x = wt_signal, y = mut_signal)) +
    geom_point(alpha = 0.35, color = chipseq_palette$navy, size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = chipseq_palette$coral) +
    geom_smooth(method = "lm", se = FALSE, color = chipseq_palette$teal, linewidth = 0.9) +
    theme_classic(base_size = 13) +
    annotate("text", x = Inf, y = -Inf, label = paste0("Pearson r = ", round(cor_val, 3)), hjust = 1.05, vjust = -0.6, size = 4.2, fontface = "bold") +
    labs(
      title = "Signal comparison at maintained CTCF peaks",
      subtitle = "Each point pairs a WT maintained peak with its best-overlap mutant peak",
      x = "WT signalValue",
      y = "MUT signalValue"
    )
  chipseq_save_plot(p_scatter, "chipseq_02_maintained_peak_signal_scatter", subdir = "peak_overview", w = 6.4, h = 5.8)
}

chipseq_msg("Script 02 complete: classified peak sets, BED exports, and overlap summaries written.")
