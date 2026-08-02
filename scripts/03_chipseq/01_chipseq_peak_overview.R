# ============================================================
# 01_chipseq_peak_overview.R — load peaks and global overview
# ============================================================

source(file.path("scripts", "chipseq", "00_chipseq_setup.R"))
chipseq_require(c("GenomicRanges", "rtracklayer", "ggplot2", "dplyr", "readr", "tidyr", "scales"))

chipseq_msg("Loading CTCF narrowPeak files.")
wt_peaks <- standardize_peak_gr(load_narrowpeak(CHIP_FILES$wt_peaks), sample_label = "WT")
mut_peaks <- standardize_peak_gr(load_narrowpeak(CHIP_FILES$mut_peaks), sample_label = "MUT")

chipseq_save_rds(wt_peaks, "chipseq_01_wt_peaks.rds")
chipseq_save_rds(mut_peaks, "chipseq_01_mut_peaks.rds")
chipseq_save_rds(list(wt = wt_peaks, mut = mut_peaks), "chipseq_01_peak_objects.rds")

peak_stats <- dplyr::bind_rows(
  granges_to_df(wt_peaks) %>% dplyr::mutate(sample = "WT"),
  granges_to_df(mut_peaks) %>% dplyr::mutate(sample = "MUT")
) %>%
  dplyr::mutate(
    width = end - start + 1,
    neg_log10_q = ifelse(!is.na(qValue), qValue / 10, NA_real_)
  )

summary_stats <- peak_stats %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    n_peaks = dplyr::n(),
    median_width = stats::median(width, na.rm = TRUE),
    median_signalValue = stats::median(signalValue, na.rm = TRUE),
    mean_signalValue = mean(signalValue, na.rm = TRUE),
    median_qValue = stats::median(qValue, na.rm = TRUE),
    median_neg_log10_q = stats::median(neg_log10_q, na.rm = TRUE),
    mean_neg_log10_q = mean(neg_log10_q, na.rm = TRUE),
    .groups = "drop"
  )
chipseq_write_csv(summary_stats, "chipseq_01_peak_summary_stats.csv")

metric_long <- peak_stats %>%
  dplyr::select(sample, width, signalValue, neg_log10_q) %>%
  tidyr::pivot_longer(cols = c(width, signalValue, neg_log10_q), names_to = "metric", values_to = "value")

p_violin <- ggplot(metric_long, aes(x = sample, y = value, fill = sample)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", alpha = 0.9) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("WT" = chipseq_palette$wt, "MUT" = chipseq_palette$mut)) +
  theme_classic(base_size = 13) +
  labs(
    title = "Genome-wide CTCF peak property overview",
    subtitle = "Peak width, signal, and q-value distributions across WT and R567W homozygous brain ChIP-seq peaks",
    x = NULL,
    y = NULL
  )
chipseq_save_plot(p_violin, "chipseq_01_peak_metric_violin", subdir = "peak_overview", w = 11, h = 4.8)

p_hist_signal <- ggplot(peak_stats, aes(x = signalValue, color = sample, fill = sample)) +
  geom_histogram(position = "identity", bins = 80, alpha = 0.28) +
  scale_color_manual(values = c("WT" = chipseq_palette$wt, "MUT" = chipseq_palette$mut)) +
  scale_fill_manual(values = c("WT" = chipseq_palette$wt, "MUT" = chipseq_palette$mut)) +
  theme_classic(base_size = 13) +
  labs(
    title = "CTCF peak signal distributions",
    subtitle = "signalValue from MACS2 narrowPeak output",
    x = "signalValue",
    y = "Peak count"
  )
chipseq_save_plot(p_hist_signal, "chipseq_01_signal_histogram", subdir = "peak_overview", w = 8, h = 5)

p_hist_q <- ggplot(peak_stats, aes(x = neg_log10_q, color = sample, fill = sample)) +
  geom_histogram(position = "identity", bins = 80, alpha = 0.28) +
  scale_color_manual(values = c("WT" = chipseq_palette$wt, "MUT" = chipseq_palette$mut)) +
  scale_fill_manual(values = c("WT" = chipseq_palette$wt, "MUT" = chipseq_palette$mut)) +
  theme_classic(base_size = 13) +
  labs(
    title = "CTCF peak confidence distributions",
    subtitle = "qValue converted to -log10 scale (MACS2 narrowPeak convention)",
    x = "-log10(q)",
    y = "Peak count"
  )
chipseq_save_plot(p_hist_q, "chipseq_01_qvalue_histogram", subdir = "peak_overview", w = 8, h = 5)

counts_df <- summary_stats %>%
  dplyr::select(sample, n_peaks) %>%
  dplyr::mutate(sample = factor(sample, levels = c("WT", "MUT"))) %>%
  dplyr::arrange(sample)
p_counts <- ggplot(counts_df, aes(x = sample, y = n_peaks, fill = sample)) +
  geom_col(width = 0.46, color = "white", linewidth = 0.45) +
  geom_text(aes(label = scales::comma(n_peaks)), vjust = -0.45, fontface = "bold", size = 4.6, color = chipseq_palette$navy) +
  scale_fill_manual(values = c("WT" = chipseq_palette$wt, "MUT" = chipseq_palette$mut)) +
  scale_y_continuous(labels = scales::comma_format(), expand = expansion(mult = c(0, 0.14))) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", color = chipseq_palette$navy, size = 15),
    plot.subtitle = element_text(color = "#4B5563", size = 10.5),
    axis.text.x = element_text(face = "bold", color = chipseq_palette$navy, size = 12),
    axis.text.y = element_text(color = chipseq_palette$navy),
    axis.title.y = element_text(color = chipseq_palette$navy),
    axis.line = element_line(color = chipseq_palette$navy, linewidth = 0.45),
    axis.ticks = element_line(color = chipseq_palette$navy, linewidth = 0.35)
  ) +
  labs(
    title = "Total CTCF peaks",
    subtitle = "MACS2 narrowPeak calls in E18.5 brain",
    x = NULL,
    y = "Peak count"
  )
chipseq_save_plot(p_counts, "chipseq_01_total_peak_counts", subdir = "peak_overview", w = 4.4, h = 4.6)

chipseq_msg("Script 01 complete: peak overview objects, tables, and plots written to results/chipseq.")
