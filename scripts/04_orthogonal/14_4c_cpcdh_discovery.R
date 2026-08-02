# ============================================================
# 14_4c_cpcdh_discovery.R — replicate-aware cPcdh locus discovery analysis from 4C bigWigs
# ============================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])[1]
file_arg <- gsub("~\\+~", " ", file_arg)
SCRIPT_DIR <- dirname(normalizePath(ifelse(length(file_arg) == 0 || is.na(file_arg), ".", file_arg)))
source(file.path(SCRIPT_DIR, "12_discovery_setup.R"))
discovery_require(c(
  "AnnotationDbi", "dplyr", "GenomicRanges", "ggplot2", "patchwork",
  "readr", "rtracklayer", "scales", "stringr", "tidyr",
  "org.Mm.eg.db", "TxDb.Mmusculus.UCSC.mm10.knownGene"
))
discovery_append_log("fourc", "14_4c_cpcdh_discovery")

pcdh_symbols <- sort(unique(AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "SYMBOL")[grepl("^Pcdh(a|b|g)", AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "SYMBOL"))]))
pcdh_gr <- get_mm10_gene_ranges(pcdh_symbols)
if (length(pcdh_gr) == 0) {
  stop("Could not map cPcdh genes in mm10. Cannot run 4C discovery analysis.", call. = FALSE)
}

pcdh_df <- as.data.frame(pcdh_gr) %>%
  dplyr::mutate(
    gene = pcdh_gr$gene_symbol,
    family = dplyr::case_when(
      grepl("^Pcdha", gene) ~ "Pcdha",
      grepl("^Pcdhb", gene) ~ "Pcdhb",
      grepl("^Pcdhg", gene) ~ "Pcdhg",
      TRUE ~ "Other"
    )
  )

family_ranges <- pcdh_df %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(
    seqnames = dplyr::first(seqnames),
    start = min(start),
    end = max(end),
    .groups = "drop"
  )

cluster_chr <- unique(family_ranges$seqnames)
if (length(cluster_chr) != 1) {
  stop("cPcdh family genes mapped to multiple chromosomes; expected one locus.", call. = FALSE)
}

bin_size <- 10000L
cluster_start <- max(1L, floor((min(family_ranges$start) - 250000L) / bin_size) * bin_size)
cluster_end <- ceiling((max(family_ranges$end) + 250000L) / bin_size) * bin_size
cluster_region <- GenomicRanges::GRanges(cluster_chr, IRanges::IRanges(cluster_start, cluster_end))
n_bins <- as.integer(ceiling(GenomicRanges::width(cluster_region) / bin_size))

fourc_files <- list.files(FOURC_RAW, pattern = "\\.bw$", full.names = TRUE)
file_meta <- data.frame(file = fourc_files, stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    basename = basename(file),
    genotype = dplyr::case_when(
      stringr::str_detect(basename, "_homo_") ~ "HOMO",
      stringr::str_detect(basename, "_wt_") ~ "WT",
      TRUE ~ NA_character_
    ),
    viewpoint = toupper(stringr::str_match(basename, "_(bc|f)_rep")[, 2]),
    replicate = paste0("rep", stringr::str_match(basename, "_rep(\\d+)\\.bw$")[, 2]),
    file_id = tools::file_path_sans_ext(basename)
  ) %>%
  dplyr::filter(!is.na(genotype), !is.na(viewpoint))

write_csv_discovery("fourc", file_meta, "14_4c_file_manifest.csv")

extract_profile <- function(file_path) {
  gr <- rtracklayer::summary(
    rtracklayer::BigWigFile(file_path),
    which = cluster_region,
    size = n_bins,
    type = "mean"
  )[[1]]
  data.frame(
    chr = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    bin_id = seq_along(gr),
    signal = gr$score,
    stringsAsFactors = FALSE
  )
}

profiles <- dplyr::bind_rows(lapply(seq_len(nrow(file_meta)), function(i) {
  cbind(file_meta[i, , drop = FALSE], extract_profile(file_meta$file[i]))
})) %>%
  dplyr::mutate(
    signal = dplyr::coalesce(signal, 0),
    log_signal = log1p(signal),
    position_mb = ((start + end) / 2) / 1e6
  ) %>%
  dplyr::group_by(viewpoint, file_id) %>%
  dplyr::mutate(z_signal = as.numeric(scale(log_signal))) %>%
  dplyr::ungroup()

write_csv_discovery("fourc", profiles, "14_4c_binned_profiles.csv")

replicate_correlations <- dplyr::bind_rows(lapply(unique(file_meta$viewpoint), function(vp) {
  dplyr::bind_rows(lapply(c("WT", "HOMO"), function(gt) {
    gt_files <- file_meta %>% dplyr::filter(viewpoint == vp, genotype == gt) %>% dplyr::pull(file_id)
    if (length(gt_files) < 2) {
      return(data.frame(viewpoint = vp, genotype = gt, spearman_cor = NA_real_, stringsAsFactors = FALSE))
    }
    wide <- profiles %>%
      dplyr::filter(viewpoint == vp, file_id %in% gt_files) %>%
      dplyr::select(bin_id, file_id, log_signal) %>%
      tidyr::pivot_wider(names_from = file_id, values_from = log_signal)
    data.frame(
      viewpoint = vp,
      genotype = gt,
      spearman_cor = safe_cor(wide[[2]], wide[[3]], method = "spearman"),
      stringsAsFactors = FALSE
    )
  }))
}))
write_csv_discovery("fourc", replicate_correlations, "14_4c_replicate_correlations.csv")

bin_wide <- profiles %>%
  dplyr::mutate(track_id = paste(genotype, replicate, sep = "_")) %>%
  dplyr::select(viewpoint, bin_id, chr, start, end, position_mb, track_id, log_signal) %>%
  tidyr::pivot_wider(names_from = track_id, values_from = log_signal)

bin_thresholds <- bin_wide %>%
  dplyr::group_by(viewpoint) %>%
  dplyr::mutate(mean_delta = rowMeans(cbind(HOMO_rep1, HOMO_rep2), na.rm = TRUE) - rowMeans(cbind(WT_rep1, WT_rep2), na.rm = TRUE)) %>%
  dplyr::summarise(delta_threshold = max(0.25, stats::quantile(abs(mean_delta), probs = 0.9, na.rm = TRUE) * 0.5), .groups = "drop")

bin_wide <- bin_wide %>%
  dplyr::left_join(bin_thresholds, by = "viewpoint") %>%
  dplyr::mutate(
    wt_min = pmin(WT_rep1, WT_rep2, na.rm = TRUE),
    wt_max = pmax(WT_rep1, WT_rep2, na.rm = TRUE),
    homo_min = pmin(HOMO_rep1, HOMO_rep2, na.rm = TRUE),
    homo_max = pmax(HOMO_rep1, HOMO_rep2, na.rm = TRUE),
    mean_delta = rowMeans(cbind(HOMO_rep1, HOMO_rep2), na.rm = TRUE) - rowMeans(cbind(WT_rep1, WT_rep2), na.rm = TRUE),
    direction = dplyr::case_when(
      homo_min > wt_max + 0.15 & mean_delta > delta_threshold ~ "higher_in_homo",
      homo_max < wt_min - 0.15 & mean_delta < -delta_threshold ~ "higher_in_wt",
      TRUE ~ "stable"
    )
  )

merge_altered_intervals <- function(df) {
  df <- df %>% dplyr::filter(direction != "stable") %>% dplyr::arrange(start)
  if (nrow(df) == 0) {
    return(data.frame())
  }
  breaks <- c(TRUE, (diff(df$start) > bin_size) | (df$direction[-1] != df$direction[-nrow(df)]))
  df$interval_group <- cumsum(breaks)
  df %>%
    dplyr::group_by(viewpoint, interval_group, direction) %>%
    dplyr::summarise(
      chr = dplyr::first(chr),
      start = min(start),
      end = max(end),
      n_bins = dplyr::n(),
      mean_delta = mean(mean_delta, na.rm = TRUE),
      peak_abs_delta = max(abs(mean_delta), na.rm = TRUE),
      .groups = "drop"
    )
}

altered_intervals <- dplyr::bind_rows(lapply(split(bin_wide, bin_wide$viewpoint), merge_altered_intervals))

family_gr <- GenomicRanges::GRanges(
  seqnames = family_ranges$seqnames,
  ranges = IRanges::IRanges(family_ranges$start, family_ranges$end)
)
family_gr$family <- family_ranges$family

if (nrow(altered_intervals) > 0) {
  interval_gr <- GenomicRanges::GRanges(
    seqnames = altered_intervals$chr,
    ranges = IRanges::IRanges(altered_intervals$start, altered_intervals$end)
  )
  hits <- GenomicRanges::findOverlaps(interval_gr, family_gr, ignore.strand = TRUE)
  interval_family_map <- data.frame(
    interval_group = altered_intervals$interval_group[queryHits(hits)],
    family = family_gr$family[subjectHits(hits)],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(interval_group) %>%
    dplyr::summarise(family_overlap = paste(sort(unique(family)), collapse = ";"), .groups = "drop")
  altered_intervals <- altered_intervals %>%
    dplyr::left_join(interval_family_map, by = "interval_group") %>%
    dplyr::mutate(family_overlap = dplyr::coalesce(family_overlap, "inter-family"))
}

write_csv_discovery("fourc", bin_wide, "14_4c_bin_deltas.csv")
write_csv_discovery("fourc", altered_intervals, "14_4c_altered_intervals.csv")

bulk_candidate_path <- file.path(discovery_dirs_for("bulk")$tables, "13_brain_candidate_gene_dosage_summary.csv")
bulk_cpcdh_summary <- if (file.exists(bulk_candidate_path)) {
  readr::read_csv(bulk_candidate_path, show_col_types = FALSE) %>%
    dplyr::filter(grepl("^Pcdh(a|b|g)", gene)) %>%
    dplyr::mutate(
      family = dplyr::case_when(
        grepl("^Pcdha", gene) ~ "Pcdha",
        grepl("^Pcdhb", gene) ~ "Pcdhb",
        grepl("^Pcdhg", gene) ~ "Pcdhg",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::group_by(family) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      bulk_mean_HOMO_vs_WT_logFC = mean(HOMO_vs_WT_logFC, na.rm = TRUE),
      dominant_pattern = paste(sort(unique(trend_class)), collapse = "; "),
      .groups = "drop"
    )
} else {
  data.frame()
}
write_csv_discovery("fourc", bulk_cpcdh_summary, "14_4c_bulk_cpcdh_integration.csv")

mean_profiles <- profiles %>%
  dplyr::group_by(viewpoint, genotype, bin_id, position_mb, start, end) %>%
  dplyr::summarise(mean_log_signal = mean(log_signal, na.rm = TRUE), .groups = "drop")

p_track <- ggplot() +
  geom_rect(
    data = family_ranges,
    aes(
      xmin = start / 1e6,
      xmax = end / 1e6,
      ymin = -Inf,
      ymax = Inf,
      fill = family
    ),
    inherit.aes = FALSE,
    alpha = 0.08
  ) +
  geom_line(
    data = profiles,
    aes(x = position_mb, y = log_signal, group = file_id, color = genotype),
    linewidth = 0.5,
    alpha = 0.35
  ) +
  geom_line(
    data = mean_profiles,
    aes(x = position_mb, y = mean_log_signal, color = genotype),
    linewidth = 1.15
  ) +
  facet_wrap(~viewpoint, ncol = 1, scales = "free_y") +
  scale_color_manual(values = DISCOVERY_COLORS[c("WT", "HOMO")]) +
  scale_fill_manual(values = c(Pcdha = "#8ecae6", Pcdhb = "#ffb703", Pcdhg = "#90be6d", Other = "#dddddd")) +
  theme_classic(base_size = 12) +
  labs(
    title = "cPcdh 4C profiles across WT and HOMO E18.5 brain",
    subtitle = "Thin lines are individual replicates; thick lines are genotype means. Shaded bands mark cPcdh family sub-clusters.",
    x = paste0(cluster_chr, " position (Mb, mm10)"),
    y = "log1p(4C signal)",
    color = "Genotype",
    fill = "Family"
  )
save_plot_discovery("fourc", p_track, "14_4c_cpcdh_structural_map", w = 11.5, h = 7.6)

if (nrow(altered_intervals) > 0) {
  p_intervals <- ggplot(altered_intervals, aes(x = start / 1e6, xend = end / 1e6, y = viewpoint, yend = viewpoint, color = direction)) +
    geom_segment(linewidth = 4, alpha = 0.9) +
    geom_text(aes(x = (start + end) / 2e6, label = family_overlap), vjust = -0.9, size = 3.1, color = "black") +
    scale_color_manual(values = c(higher_in_homo = "#d95f02", higher_in_wt = "#1f77b4")) +
    theme_classic(base_size = 12) +
    labs(
      title = "Reproducible altered 4C intervals within the cPcdh locus",
      subtitle = "Intervals require consistent direction across both replicates within each viewpoint.",
      x = paste0(cluster_chr, " position (Mb, mm10)"),
      y = "Viewpoint",
      color = "Direction"
    )
  save_plot_discovery("fourc", p_intervals, "14_4c_altered_interval_summary", w = 10.5, h = 4.6)
}

p_corr <- ggplot(replicate_correlations, aes(x = viewpoint, y = spearman_cor, fill = genotype)) +
  geom_col(position = position_dodge(width = 0.6), width = 0.55) +
  geom_text(aes(label = sprintf("%.2f", spearman_cor)), position = position_dodge(width = 0.6), vjust = -0.5, size = 3.4) +
  scale_fill_manual(values = DISCOVERY_COLORS[c("WT", "HOMO")]) +
  coord_cartesian(ylim = c(0, 1.05)) +
  theme_classic(base_size = 12) +
  labs(
    title = "4C replicate concordance by viewpoint and genotype",
    x = "Viewpoint",
    y = "Spearman correlation across 10 kb bins",
    fill = "Genotype"
  )
save_plot_discovery("fourc", p_corr, "14_4c_replicate_concordance", w = 7.2, h = 4.8)

summary_lines <- c(
  "# 4C cPcdh discovery summary",
  "",
  paste0("- Locus window analysed: ", cluster_chr, ":", scales::comma(cluster_start), "-", scales::comma(cluster_end), " (", n_bins, " bins at 10 kb)."),
  paste0("- WT replicate correlations: ", paste(sprintf("%.2f", replicate_correlations$spearman_cor[replicate_correlations$genotype == "WT"]), collapse = ", ")),
  paste0("- HOMO replicate correlations: ", paste(sprintf("%.2f", replicate_correlations$spearman_cor[replicate_correlations$genotype == "HOMO"]), collapse = ", ")),
  paste0("- Reproducible altered intervals called: ", nrow(altered_intervals)),
  if (nrow(altered_intervals) > 0) {
    paste0("- Strongest interval: ", altered_intervals$viewpoint[1], " ", altered_intervals$direction[1], " at ", altered_intervals$chr[1], ":", scales::comma(altered_intervals$start[1]), "-", scales::comma(altered_intervals$end[1]), " (", altered_intervals$family_overlap[1], ").")
  } else {
    "- No intervals passed the cross-replicate consistency threshold."
  }
)
write_status_note("fourc", "14_4c_summary.md", summary_lines)

discovery_msg("4C cPcdh discovery analysis complete.")
