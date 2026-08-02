#!/usr/bin/env Rscript

# ============================================================
# 09j_pseudotime_variableN_sensitivity.R
# Variable-N balanced subsampling sensitivity for lineage
# pseudotime shifts:
# - delta median pseudotime (MUT - WT)
# - delta late-fraction (MUT - WT)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

source(switch(
  Sys.getenv("PROJECTION_CONFIG", "default"),
  bandler_ruan    = "scripts/09_config_projection_bandler_ruan.R",
  bandler_telley  = "scripts/09_config_projection_bandler_telley.R",
  bandler_dibella = "scripts/09_config_projection_bandler_dibella.R",
  mayer_dibella   = "scripts/09_config_projection_mayer_dibella.R",
  rooted_v2       = "scripts/09_config_projection_rooted_v2.R",
  "scripts/09_config_projection.R"
))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

save_pdf <- function(plot_obj, path_no_ext, w = 8.5, h = 6.0) {
  ggsave(paste0(path_no_ext, ".pdf"), plot_obj, width = w, height = h, useDingbats = FALSE)
}

canon_genotype <- function(x, wt_values, mut_values) {
  x <- as.character(x)
  x_low <- tolower(trimws(x))
  wt_low <- tolower(trimws(as.character(wt_values)))
  mut_low <- tolower(trimws(as.character(mut_values)))
  out <- rep(NA_character_, length(x))
  out[x_low %in% wt_low] <- "WT"
  out[x_low %in% mut_low] <- "MUT"
  out[is.na(out) & grepl("(^wt$|wild|control|ctrl)", x_low)] <- "WT"
  out[is.na(out) & grepl("(mut|r567w|homo|hom|ko|case)", x_low)] <- "MUT"
  out
}

pick_genotype_col <- function(md, preferred) {
  nms <- colnames(md)
  if (!is.null(preferred) && preferred %in% nms) return(preferred)
  fallbacks <- c("condition", "Condition", "group", "Group", "genotype", "Genotype", "genotype2")
  hit <- fallbacks[fallbacks %in% nms]
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

run_variable_n_pt <- function(df, n_values, n_iter = 50L, late_q = 0.80, seed = 1L) {
  set.seed(seed)
  df <- copy(df)[!is.na(genotype2) & is.finite(pseudotime_transfer)]
  wt <- df[genotype2 == "WT"]
  mut <- df[genotype2 == "MUT"]
  out <- vector("list", length(n_values))

  for (i in seq_along(n_values)) {
    n_equal <- n_values[i]
    iter_list <- vector("list", n_iter)
    for (b in seq_len(n_iter)) {
      wi <- wt[sample(.N, n_equal, replace = FALSE)]
      mi <- mut[sample(.N, n_equal, replace = FALSE)]
      sub <- rbindlist(list(wi, mi), fill = TRUE)
      late_thr <- as.numeric(quantile(sub$pseudotime_transfer, late_q, na.rm = TRUE))
      iter_list[[b]] <- data.table(
        N = n_equal,
        iter = b,
        delta_median = median(mi$pseudotime_transfer, na.rm = TRUE) - median(wi$pseudotime_transfer, na.rm = TRUE),
        delta_latefrac = mean(mi$pseudotime_transfer >= late_thr, na.rm = TRUE) - mean(wi$pseudotime_transfer >= late_thr, na.rm = TRUE)
      )
    }
    out[[i]] <- rbindlist(iter_list)
  }

  rbindlist(out)
}

run_variable_n_bins <- function(df, n_values, n_iter = 50L, n_bins = 10L, seed = 1L) {
  set.seed(seed)
  df <- copy(df)[!is.na(genotype2) & is.finite(pseudotime_transfer)]
  wt <- df[genotype2 == "WT"]
  mut <- df[genotype2 == "MUT"]
  out <- vector("list", length(n_values))

  for (i in seq_along(n_values)) {
    n_equal <- n_values[i]
    iter_list <- vector("list", n_iter)
    for (b in seq_len(n_iter)) {
      wi <- wt[sample(.N, n_equal, replace = FALSE)]
      mi <- mut[sample(.N, n_equal, replace = FALSE)]
      sub <- rbindlist(list(wi, mi), fill = TRUE)
      probs <- seq(0, 1, length.out = n_bins + 1L)
      brks <- as.numeric(quantile(sub$pseudotime_transfer, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
      brks <- unique(brks)
      if (length(brks) < 3L) next

      sub[, pt_bin := cut(
        pseudotime_transfer,
        breaks = brks,
        include.lowest = TRUE,
        labels = paste0("D", sprintf("%02d", seq_len(length(brks) - 1L)))
      )]
      occ <- sub[, .(n_cells = .N), by = .(genotype2, pt_bin)]
      occ[, frac := n_cells / sum(n_cells), by = genotype2]
      occ_w <- dcast(occ, pt_bin ~ genotype2, value.var = "frac", fill = 0)
      if (!"WT" %in% colnames(occ_w)) occ_w[, WT := 0]
      if (!"MUT" %in% colnames(occ_w)) occ_w[, MUT := 0]
      occ_delta <- occ_w[, .(
        N = n_equal,
        iter = b,
        pt_bin,
        delta_frac = MUT - WT
      )]
      occ_geno <- occ[, .(N = n_equal, iter = b, genotype2, pt_bin, frac)]
      iter_list[[b]] <- merge(
        occ_geno,
        occ_delta,
        by = c("N", "iter", "pt_bin"),
        all = TRUE,
        allow.cartesian = FALSE
      )
    }
    out[[i]] <- rbindlist(iter_list, fill = TRUE)
  }

  rbindlist(out, fill = TRUE)
}

summarize_metric <- function(df, metric) {
  stopifnot(metric %in% colnames(df))
  as.data.table(df)[, .(
    mean_value = mean(get(metric), na.rm = TRUE),
    ci_low = quantile(get(metric), 0.025, na.rm = TRUE),
    ci_high = quantile(get(metric), 0.975, na.rm = TRUE)
  ), by = N][order(N)]
}

plot_metric_vs_n <- function(summary_tbl, x_col, lo_col, hi_col, title, x_lab) {
  ggplot(summary_tbl, aes(y = N, x = .data[[x_col]])) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = .data[[lo_col]], xmax = .data[[hi_col]]), height = 120, color = "#9ecae1", linewidth = 1.0) +
    geom_path(color = "#2171b5", linewidth = 0.8) +
    geom_point(color = "#08519c", size = 1.8) +
    theme_classic(12) +
    labs(title = title, x = x_lab, y = "Equal cells sampled per condition (N)")
}

plot_lollipop <- function(tbl, value_col, lo_col, hi_col, title, subtitle, x_lab) {
  d <- copy(tbl)
  d[, direction := ifelse(get(value_col) >= 0, "Enriched in mut", "Enriched in wt")]
  d <- d[order(get(value_col))]
  d[, lineage_name := factor(lineage_name, levels = lineage_name)]

  ggplot(d, aes(x = .data[[value_col]], y = lineage_name, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = 0, xend = .data[[value_col]], y = lineage_name, yend = lineage_name), linewidth = 1.0) +
    geom_errorbarh(aes(xmin = .data[[lo_col]], xmax = .data[[hi_col]]), height = 0.18, linewidth = 0.8) +
    geom_point(size = 2.8) +
    scale_color_manual(values = c("Enriched in mut" = "#d95f02", "Enriched in wt" = "#1f77b4")) +
    theme_classic(12) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = NULL, color = NULL)
}

summarize_bins <- function(df) {
  as.data.table(df)[, .(
    mean_value = mean(delta_frac, na.rm = TRUE),
    ci_low = quantile(delta_frac, 0.025, na.rm = TRUE),
    ci_high = quantile(delta_frac, 0.975, na.rm = TRUE)
  ), by = .(N, pt_bin)][order(pt_bin, N)]
}

summarize_bins_by_genotype <- function(df) {
  as.data.table(df)[!is.na(genotype2), .(
    frac_mean = mean(frac, na.rm = TRUE),
    frac_ci_low = quantile(frac, 0.025, na.rm = TRUE),
    frac_ci_high = quantile(frac, 0.975, na.rm = TRUE)
  ), by = .(N, genotype2, pt_bin)][order(genotype2, pt_bin, N)]
}

plot_bins_vs_n <- function(summary_tbl, title) {
  pdat <- copy(summary_tbl)
  pdat[, `:=`(
    mean_value_pct = 100 * mean_value,
    ci_low_pct = 100 * ci_low,
    ci_high_pct = 100 * ci_high
  )]

  ggplot(pdat, aes(y = N, x = mean_value_pct)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = ci_low_pct, xmax = ci_high_pct), height = 120, color = "#9ecae1", linewidth = 0.9) +
    geom_path(color = "#2171b5", linewidth = 0.8) +
    geom_point(color = "#08519c", size = 1.5) +
    facet_wrap(~ pt_bin, ncol = 5, scales = "free_x") +
    theme_classic(11) +
    labs(
      title = title,
      x = "Delta occupancy (percentage points, MUT - WT)",
      y = "Equal cells sampled per condition (N)"
    )
}

plot_bins_pooled <- function(tbl, title, subtitle) {
  d <- copy(tbl)
  d[, direction := ifelse(mean_value_pct >= 0, "Enriched in mut", "Enriched in wt")]
  d <- d[order(mean_value_pct)]
  d[, pt_bin := factor(pt_bin, levels = pt_bin)]

  ggplot(d, aes(x = mean_value_pct, y = pt_bin, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = 0, xend = mean_value_pct, y = pt_bin, yend = pt_bin), linewidth = 0.9) +
    geom_errorbarh(aes(xmin = ci_low_pct, xmax = ci_high_pct), height = 0.18, linewidth = 0.8) +
    geom_point(size = 2.6) +
    scale_color_manual(values = c("Enriched in mut" = "#d95f02", "Enriched in wt" = "#1f77b4")) +
    theme_classic(11) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Average delta occupancy (percentage points, MUT - WT)",
      y = "Pseudotime decile",
      color = NULL
    )
}

plot_occupancy_by_bin_ci <- function(summary_tbl, title, subtitle) {
  pdat <- copy(summary_tbl)
  pdat[, pt_bin_num := as.integer(sub("^D", "", pt_bin))]

  ggplot(pdat, aes(pt_bin_num, frac_mean, color = genotype2, fill = genotype2, group = genotype2)) +
    geom_ribbon(aes(ymin = frac_ci_low, ymax = frac_ci_high), alpha = 0.18, linewidth = 0) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.5) +
    scale_x_continuous(breaks = seq_len(max(pdat$pt_bin_num, na.rm = TRUE))) +
    theme_classic(base_size = 12) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Pseudotime decile (low -> high)",
      y = "Cell fraction"
    )
}

plot_occupancy_delta_by_bin_ci <- function(summary_tbl, title, subtitle) {
  pdat <- copy(summary_tbl)
  pdat[, `:=`(
    mean_value_pct = 100 * mean_value,
    ci_low_pct = 100 * ci_low,
    ci_high_pct = 100 * ci_high,
    pt_bin_num = as.integer(sub("^D", "", pt_bin))
  )]

  ggplot(pdat, aes(pt_bin_num, mean_value_pct, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_ribbon(aes(ymin = ci_low_pct, ymax = ci_high_pct), fill = "#9ecae1", alpha = 0.22, linewidth = 0) +
    geom_line(color = "#2171b5", linewidth = 1.0) +
    geom_point(color = "#08519c", size = 1.6) +
    scale_x_continuous(breaks = seq_len(max(pdat$pt_bin_num, na.rm = TRUE))) +
    theme_classic(base_size = 12) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Pseudotime decile (low -> high)",
      y = "Delta occupancy (percentage points, MUT - WT)"
    )
}

compute_observed_stats <- function(df, late_q = 0.80) {
  df <- copy(df)[!is.na(genotype2) & is.finite(pseudotime_transfer)]
  wt <- df[genotype2 == "WT"]
  mut <- df[genotype2 == "MUT"]
  if (nrow(wt) == 0 || nrow(mut) == 0) return(NULL)

  late_thr <- as.numeric(quantile(df$pseudotime_transfer, late_q, na.rm = TRUE))
  wt_vals <- wt$pseudotime_transfer
  mut_vals <- mut$pseudotime_transfer

  wilcox_p <- tryCatch(
    wilcox.test(mut_vals, wt_vals, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  data.table(
    n_wt = nrow(wt),
    n_mut = nrow(mut),
    observed_delta_median = median(mut_vals, na.rm = TRUE) - median(wt_vals, na.rm = TRUE),
    observed_delta_mean = mean(mut_vals, na.rm = TRUE) - mean(wt_vals, na.rm = TRUE),
    observed_late_thr = late_thr,
    observed_delta_latefrac = mean(mut_vals >= late_thr, na.rm = TRUE) - mean(wt_vals >= late_thr, na.rm = TRUE),
    observed_wilcox_p = wilcox_p
  )
}

sanitize_id <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

analyze_one_group <- function(df, analysis_id, display_name, B, late_q, n_start, n_step, out_tbl, out_plt) {
  wt_n <- sum(df$genotype2 == "WT")
  mut_n <- sum(df$genotype2 == "MUT")
  n_min <- min(wt_n, mut_n)
  if (!is.finite(n_min) || n_min < 50L) {
    msg("  skip ", analysis_id, " (n_min too small: ", n_min, ")")
    return(NULL)
  }

  local_step <- if (n_min < n_start) {
    if (n_min >= 200L) 100L else 50L
  } else {
    n_step
  }
  local_start <- if (n_min < n_start) local_step else n_start

  n_values <- seq.int(from = local_start, to = n_min, by = local_step)
  if (tail(n_values, 1) != n_min) n_values <- unique(c(n_values, n_min))

  observed_stats <- compute_observed_stats(df, late_q = late_q)
  resamp_df <- run_variable_n_pt(df, n_values = n_values, n_iter = B, late_q = late_q, seed = 1L)
  fwrite(resamp_df, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_resamples.csv")))
  bin_df <- run_variable_n_bins(df, n_values = n_values, n_iter = B, n_bins = 10L, seed = 1L)
  fwrite(bin_df, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_decile_resamples.csv")))

  median_tbl <- summarize_metric(resamp_df, "delta_median")
  median_tbl[, `:=`(analysis_id = analysis_id, display_name = display_name)]
  fwrite(median_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_delta_median_summary.csv")))

  late_tbl <- summarize_metric(resamp_df, "delta_latefrac")
  late_tbl[, `:=`(
    analysis_id = analysis_id,
    display_name = display_name,
    mean_value_pct = 100 * mean_value,
    ci_low_pct = 100 * ci_low,
    ci_high_pct = 100 * ci_high
  )]
  fwrite(late_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_delta_latefrac_summary.csv")))

  p_med <- plot_metric_vs_n(
    summary_tbl = median_tbl,
    x_col = "mean_value",
    lo_col = "ci_low",
    hi_col = "ci_high",
    title = paste0(display_name, ": variable-N pseudotime shift"),
    x_lab = "Delta median pseudotime (MUT - WT)"
  )
  save_pdf(p_med, file.path(out_plt, paste0("09j_", analysis_id, "_variableN_delta_median")), w = 8.2, h = 5.2)

  p_late <- plot_metric_vs_n(
    summary_tbl = late_tbl,
    x_col = "mean_value_pct",
    lo_col = "ci_low_pct",
    hi_col = "ci_high_pct",
    title = paste0(display_name, ": variable-N late-state shift"),
    x_lab = "Delta late fraction (percentage points, MUT - WT)"
  )
  save_pdf(p_late, file.path(out_plt, paste0("09j_", analysis_id, "_variableN_delta_latefrac")), w = 8.2, h = 5.2)

  if (nrow(bin_df) > 0) {
    bin_tbl <- summarize_bins(bin_df)
    fwrite(bin_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_decile_summary.csv")))
    bin_geno_tbl <- summarize_bins_by_genotype(bin_df)
    fwrite(bin_geno_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_decile_by_genotype_summary.csv")))

    p_bin <- plot_bins_vs_n(
      summary_tbl = bin_tbl,
      title = paste0(display_name, ": variable-N decile occupancy shift")
    )
    save_pdf(p_bin, file.path(out_plt, paste0("09j_", analysis_id, "_variableN_decile_occupancy")), w = 11.5, h = 7.5)

    pooled_bin_tbl <- bin_df[, .(
      mean_value_pct = 100 * mean(delta_frac, na.rm = TRUE),
      ci_low_pct = 100 * quantile(delta_frac, 0.025, na.rm = TRUE),
      ci_high_pct = 100 * quantile(delta_frac, 0.975, na.rm = TRUE)
    ), by = pt_bin][order(pt_bin)]
    fwrite(pooled_bin_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_decile_pooled_summary.csv")))

    p_bin_pool <- plot_bins_pooled(
      tbl = pooled_bin_tbl,
      title = paste0(display_name, ": pooled decile occupancy shift"),
      subtitle = paste0(B, " iterations per N across all tested N values")
    )
    save_pdf(p_bin_pool, file.path(out_plt, paste0("09j_", analysis_id, "_variableN_decile_pooled_lollipop")), w = 9.5, h = 6.5)

    pooled_occ_tbl <- bin_df[!is.na(genotype2), .(
      frac_mean = mean(frac, na.rm = TRUE),
      frac_ci_low = quantile(frac, 0.025, na.rm = TRUE),
      frac_ci_high = quantile(frac, 0.975, na.rm = TRUE)
    ), by = .(genotype2, pt_bin)][order(genotype2, pt_bin)]
    fwrite(pooled_occ_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_occupancy_by_bin_ci_summary.csv")))

    p_occ_pool <- plot_occupancy_by_bin_ci(
      summary_tbl = pooled_occ_tbl,
      title = paste0(display_name, ": pseudotime-bin occupancy by genotype"),
      subtitle = paste0("Variable-N balanced occupancy (", B, " iterations per N, mean ±95% CI)")
    )
    save_pdf(p_occ_pool, file.path(out_plt, paste0("09j_", analysis_id, "_variableN_occupancy_by_bin_ci")), w = 8.2, h = 4.9)

    pooled_occ_delta_tbl <- bin_df[, .(
      mean_value = mean(delta_frac, na.rm = TRUE),
      ci_low = quantile(delta_frac, 0.025, na.rm = TRUE),
      ci_high = quantile(delta_frac, 0.975, na.rm = TRUE)
    ), by = pt_bin][order(pt_bin)]
    fwrite(pooled_occ_delta_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_variableN_occupancy_delta_by_bin_ci_summary.csv")))

    p_occ_delta_pool <- plot_occupancy_delta_by_bin_ci(
      summary_tbl = pooled_occ_delta_tbl,
      title = paste0(display_name, ": pseudotime-bin occupancy delta"),
      subtitle = paste0("Variable-N balanced occupancy (", B, " iterations per N, mean ±95% CI)")
    )
    save_pdf(p_occ_delta_pool, file.path(out_plt, paste0("09j_", analysis_id, "_variableN_occupancy_delta_by_bin_ci")), w = 8.2, h = 4.9)

    endpoint_occ_tbl <- bin_df[N == max(N) & !is.na(genotype2), .(
      frac_mean = mean(frac, na.rm = TRUE),
      frac_ci_low = quantile(frac, 0.025, na.rm = TRUE),
      frac_ci_high = quantile(frac, 0.975, na.rm = TRUE)
    ), by = .(genotype2, pt_bin)][order(genotype2, pt_bin)]
    fwrite(endpoint_occ_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_endpoint_occupancy_by_bin_ci_summary.csv")))

    p_occ_endpoint <- plot_occupancy_by_bin_ci(
      summary_tbl = endpoint_occ_tbl,
      title = paste0(display_name, ": pseudotime-bin occupancy by genotype"),
      subtitle = paste0("Maximum balanced N only (", B, " iterations, mean ±95% CI)")
    )
    save_pdf(p_occ_endpoint, file.path(out_plt, paste0("09j_", analysis_id, "_endpoint_occupancy_by_bin_ci")), w = 8.2, h = 4.9)

    endpoint_occ_delta_tbl <- bin_df[N == max(N), .(
      mean_value = mean(delta_frac, na.rm = TRUE),
      ci_low = quantile(delta_frac, 0.025, na.rm = TRUE),
      ci_high = quantile(delta_frac, 0.975, na.rm = TRUE)
    ), by = pt_bin][order(pt_bin)]
    fwrite(endpoint_occ_delta_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_endpoint_occupancy_delta_by_bin_ci_summary.csv")))

    p_occ_delta_endpoint <- plot_occupancy_delta_by_bin_ci(
      summary_tbl = endpoint_occ_delta_tbl,
      title = paste0(display_name, ": pseudotime-bin occupancy delta"),
      subtitle = paste0("Maximum balanced N only (", B, " iterations, mean ±95% CI)")
    )
    save_pdf(p_occ_delta_endpoint, file.path(out_plt, paste0("09j_", analysis_id, "_endpoint_occupancy_delta_by_bin_ci")), w = 8.2, h = 4.9)
  }

  pooled_med <- data.table(
    analysis_id = analysis_id,
    display_name = display_name,
    mean_value = mean(resamp_df$delta_median, na.rm = TRUE),
    ci_low = quantile(resamp_df$delta_median, 0.025, na.rm = TRUE),
    ci_high = quantile(resamp_df$delta_median, 0.975, na.rm = TRUE),
    n_total_resamples = nrow(resamp_df)
  )

  pooled_late <- data.table(
    analysis_id = analysis_id,
    display_name = display_name,
    mean_value_pct = 100 * mean(resamp_df$delta_latefrac, na.rm = TRUE),
    ci_low_pct = 100 * quantile(resamp_df$delta_latefrac, 0.025, na.rm = TRUE),
    ci_high_pct = 100 * quantile(resamp_df$delta_latefrac, 0.975, na.rm = TRUE),
    n_total_resamples = nrow(resamp_df)
  )

  endpoint_med <- copy(median_tbl[N == max(N)][1])
  endpoint_late <- copy(late_tbl[N == max(N)][1])

  stats_tbl <- cbind(
    data.table(
      analysis_id = analysis_id,
      display_name = display_name
    ),
    observed_stats,
    data.table(
      pooled_delta_median_mean = pooled_med$mean_value,
      pooled_delta_median_ci_low = pooled_med$ci_low,
      pooled_delta_median_ci_high = pooled_med$ci_high,
      pooled_delta_median_ci_excludes_zero = (pooled_med$ci_low > 0) || (pooled_med$ci_high < 0),
      pooled_delta_median_same_sign_frac = mean(sign(resamp_df$delta_median) == sign(mean(resamp_df$delta_median, na.rm = TRUE)), na.rm = TRUE),
      pooled_delta_latefrac_mean_pct = pooled_late$mean_value_pct,
      pooled_delta_latefrac_ci_low_pct = pooled_late$ci_low_pct,
      pooled_delta_latefrac_ci_high_pct = pooled_late$ci_high_pct,
      pooled_delta_latefrac_ci_excludes_zero = (pooled_late$ci_low_pct > 0) || (pooled_late$ci_high_pct < 0),
      pooled_delta_latefrac_same_sign_frac = mean(sign(resamp_df$delta_latefrac) == sign(mean(resamp_df$delta_latefrac, na.rm = TRUE)), na.rm = TRUE),
      endpoint_N = max(median_tbl$N, na.rm = TRUE),
      endpoint_delta_median_mean = endpoint_med$mean_value,
      endpoint_delta_median_ci_low = endpoint_med$ci_low,
      endpoint_delta_median_ci_high = endpoint_med$ci_high,
      endpoint_delta_median_ci_excludes_zero = (endpoint_med$ci_low > 0) || (endpoint_med$ci_high < 0),
      endpoint_delta_latefrac_mean_pct = endpoint_late$mean_value_pct,
      endpoint_delta_latefrac_ci_low_pct = endpoint_late$ci_low_pct,
      endpoint_delta_latefrac_ci_high_pct = endpoint_late$ci_high_pct,
      endpoint_delta_latefrac_ci_excludes_zero = (endpoint_late$ci_low_pct > 0) || (endpoint_late$ci_high_pct < 0)
    )
  )
  fwrite(stats_tbl, file.path(out_tbl, paste0("09j_", analysis_id, "_stats_summary.csv")))

  run_summary <- data.table(
    analysis_id = analysis_id,
    display_name = display_name,
    n_wt = wt_n,
    n_mut = mut_n,
    n_min = n_min,
    n_start_used = local_start,
    n_step_used = local_step,
    n_values_tested = length(n_values),
    iterations_per_N = B
  )

  list(
    run_summary = run_summary,
    stats_tbl = stats_tbl,
    pooled_med = pooled_med,
    pooled_late = pooled_late,
    endpoint_med = endpoint_med,
    endpoint_late = endpoint_late
  )
}

out_base <- file.path(CFG$paths$out_dir, "results_09j")
out_tbl <- file.path(out_base, "tables")
out_plt <- file.path(out_base, "plots")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plt, recursive = TRUE, showWarnings = FALSE)

n_start <- 200L
n_step <- 200L
B <- 50L
late_q <- CFG$analysis$late_quantile %||% 0.80

summary_rows <- list()
sr_i <- 0L
stats_rows <- list()
st_i <- 0L
pooled_median_rows <- list()
pm_i <- 0L
pooled_late_rows <- list()
pl_i <- 0L
endpoint_median_rows <- list()
em_i <- 0L
endpoint_late_rows <- list()
el_i <- 0L

for (lineage_name in names(CFG$lineages)) {
  msg("09j: ", lineage_name)
  in_obj <- file.path(CFG$paths$out_dir, "mapping", "objects", paste0("09b_", lineage_name, "_mapped_query.rds"))
  if (!file.exists(in_obj)) {
    msg("  skip (missing mapped object): ", in_obj)
    next
  }

  seu <- readRDS(in_obj)
  md0 <- seu@meta.data
  if ("cell_id" %in% colnames(md0)) md0$cell_id <- NULL
  md <- as.data.table(md0, keep.rownames = "cell_id")
  gcol <- pick_genotype_col(md, CFG$cols$genotype_col)
  if (is.na(gcol) || !"pseudotime_transfer" %in% colnames(md)) {
    msg("  skip (missing genotype or pseudotime_transfer)")
    next
  }

  md[, genotype2 := canon_genotype(get(gcol), CFG$genotype$wt_values, CFG$genotype$mut_values)]
  df <- md[!is.na(genotype2) & is.finite(pseudotime_transfer)]
  res <- analyze_one_group(
    df = df,
    analysis_id = lineage_name,
    display_name = lineage_name,
    B = B,
    late_q = late_q,
    n_start = n_start,
    n_step = n_step,
    out_tbl = out_tbl,
    out_plt = out_plt
  )
  if (is.null(res)) next

  sr_i <- sr_i + 1L
  summary_rows[[sr_i]] <- copy(res$run_summary)[, lineage_name := analysis_id]
  st_i <- st_i + 1L
  stats_rows[[st_i]] <- copy(res$stats_tbl)[, lineage_name := analysis_id]
  pm_i <- pm_i + 1L
  pooled_median_rows[[pm_i]] <- copy(res$pooled_med)[, lineage_name := analysis_id]
  pl_i <- pl_i + 1L
  pooled_late_rows[[pl_i]] <- copy(res$pooled_late)[, lineage_name := analysis_id]
  em_i <- em_i + 1L
  endpoint_median_rows[[em_i]] <- copy(res$endpoint_med)[, lineage_name := analysis_id]
  el_i <- el_i + 1L
  endpoint_late_rows[[el_i]] <- copy(res$endpoint_late)[, lineage_name := analysis_id]
}

subgroup_specs <- list(
  list(
    parent_lineage = "dorsal_rg_ipc_exc",
    label_col = "celltype_broad",
    subgroups = list(
      dorsal_upper_layer_EN = c("Upper layer EN", "Upper layer En"),
      dorsal_deep_layer_EN = c("Deep layer EN", "Deep layer En"),
      dorsal_immature_astrocytes = c("Immature astrocytes", "Immature Astrocytes")
    )
  ),
  list(
    parent_lineage = "ventral_lge_spn",
    label_col = "celltype_broad",
    subgroups = list(
      ventral_SPNs = c("SPNs"),
      ventral_LGE_derived_IN = c("LGE-IN prec", "LGE-IN precursors")
    )
  ),
  list(
    parent_lineage = "ventral_lge_spn",
    label_col = "celltype_label",
    subgroups = list(
      ventral_SPN_D1 = c("SPN-D1"),
      ventral_SPN_D2 = c("SPN-D2")
    )
  )
)

for (spec in subgroup_specs) {
  parent_lineage <- spec$parent_lineage
  in_obj <- file.path(CFG$paths$out_dir, "mapping", "objects", paste0("09b_", parent_lineage, "_mapped_query.rds"))
  if (!file.exists(in_obj)) {
    msg("09j subgroup: skip missing parent object ", in_obj)
    next
  }
  msg("09j subgroup parent: ", parent_lineage)
  seu <- readRDS(in_obj)
  md0 <- seu@meta.data
  if ("cell_id" %in% colnames(md0)) md0$cell_id <- NULL
  md <- as.data.table(md0, keep.rownames = "cell_id")
  gcol <- pick_genotype_col(md, CFG$cols$genotype_col)
  if (is.na(gcol) || !"pseudotime_transfer" %in% colnames(md) || !spec$label_col %in% colnames(md)) {
    msg("  skip subgroup set for ", parent_lineage, " (missing required columns)")
    next
  }
  md[, genotype2 := canon_genotype(get(gcol), CFG$genotype$wt_values, CFG$genotype$mut_values)]
  md <- md[!is.na(genotype2) & is.finite(pseudotime_transfer)]

  for (sub_name in names(spec$subgroups)) {
    keep_labels <- spec$subgroups[[sub_name]]
    sdf <- md[get(spec$label_col) %in% keep_labels]
    analysis_id <- sanitize_id(sub_name)
    display_name <- gsub("_", " ", sub_name)
    msg("  subgroup: ", display_name)
    analyze_one_group(
      df = sdf,
      analysis_id = analysis_id,
      display_name = display_name,
      B = B,
      late_q = late_q,
      n_start = n_start,
      n_step = n_step,
      out_tbl = out_tbl,
      out_plt = out_plt
    )
  }
}

if (length(summary_rows) > 0) {
  fwrite(rbindlist(summary_rows, fill = TRUE), file.path(out_tbl, "09j_variableN_run_summary.csv"))
}

if (length(stats_rows) > 0) {
  fwrite(rbindlist(stats_rows, fill = TRUE), file.path(out_tbl, "09j_variableN_stats_summary.csv"))
}

if (length(pooled_median_rows) > 0) {
  pooled_med_tbl <- rbindlist(pooled_median_rows, fill = TRUE)
  fwrite(pooled_med_tbl, file.path(out_tbl, "09j_variableN_pooled_delta_median_summary.csv"))
  p <- plot_lollipop(
    tbl = pooled_med_tbl,
    value_col = "mean_value",
    lo_col = "ci_low",
    hi_col = "ci_high",
    title = "Pseudotime shift pooled across all tested N values",
    subtitle = paste0(B, " iterations per N"),
    x_lab = "Average delta median pseudotime (MUT - WT)"
  )
  save_pdf(p, file.path(out_plt, "09j_variableN_pooled_delta_median_lollipop"), w = 9.0, h = 6.5)
}

if (length(pooled_late_rows) > 0) {
  pooled_late_tbl <- rbindlist(pooled_late_rows, fill = TRUE)
  fwrite(pooled_late_tbl, file.path(out_tbl, "09j_variableN_pooled_delta_latefrac_summary.csv"))
  p <- plot_lollipop(
    tbl = pooled_late_tbl,
    value_col = "mean_value_pct",
    lo_col = "ci_low_pct",
    hi_col = "ci_high_pct",
    title = "Late-state shift pooled across all tested N values",
    subtitle = paste0(B, " iterations per N"),
    x_lab = "Average delta late fraction (percentage points, MUT - WT)"
  )
  save_pdf(p, file.path(out_plt, "09j_variableN_pooled_delta_latefrac_lollipop"), w = 9.0, h = 6.5)
}

if (length(endpoint_median_rows) > 0) {
  endpoint_med_tbl <- rbindlist(endpoint_median_rows, fill = TRUE)
  fwrite(endpoint_med_tbl, file.path(out_tbl, "09j_variableN_endpoint_delta_median_summary.csv"))
  p <- plot_lollipop(
    tbl = endpoint_med_tbl,
    value_col = "mean_value",
    lo_col = "ci_low",
    hi_col = "ci_high",
    title = "Pseudotime shift at maximum balanced N",
    subtitle = paste0(B, " iterations at N = n_min"),
    x_lab = "Delta median pseudotime (MUT - WT)"
  )
  save_pdf(p, file.path(out_plt, "09j_variableN_endpoint_delta_median_lollipop"), w = 9.0, h = 6.5)
}

if (length(endpoint_late_rows) > 0) {
  endpoint_late_tbl <- rbindlist(endpoint_late_rows, fill = TRUE)
  fwrite(endpoint_late_tbl, file.path(out_tbl, "09j_variableN_endpoint_delta_latefrac_summary.csv"))
  p <- plot_lollipop(
    tbl = endpoint_late_tbl,
    value_col = "mean_value_pct",
    lo_col = "ci_low_pct",
    hi_col = "ci_high_pct",
    title = "Late-state shift at maximum balanced N",
    subtitle = paste0(B, " iterations at N = n_min"),
    x_lab = "Delta late fraction (percentage points, MUT - WT)"
  )
  save_pdf(p, file.path(out_plt, "09j_variableN_endpoint_delta_latefrac_lollipop"), w = 9.0, h = 6.5)
}

msg("09j complete. Outputs under: ", out_base)
