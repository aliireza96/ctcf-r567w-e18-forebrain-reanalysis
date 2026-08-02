# ============================================================
# 04b_composition_variableN_sensitivity.R
# Variable-N equal-cell downsampling sensitivity analysis for
# cell-type composition differences (MUT - WT).
# ============================================================

source("scripts/00_setup.R")

use_kept_only <- TRUE
use_cortical_striatal_scope <- TRUE
analysis_scope_label <- "cortical/striatal analysis scope"

# Any metadata annotation columns to analyze
celltype_cols <- c("celltype_broad", "celltype_label")

# Variable-N settings
n_start <- 500L
n_step <- 500L
# Use repeated balanced subsampling at each fixed N.
# 50 iterations is usually enough for stable CIs without making the plot pipeline heavy.
B <- 50L
set.seed(1)

in_obj_full <- file.path(OUT$objects, "03_seu_annotated.rds")
in_obj_kept <- file.path(OUT$objects, "03_seu_annotated_KEPT.rds")
in_obj_cortstr <- file.path(OUT$objects, "03_seu_annotated_CORTSTR.rds")
in_obj_tel  <- file.path(OUT$objects, "03_seu_annotated_TEL.rds")
in_obj <- if (use_cortical_striatal_scope && file.exists(in_obj_cortstr)) {
  in_obj_cortstr
} else if (use_cortical_striatal_scope && file.exists(in_obj_tel)) {
  in_obj_tel
} else if (use_kept_only) {
  in_obj_kept
} else {
  in_obj_full
}
stopifnot(file.exists(in_obj))

seu <- readRDS(in_obj)
message("Loaded composition object for variable-N analysis: ", in_obj)

run_variable_n_composition <- function(celltype_col) {
  stopifnot(celltype_col %in% colnames(seu@meta.data))

  meta <- seu@meta.data %>%
    mutate(
      celltype = seu[[celltype_col, drop = TRUE]],
      condition = factor(condition, levels = c("wt", "mut"))
    ) %>%
    filter(!is.na(celltype))

  counts <- meta %>%
    count(condition, celltype, name = "n")

  totals <- meta %>%
    count(condition, name = "total") %>%
    tidyr::pivot_wider(names_from = condition, values_from = total)

  counts_wide <- counts %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0)

  mut_total <- totals$mut
  wt_total  <- totals$wt
  eps <- 0.5 / min(mut_total, wt_total)

  calc_p <- function(m, w) {
    mat <- matrix(c(m, mut_total - m, w, wt_total - w), nrow = 2, byrow = TRUE)
    fisher.test(mat)$p.value
  }

  stats_tbl <- counts_wide %>%
    mutate(
      mut_prop = mut / mut_total,
      wt_prop  = wt / wt_total,
      log2fc   = log2((mut_prop + eps) / (wt_prop + eps)),
      p_val    = mapply(calc_p, mut, wt),
      p_adj    = p.adjust(p_val, method = "BH"),
      signif = dplyr::case_when(
        p_adj < 0.01 ~ "FDR < 0.01",
        p_adj < 0.05 ~ "FDR 0.01-0.05",
        TRUE ~ "ns"
      ),
      sig_label = dplyr::case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01 ~ "**",
        p_adj < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  counts_by_cond <- table(meta$condition)
  n_min <- min(counts_by_cond)
  n_values <- seq.int(from = n_start, to = n_min, by = n_step)
  if (tail(n_values, 1) != n_min) {
    n_values <- unique(c(n_values, n_min))
  }
  n_values <- n_values[n_values > 0]

  if (length(n_values) == 0) {
    stop("No valid N values generated for ", celltype_col, ".")
  }

  suffix <- celltype_col

  resamp <- lapply(n_values, function(n_equal) {
    iter_list <- lapply(seq_len(B), function(b) {
      sub <- meta %>%
        group_by(condition) %>%
        slice_sample(n = n_equal, replace = FALSE) %>%
        ungroup()

      sub %>%
        count(condition, celltype, name = "n") %>%
        group_by(condition) %>%
        mutate(prop = n / sum(n)) %>%
        ungroup() %>%
        dplyr::select(condition, celltype, prop) %>%
        tidyr::pivot_wider(names_from = condition, values_from = prop, values_fill = 0) %>%
        mutate(
          delta = mut - wt,
          N = n_equal,
          iter = b
        )
    })
    bind_rows(iter_list)
  })

  resamp_df <- bind_rows(resamp)

  write_csv(
    resamp_df,
    file.path(OUT$tables, paste0("04b_variableN_resamples_", suffix, ".csv"))
  )

  summary_tbl <- resamp_df %>%
    group_by(celltype, N) %>%
    summarise(
      mean_delta = mean(delta, na.rm = TRUE),
      ci_low = quantile(delta, 0.025, na.rm = TRUE),
      ci_high = quantile(delta, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean_delta_pct = 100 * mean_delta,
      ci_low_pct = 100 * ci_low,
      ci_high_pct = 100 * ci_high
    ) %>%
    arrange(celltype, N)

  write_csv(
    summary_tbl,
    file.path(OUT$tables, paste0("04b_variableN_delta_summary_", suffix, ".csv"))
  )

  p_faceted <- summary_tbl %>%
    ggplot(aes(y = N, x = mean_delta_pct, group = celltype)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = ci_low_pct, xmax = ci_high_pct), height = 120, color = "#9ecae1", linewidth = 1.0) +
    geom_path(color = "#2171b5", linewidth = 0.8) +
    geom_point(color = "#08519c", size = 1.8) +
    facet_wrap(~ celltype, scales = "free_x") +
    theme_classic(12) +
    labs(
      title = paste0("Variable-N composition sensitivity (", celltype_col, ")"),
      subtitle = paste0(analysis_scope_label, " | Mean delta in percentage points (mut - wt), N = ", min(n_values), " to ", max(n_values), ", ", B, " iterations per N"),
      x = "Delta fraction (percentage points)",
      y = "Equal cells sampled per condition (N)"
    )
  save_plot(p_faceted, paste0("composition/Composition_variableN_delta_faceted_", suffix), w = 12, h = 10)

  p_overlay <- summary_tbl %>%
    mutate(celltype = forcats::fct_inorder(celltype)) %>%
    ggplot(aes(x = N, y = mean_delta_pct, color = celltype, fill = celltype, group = celltype)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_ribbon(aes(ymin = ci_low_pct, ymax = ci_high_pct), alpha = 0.10, linewidth = 0) +
    geom_line(linewidth = 0.9) +
    theme_classic(12) +
    labs(
      title = paste0("Variable-N composition sensitivity overlay (", celltype_col, ")"),
      subtitle = analysis_scope_label,
      x = "Equal cells sampled per condition (N)",
      y = "Delta fraction (percentage points)",
      color = NULL,
      fill = NULL
    )
  save_plot(p_overlay, paste0("composition/Composition_variableN_delta_overlay_", suffix), w = 11, h = 7)

  # Pooled summary across all tested N values and all iterations.
  pooled_tbl <- resamp_df %>%
    group_by(celltype) %>%
    summarise(
      mean_delta = mean(delta, na.rm = TRUE),
      ci_low = quantile(delta, 0.025, na.rm = TRUE),
      ci_high = quantile(delta, 0.975, na.rm = TRUE),
      n_total_resamples = dplyr::n(),
      .groups = "drop"
    ) %>%
    mutate(
      mean_delta_pct = 100 * mean_delta,
      ci_low_pct = 100 * ci_low,
      ci_high_pct = 100 * ci_high,
      direction = ifelse(mean_delta_pct >= 0, "Enriched in mut", "Enriched in wt")
    ) %>%
    left_join(
      stats_tbl %>% dplyr::select(celltype, mut, wt, mut_prop, wt_prop, log2fc, p_val, p_adj, signif, sig_label),
      by = "celltype"
    ) %>%
    arrange(mean_delta_pct) %>%
    mutate(celltype = factor(celltype, levels = celltype))

  write_csv(
    pooled_tbl,
    file.path(OUT$tables, paste0("04b_variableN_pooled_summary_", suffix, ".csv"))
  )

  p_pooled <- ggplot(pooled_tbl, aes(x = mean_delta_pct, y = celltype, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = 0, xend = mean_delta_pct, y = celltype, yend = celltype), linewidth = 1.0) +
    geom_errorbarh(aes(xmin = ci_low_pct, xmax = ci_high_pct), height = 0.18, linewidth = 0.8) +
    geom_point(size = 2.8) +
    geom_text(
      aes(
        x = ifelse(mean_delta_pct >= 0, ci_high_pct + 0.6, ci_low_pct - 0.6),
        y = celltype,
        label = sig_label
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3.7,
      fontface = "bold"
    ) +
    scale_color_manual(values = c("Enriched in mut" = "#d95f02", "Enriched in wt" = "#1f77b4")) +
    theme_classic(12) +
    labs(
      title = paste0("Composition shift pooled across all tested N values (", celltype_col, ")"),
      subtitle = paste0(analysis_scope_label, " | ", B, " iterations per N, pooled over N = ", min(n_values), " to ", max(n_values)),
      x = "Average delta fraction (percentage points, mut - wt)",
      y = NULL,
      color = NULL,
      caption = "Points/ribbons: variable-N balanced downsampling effect size. Asterisks: Fisher exact test on full counts with BH correction (* <0.05, ** <0.01, *** <0.001)."
    ) +
    coord_cartesian(
      xlim = c(min(pooled_tbl$ci_low_pct, na.rm = TRUE) - 1.0, max(pooled_tbl$ci_high_pct, na.rm = TRUE) + 1.3),
      clip = "off"
    ) +
    theme(plot.margin = margin(5.5, 24, 5.5, 5.5))
  save_plot(p_pooled, paste0("composition/Composition_variableN_pooled_lollipop_", suffix), w = 9, h = 7)

  # Endpoint summary at the maximum balanced N is easier to read in talks than the full ribbon grid.
  endpoint_tbl <- summary_tbl %>%
    filter(N == max(N)) %>%
    left_join(
      stats_tbl %>% dplyr::select(celltype, mut, wt, mut_prop, wt_prop, log2fc, p_val, p_adj, signif, sig_label),
      by = "celltype"
    ) %>%
    arrange(mean_delta_pct) %>%
    mutate(
      celltype = factor(celltype, levels = celltype),
      direction = ifelse(mean_delta_pct >= 0, "Enriched in mut", "Enriched in wt")
    )

  write_csv(
    endpoint_tbl,
    file.path(OUT$tables, paste0("04b_variableN_endpoint_summary_", suffix, ".csv"))
  )

  p_endpoint <- ggplot(endpoint_tbl, aes(x = mean_delta_pct, y = celltype, color = direction)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = 0, xend = mean_delta_pct, y = celltype, yend = celltype), linewidth = 1.0) +
    geom_errorbarh(aes(xmin = ci_low_pct, xmax = ci_high_pct), height = 0.18, linewidth = 0.8) +
    geom_point(size = 2.6) +
    geom_text(
      aes(
        x = ifelse(mean_delta_pct >= 0, ci_high_pct + 0.6, ci_low_pct - 0.6),
        y = celltype,
        label = sig_label
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3.7,
      fontface = "bold"
    ) +
    scale_color_manual(values = c("Enriched in mut" = "#d95f02", "Enriched in wt" = "#1f77b4")) +
    theme_classic(12) +
    labs(
      title = paste0("Composition shift at maximum balanced N (", celltype_col, ")"),
      subtitle = paste0(analysis_scope_label, " | N = ", max(n_values), ", ", B, " iterations"),
      x = "Delta fraction (percentage points, mut - wt)",
      y = NULL,
      color = NULL,
      caption = "Points/ribbons: variable-N balanced downsampling effect size. Asterisks: Fisher exact test on full counts with BH correction (* <0.05, ** <0.01, *** <0.001)."
    ) +
    coord_cartesian(
      xlim = c(min(endpoint_tbl$ci_low_pct, na.rm = TRUE) - 1.0, max(endpoint_tbl$ci_high_pct, na.rm = TRUE) + 1.3),
      clip = "off"
    ) +
    theme(plot.margin = margin(5.5, 24, 5.5, 5.5))
  save_plot(p_endpoint, paste0("composition/Composition_variableN_endpoint_lollipop_", suffix), w = 9, h = 7)
}

for (ct in celltype_cols) {
  run_variable_n_composition(ct)
}
