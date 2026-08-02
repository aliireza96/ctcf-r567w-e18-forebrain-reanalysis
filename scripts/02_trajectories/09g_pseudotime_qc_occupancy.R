#!/usr/bin/env Rscript

# ============================================================
# 09g_pseudotime_qc_occupancy.R
# Additional pseudotime diagnostics:
# - Pseudotime-bin occupancy by genotype with downsampling CI
# - Mapping confidence QC by genotype/lineage
# - Mapping coverage summary versus query target cells
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

assign_bins <- function(x, n_bins = 10L) {
  qs <- unique(as.numeric(quantile(x, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE)))
  if (length(qs) < 3) return(rep("D01", length(x)))
  b <- cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE)
  b[is.na(b)] <- 1L
  paste0("D", sprintf("%02d", b))
}

downsample_bin_occupancy <- function(df, n_iter = 300L, seed = 1L, n_bins = 10L) {
  set.seed(seed)
  d <- copy(df)
  wt <- d[genotype2 == "WT"]
  mut <- d[genotype2 == "MUT"]
  n <- min(nrow(wt), nrow(mut))
  if (n < 10) return(data.table())

  out <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    wi <- wt[sample(.N, n)]
    mi <- mut[sample(.N, n)]
    r <- rbindlist(list(wi, mi), fill = TRUE)
    r[, pt_bin := assign_bins(pseudotime_transfer, n_bins = n_bins)]
    tab <- r[, .(n_cells = .N), by = .(genotype2, pt_bin)]
    tab[, frac := n_cells / sum(n_cells), by = genotype2]
    tab[, iter := i]
    out[[i]] <- tab
  }
  rbindlist(out, fill = TRUE)
}

out_base <- file.path(CFG$paths$out_dir, "results_09g")
out_tbl <- file.path(out_base, "tables")
out_plt <- file.path(out_base, "plots")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plt, recursive = TRUE, showWarnings = FALSE)

query <- readRDS(CFG$paths$query_e18_rds)
qmeta <- as.data.table(query@meta.data)
broad_col <- CFG$cols$broad_col
if (!broad_col %in% colnames(qmeta)) broad_col <- "celltype_broad"

lineage_summary <- list()
ls_i <- 0L

for (lineage_name in names(CFG$lineages)) {
  msg("09g: ", lineage_name)
  obj_path <- file.path(CFG$paths$out_dir, "mapping", "objects", paste0("09b_", lineage_name, "_mapped_query.rds"))
  if (!file.exists(obj_path)) {
    msg("  skip (missing): ", obj_path)
    next
  }

  seu <- readRDS(obj_path)
  md0 <- seu@meta.data
  if ("cell_id" %in% colnames(md0)) md0$cell_id <- NULL
  md <- as.data.table(md0, keep.rownames = "cell_id")

  gcol <- pick_genotype_col(md, CFG$cols$genotype_col)
  if (is.na(gcol)) {
    msg("  skip (genotype column not found)")
    next
  }

  req <- c("pseudotime_transfer", "branch_transfer_score", "ref_label_transfer_score")
  for (x in req) if (!x %in% colnames(md)) md[, (x) := NA_real_]

  md[, genotype2 := canon_genotype(get(gcol), CFG$genotype$wt_values, CFG$genotype$mut_values)]
  df <- md[!is.na(genotype2) & is.finite(pseudotime_transfer)]
  if (nrow(df) == 0) {
    msg("  skip (no usable cells)")
    next
  }

  # Coverage summary
  keep_broad <- CFG$lineages[[lineage_name]]$query_keep_broad
  n_target <- if (!is.null(keep_broad) && broad_col %in% colnames(qmeta)) {
    sum(qmeta[[broad_col]] %in% keep_broad, na.rm = TRUE)
  } else {
    NA_integer_
  }

  ls_i <- ls_i + 1L
  lineage_summary[[ls_i]] <- data.table(
    lineage_name = lineage_name,
    genotype_source_col = gcol,
    n_target_query_cells = n_target,
    n_mapped_cells = nrow(df),
    mapped_fraction_of_target = if (is.na(n_target) || n_target == 0) NA_real_ else nrow(df) / n_target,
    n_wt = sum(df$genotype2 == "WT"),
    n_mut = sum(df$genotype2 == "MUT")
  )

  # Occupancy bins (observed)
  df[, pt_bin := assign_bins(pseudotime_transfer, n_bins = 10L)]
  occ <- df[, .(n_cells = .N), by = .(genotype2, pt_bin)]
  occ[, frac := n_cells / sum(n_cells), by = genotype2]
  occ[, lineage_name := lineage_name]
  fwrite(occ, file.path(out_tbl, paste0("09g_", lineage_name, "_occupancy_observed.csv")))

  # Occupancy bins (downsample CI)
  ds <- downsample_bin_occupancy(
    df,
    n_iter = as.integer(CFG$analysis$downsample_iterations %||% 300L),
    seed = 1L,
    n_bins = 10L
  )
  if (nrow(ds) > 0) {
    fwrite(ds, file.path(out_tbl, paste0("09g_", lineage_name, "_occupancy_downsample.csv")))
    ds_sum <- ds[, .(
      frac_mean = mean(frac, na.rm = TRUE),
      frac_ci_low = quantile(frac, 0.025, na.rm = TRUE),
      frac_ci_high = quantile(frac, 0.975, na.rm = TRUE)
    ), by = .(genotype2, pt_bin)]
    ds_sum[, lineage_name := lineage_name]
    fwrite(ds_sum, file.path(out_tbl, paste0("09g_", lineage_name, "_occupancy_downsample_summary.csv")))

    ds_sum[, pt_bin_num := as.integer(sub("^D", "", pt_bin))]
    p_occ <- ggplot(ds_sum, aes(pt_bin_num, frac_mean, color = genotype2, fill = genotype2, group = genotype2)) +
      geom_ribbon(aes(ymin = frac_ci_low, ymax = frac_ci_high), alpha = 0.18, linewidth = 0) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 1.5) +
      scale_x_continuous(breaks = seq_len(max(ds_sum$pt_bin_num, na.rm = TRUE))) +
      theme_classic(base_size = 12) +
      labs(
        title = paste0(lineage_name, ": pseudotime-bin occupancy by genotype"),
        subtitle = "Downsampled equal-cell occupancy (mean ±95% CI)",
        x = "Pseudotime decile (low -> high)",
        y = "Cell fraction"
      )
    save_pdf(p_occ, file.path(out_plt, paste0("09g_", lineage_name, "_occupancy_by_bin_ci")), w = 8.2, h = 4.9)
  }

  # Confidence QC by genotype
  conf_long <- melt(
    df[, .(genotype2, branch_transfer_score, ref_label_transfer_score, pseudotime_transfer_score)],
    id.vars = "genotype2",
    variable.name = "metric",
    value.name = "score"
  )
  conf_long <- conf_long[!is.na(score)]
  if (nrow(conf_long) > 0) {
    conf_sum <- conf_long[, .(
      n = .N,
      median = median(score, na.rm = TRUE),
      q1 = quantile(score, 0.25, na.rm = TRUE),
      q3 = quantile(score, 0.75, na.rm = TRUE),
      frac_lt_0.8 = mean(score < 0.8, na.rm = TRUE),
      frac_lt_0.9 = mean(score < 0.9, na.rm = TRUE)
    ), by = .(metric, genotype2)]
    conf_sum[, lineage_name := lineage_name]
    fwrite(conf_sum, file.path(out_tbl, paste0("09g_", lineage_name, "_confidence_summary.csv")))

    p_conf <- ggplot(conf_long, aes(genotype2, score, fill = genotype2)) +
      geom_violin(trim = FALSE, alpha = 0.75) +
      geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
      facet_wrap(~ metric, scales = "free_y") +
      theme_classic(base_size = 12) +
      labs(
        title = paste0(lineage_name, ": mapping confidence by genotype"),
        x = NULL,
        y = "Score"
      )
    save_pdf(p_conf, file.path(out_plt, paste0("09g_", lineage_name, "_confidence_violin")), w = 10.5, h = 4.8)
  }
}

if (length(lineage_summary) > 0) {
  cov_tbl <- rbindlist(lineage_summary, fill = TRUE)
  fwrite(cov_tbl, file.path(out_tbl, "09g_mapping_coverage_summary.csv"))

  p_cov <- ggplot(cov_tbl, aes(lineage_name, mapped_fraction_of_target)) +
    geom_col(fill = "#4C78A8", width = 0.7) +
    geom_text(aes(label = paste0("n=", n_mapped_cells, "/", n_target_query_cells)), vjust = -0.35, size = 3.4) +
    ylim(0, max(1, max(cov_tbl$mapped_fraction_of_target, na.rm = TRUE) + 0.1)) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(
      title = "Mapping coverage by lineage",
      subtitle = "Mapped cells relative to query target pool",
      x = NULL,
      y = "Mapped fraction"
    )
  save_pdf(p_cov, file.path(out_plt, "09g_mapping_coverage_summary"), w = 8.4, h = 4.8)
}

msg("09g complete. Outputs under: ", out_base)
