#!/usr/bin/env Rscript

# ============================================================
# 09c_analyze_genotype_pseudotime.R
# Analyze genotype effects on transferred pseudotime/branches
# and gene dynamics for each lineage.
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

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

save_pdf <- function(plot_obj, path_no_ext, w = 8.5, h = 6.0) {
  ggsave(paste0(path_no_ext, ".pdf"), plot_obj, width = w, height = h, useDingbats = FALSE)
}

fetch_expr_data <- function(seu, vars, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(seu)
  out <- tryCatch(
    FetchData(seu, vars = vars, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  out <- tryCatch(
    FetchData(seu, vars = vars, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  out <- tryCatch(
    FetchData(seu, vars = vars, assay = assay),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  stop("Unable to fetch expression data for assay=", assay, " and vars=", paste(vars, collapse = ","))
}

canon_genotype <- function(x, wt_values, mut_values) {
  x <- as.character(x)
  x_trim <- trimws(x)
  x_low <- tolower(x_trim)
  wt_low <- tolower(trimws(as.character(wt_values)))
  mut_low <- tolower(trimws(as.character(mut_values)))
  out <- rep(NA_character_, length(x))
  out[x_low %in% wt_low] <- "WT"
  out[x_low %in% mut_low] <- "MUT"
  # Regex fallback for common naming variants.
  out[is.na(out) & grepl("(^wt$|wild|control|ctrl)", x_low)] <- "WT"
  out[is.na(out) & grepl("(mut|r567w|homo|hom|ko|case)", x_low)] <- "MUT"
  out
}

pick_genotype_col <- function(md, preferred) {
  nms <- colnames(md)
  if (preferred %in% nms) return(preferred)
  fallbacks <- c("condition", "Condition", "group", "Group", "genotype", "Genotype")
  hit <- fallbacks[fallbacks %in% nms]
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

downsample_pt <- function(df, n_iter = 300L, late_q = 0.80, seed = 1L) {
  set.seed(seed)
  df <- df[!is.na(genotype2) & !is.na(pseudotime_transfer)]
  wt <- df[genotype2 == "WT"]
  mut <- df[genotype2 == "MUT"]
  n <- min(nrow(wt), nrow(mut))
  if (n < 10) return(data.table())
  late_thr <- as.numeric(quantile(df$pseudotime_transfer, late_q, na.rm = TRUE))
  out <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    wi <- wt[sample(.N, n)]
    mi <- mut[sample(.N, n)]
    out[[i]] <- data.table(
      iter = i,
      delta_median_mut_minus_wt = median(mi$pseudotime_transfer, na.rm = TRUE) - median(wi$pseudotime_transfer, na.rm = TRUE),
      delta_latefrac_mut_minus_wt = mean(mi$pseudotime_transfer >= late_thr, na.rm = TRUE) - mean(wi$pseudotime_transfer >= late_thr, na.rm = TRUE)
    )
  }
  rbindlist(out)
}

downsample_branch <- function(df, n_iter = 300L, seed = 1L) {
  set.seed(seed)
  df <- df[!is.na(genotype2) & !is.na(branch_transfer)]
  wt <- df[genotype2 == "WT"]
  mut <- df[genotype2 == "MUT"]
  n <- min(nrow(wt), nrow(mut))
  if (n < 10) return(data.table())
  branches <- sort(unique(df$branch_transfer))
  out <- vector("list", n_iter * max(1, length(branches)))
  k <- 0L
  for (i in seq_len(n_iter)) {
    wi <- wt[sample(.N, n)]
    mi <- mut[sample(.N, n)]
    for (b in branches) {
      k <- k + 1L
      out[[k]] <- data.table(
        iter = i,
        branch_transfer = b,
        delta_frac_mut_minus_wt = mean(mi$branch_transfer == b) - mean(wi$branch_transfer == b)
      )
    }
  }
  rbindlist(out[seq_len(k)])
}

out_res <- file.path(CFG$paths$out_dir, "results")
out_tbl <- file.path(out_res, "tables")
out_plt <- file.path(out_res, "plots")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plt, recursive = TRUE, showWarnings = FALSE)

extra_genes <- data.table(lineage = character(), gene = character())
if (!is.na(CFG$extra_genes_csv) && file.exists(CFG$extra_genes_csv)) {
  extra_genes <- fread(CFG$extra_genes_csv)
}

summary_rows <- list()
sr_i <- 0L
skip_rows <- list()
sk_i <- 0L

for (lineage_name in names(CFG$lineages)) {
  msg("Analyzing lineage: ", lineage_name)
  map_obj_path <- file.path(CFG$paths$out_dir, "mapping", "objects", paste0("09b_", lineage_name, "_mapped_query.rds"))
  if (!file.exists(map_obj_path)) {
    msg("  skip (missing mapped object): ", map_obj_path)
    next
  }

  seu <- readRDS(map_obj_path)
  md <- as.data.table(seu@meta.data, keep.rownames = "cell_id")

  geno_col <- pick_genotype_col(md, CFG$cols$genotype_col)
  req_cols <- c("pseudotime_transfer", "branch_transfer")
  if (!is.na(geno_col)) req_cols <- c(geno_col, req_cols)
  miss <- req_cols[!req_cols %in% colnames(md)]
  if (is.na(geno_col) || length(miss) > 0) {
    msg("  skip (missing columns): ", paste(miss, collapse = ", "))
    sk_i <- sk_i + 1L
    skip_rows[[sk_i]] <- data.table(lineage_name = lineage_name, reason = paste0("missing columns: ", paste(miss, collapse = ",")))
    next
  }

  md[, genotype2 := canon_genotype(get(geno_col), CFG$genotype$wt_values, CFG$genotype$mut_values)]
  df <- md[!is.na(genotype2)]
  if (nrow(df) == 0) {
    msg("  skip (no WT/MUT recognized after canonicalization from column ", geno_col, ").")
    sk_i <- sk_i + 1L
    skip_rows[[sk_i]] <- data.table(lineage_name = lineage_name, reason = paste0("no canonical WT/MUT from column: ", geno_col))
    next
  }

  n_by_g <- df[, .N, by = genotype2]
  if (nrow(n_by_g) < 2 || min(n_by_g$N, na.rm = TRUE) < CFG$analysis$min_cells_per_genotype) {
    msg("  warning: low cells in one genotype for lineage ", lineage_name)
  }

  # Pseudotime effect sizes
  med_mut <- median(df[genotype2 == "MUT"]$pseudotime_transfer, na.rm = TRUE)
  med_wt <- median(df[genotype2 == "WT"]$pseudotime_transfer, na.rm = TRUE)
  late_thr <- as.numeric(quantile(df$pseudotime_transfer, CFG$analysis$late_quantile, na.rm = TRUE))
  late_mut <- mean(df[genotype2 == "MUT"]$pseudotime_transfer >= late_thr, na.rm = TRUE)
  late_wt <- mean(df[genotype2 == "WT"]$pseudotime_transfer >= late_thr, na.rm = TRUE)
  wt_p <- tryCatch(wilcox.test(pseudotime_transfer ~ genotype2, data = df)$p.value, error = function(e) NA_real_)

  eff <- data.table(
    lineage_name = lineage_name,
    genotype_source_col = geno_col,
    n_cells = nrow(df),
    n_wt = sum(df$genotype2 == "WT"),
    n_mut = sum(df$genotype2 == "MUT"),
    median_mut = med_mut,
    median_wt = med_wt,
    delta_median_mut_minus_wt = med_mut - med_wt,
    late_thr = late_thr,
    late_frac_mut = late_mut,
    late_frac_wt = late_wt,
    delta_latefrac_mut_minus_wt = late_mut - late_wt,
    wilcox_p = wt_p
  )
  if ("ref_filter_label_col" %in% colnames(df)) eff[, ref_filter_label_col := as.character(df$ref_filter_label_col[1])]
  if ("ref_slingshot_cluster_col" %in% colnames(df)) eff[, ref_slingshot_cluster_col := as.character(df$ref_slingshot_cluster_col[1])]
  if ("ref_root_label" %in% colnames(df)) eff[, ref_root_label := as.character(df$ref_root_label[1])]
  if ("ref_start_cluster" %in% colnames(df)) eff[, ref_start_cluster := as.character(df$ref_start_cluster[1])]
  fwrite(eff, file.path(out_tbl, paste0("09c_", lineage_name, "_effect_sizes.csv")))

  # Violin/density
  p_v <- ggplot(df, aes(genotype2, pseudotime_transfer, fill = genotype2)) +
    geom_violin(trim = FALSE, alpha = 0.75) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    theme_classic(base_size = 12) +
    labs(title = paste0(lineage_name, ": transferred pseudotime by genotype"), x = NULL, y = "Transferred pseudotime")
  save_pdf(p_v, file.path(out_plt, paste0("09c_", lineage_name, "_pseudotime_violin")), w = 6.2, h = 4.4)

  p_d <- ggplot(df, aes(pseudotime_transfer, color = genotype2, fill = genotype2)) +
    geom_density(alpha = 0.2, linewidth = 1.0) +
    theme_classic(base_size = 12) +
    labs(title = paste0(lineage_name, ": transferred pseudotime density"), x = "Transferred pseudotime", y = "Density")
  save_pdf(p_d, file.path(out_plt, paste0("09c_", lineage_name, "_pseudotime_density")), w = 7.0, h = 4.2)

  # Downsampling robustness
  ds <- downsample_pt(df, n_iter = CFG$analysis$downsample_iterations, late_q = CFG$analysis$late_quantile)
  if (nrow(ds) > 0) {
    fwrite(ds, file.path(out_tbl, paste0("09c_", lineage_name, "_downsample_pseudotime.csv")))
    ds_sum <- data.table(
      lineage_name = lineage_name,
      delta_median_mean = mean(ds$delta_median_mut_minus_wt, na.rm = TRUE),
      delta_median_ci_low = quantile(ds$delta_median_mut_minus_wt, 0.025, na.rm = TRUE),
      delta_median_ci_high = quantile(ds$delta_median_mut_minus_wt, 0.975, na.rm = TRUE),
      delta_latefrac_mean = mean(ds$delta_latefrac_mut_minus_wt, na.rm = TRUE),
      delta_latefrac_ci_low = quantile(ds$delta_latefrac_mut_minus_wt, 0.025, na.rm = TRUE),
      delta_latefrac_ci_high = quantile(ds$delta_latefrac_mut_minus_wt, 0.975, na.rm = TRUE)
    )
    fwrite(ds_sum, file.path(out_tbl, paste0("09c_", lineage_name, "_downsample_pseudotime_summary.csv")))

    p_ds <- ggplot(ds, aes(delta_median_mut_minus_wt)) +
      geom_histogram(bins = 40, fill = "#4C72B0", alpha = 0.75) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_classic(base_size = 12) +
      labs(title = paste0(lineage_name, ": downsampling robustness (median delta)"), x = "Delta median (MUT - WT)", y = "Count")
    save_pdf(p_ds, file.path(out_plt, paste0("09c_", lineage_name, "_downsample_median_delta")), w = 7, h = 4.2)
  }

  # Branch allocation
  br <- df[, .(n_cells = .N), by = .(genotype2, branch_transfer)]
  br[, frac := n_cells / sum(n_cells), by = genotype2]
  fwrite(br, file.path(out_tbl, paste0("09c_", lineage_name, "_branch_allocation.csv")))

  p_br <- ggplot(br, aes(genotype2, frac, fill = branch_transfer)) +
    geom_col(position = "stack", width = 0.75) +
    theme_classic(base_size = 12) +
    labs(title = paste0(lineage_name, ": branch allocation by genotype"), x = NULL, y = "Fraction")
  save_pdf(p_br, file.path(out_plt, paste0("09c_", lineage_name, "_branch_allocation_bar")), w = 6.4, h = 4.4)

  dsb <- downsample_branch(df, n_iter = CFG$analysis$downsample_iterations)
  if (nrow(dsb) > 0) {
    fwrite(dsb, file.path(out_tbl, paste0("09c_", lineage_name, "_downsample_branch.csv")))
    dsb_sum <- dsb[, .(
      mean_delta = mean(delta_frac_mut_minus_wt, na.rm = TRUE),
      ci_low = quantile(delta_frac_mut_minus_wt, 0.025, na.rm = TRUE),
      ci_high = quantile(delta_frac_mut_minus_wt, 0.975, na.rm = TRUE)
    ), by = branch_transfer]
    dsb_sum[, lineage_name := lineage_name]
    fwrite(dsb_sum, file.path(out_tbl, paste0("09c_", lineage_name, "_downsample_branch_summary.csv")))
  }

  # Gene dynamics along transferred pseudotime
  genes_cfg <- unique(CFG$lineages[[lineage_name]]$marker_genes)
  genes_extra <- extra_genes[lineage == lineage_name]$gene
  genes <- unique(c(genes_cfg, genes_extra))
  genes <- genes[genes %in% rownames(seu)]

  if (length(genes) > 0) {
    assay_use <- if ("SCT" %in% names(seu@assays)) "SCT" else "RNA"
    expr <- fetch_expr_data(seu, vars = genes, assay = assay_use)
    expr <- as.data.table(expr, keep.rownames = "cell_id")
    gx <- merge(
      df[, .(cell_id, genotype2, pseudotime_transfer)],
      expr,
      by = "cell_id",
      all.x = FALSE,
      all.y = FALSE
    )
    long <- melt(
      gx,
      id.vars = c("cell_id", "genotype2", "pseudotime_transfer"),
      variable.name = "gene",
      value.name = "expr"
    )

    p_tr <- ggplot(long, aes(pseudotime_transfer, expr, color = genotype2)) +
      geom_point(size = 0.08, alpha = 0.12) +
      geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE, linewidth = 0.9) +
      facet_wrap(~ gene, scales = "free_y", ncol = 4) +
      theme_classic(base_size = 11) +
      labs(title = paste0(lineage_name, ": gene dynamics vs transferred pseudotime"), x = "Transferred pseudotime", y = "Expression")
    save_pdf(p_tr, file.path(out_plt, paste0("09c_", lineage_name, "_gene_trends")), w = 12.5, h = 8.5)

    if (isTRUE(CFG$analysis$run_mgcv) && requireNamespace("mgcv", quietly = TRUE)) {
      gstats <- vector("list", length(genes))
      gi <- 0L
      for (g in genes) {
        gi <- gi + 1L
        dfg <- long[gene == g & !is.na(expr) & !is.na(pseudotime_transfer)]
        if (nrow(dfg) < 120 || uniqueN(dfg$genotype2) < 2) {
          gstats[[gi]] <- data.table(gene = g, p_interaction = NA_real_)
          next
        }
        fit0 <- tryCatch(mgcv::gam(expr ~ genotype2 + s(pseudotime_transfer, k = 5), data = dfg), error = function(e) NULL)
        fit1 <- tryCatch(mgcv::gam(expr ~ genotype2 + s(pseudotime_transfer, k = 5) + s(pseudotime_transfer, by = genotype2, k = 5), data = dfg), error = function(e) NULL)
        if (is.null(fit0) || is.null(fit1)) {
          gstats[[gi]] <- data.table(gene = g, p_interaction = NA_real_)
          next
        }
        a <- tryCatch(anova(fit0, fit1, test = "Chisq"), error = function(e) NULL)
        p_int <- if (!is.null(a) && nrow(a) >= 2) as.numeric(a$`Pr(>Chi)`[2]) else NA_real_
        gstats[[gi]] <- data.table(gene = g, p_interaction = p_int)
      }
      gt <- rbindlist(gstats, fill = TRUE)
      gt[, padj := p.adjust(p_interaction, method = "BH")]
      setorder(gt, padj, p_interaction)
      fwrite(gt, file.path(out_tbl, paste0("09c_", lineage_name, "_mgcv_interaction_ranked.csv")))
    }
  }

  sr_i <- sr_i + 1L
  summary_rows[[sr_i]] <- eff
}

summary_tbl <- if (length(summary_rows) > 0) rbindlist(summary_rows, fill = TRUE) else data.table()
fwrite(summary_tbl, file.path(out_tbl, "09c_lineage_effect_summary.csv"))
if (length(skip_rows) > 0) {
  fwrite(rbindlist(skip_rows, fill = TRUE), file.path(out_tbl, "09c_lineage_skips.csv"))
}

# Markdown summary report
report_path <- file.path(out_res, "summary_report.md")
con <- file(report_path, open = "wt")
writeLines("# Pseudotime Mapping Summary\n", con)
writeLines("This report summarizes genotype shifts on transferred lineage pseudotime axes.\n", con)
if (nrow(summary_tbl) == 0) {
  writeLines("No lineage outputs were available.\n", con)
} else {
  for (i in seq_len(nrow(summary_tbl))) {
    x <- summary_tbl[i]
    writeLines(paste0("## ", x$lineage_name), con)
    writeLines(paste0("- Cells (WT/MUT): ", x$n_wt, " / ", x$n_mut), con)
    writeLines(paste0("- Delta median pseudotime (MUT - WT): ", signif(x$delta_median_mut_minus_wt, 4)), con)
    writeLines(paste0("- Delta late fraction (MUT - WT): ", signif(x$delta_latefrac_mut_minus_wt, 4)), con)
    writeLines(paste0("- Wilcoxon p-value: ", signif(x$wilcox_p, 4), "\n"), con)
    if ("ref_filter_label_col" %in% colnames(summary_tbl)) {
      writeLines(paste0("- Ref label column used: ", x$ref_filter_label_col), con)
    }
    if ("ref_slingshot_cluster_col" %in% colnames(summary_tbl)) {
      writeLines(paste0("- Slingshot cluster column used: ", x$ref_slingshot_cluster_col), con)
    }
    if ("ref_root_label" %in% colnames(summary_tbl)) {
      writeLines(paste0("- Root label/start cluster: ", x$ref_root_label, " / ", x$ref_start_cluster, "\n"), con)
    }
  }
}
close(con)

msg("09c complete. Outputs under: ", out_res)
