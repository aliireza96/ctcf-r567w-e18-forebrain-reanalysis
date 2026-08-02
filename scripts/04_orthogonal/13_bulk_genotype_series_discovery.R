# ============================================================
# 13_bulk_genotype_series_discovery.R — discovery-oriented genotype-series bulk RNA analysis
# ============================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])[1]
file_arg <- gsub("~\\+~", " ", file_arg)
SCRIPT_DIR <- dirname(normalizePath(ifelse(length(file_arg) == 0 || is.na(file_arg), ".", file_arg)))
source(file.path(SCRIPT_DIR, "12_discovery_setup.R"))
discovery_require(c(
  "data.table", "dplyr", "ggplot2", "ggrepel", "readr", "tidyr", "stringr",
  "limma", "pheatmap", "patchwork", "scales", "clusterProfiler", "org.Mm.eg.db"
))
discovery_append_log("bulk", "13_bulk_genotype_series_discovery")

sc_global <- load_scrna_global_de() %>%
  dplyr::select(gene, sc_avg_log2FC = avg_log2FC, sc_p_adj = p_val_adj)

infer_sex <- function(expr) {
  xist <- if ("Xist" %in% rownames(expr)) log2(expr["Xist", ] + 1) else rep(NA_real_, ncol(expr))
  y_genes <- intersect(c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d"), rownames(expr))
  y_score <- if (length(y_genes) > 0) {
    colMeans(log2(expr[y_genes, , drop = FALSE] + 1), na.rm = TRUE)
  } else {
    rep(NA_real_, ncol(expr))
  }
  male_score <- y_score - xist
  sex_call <- ifelse(is.na(male_score), "unknown", ifelse(male_score > 0, "male", "female"))
  data.frame(
    sample = colnames(expr),
    xist_log2 = xist,
    y_marker_log2 = y_score,
    male_score = male_score,
    sex_call = sex_call,
    stringsAsFactors = FALSE
  )
}

rename_contrast_cols <- function(tbl, prefix) {
  tbl %>%
    dplyr::select(gene, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
    dplyr::rename(
      !!paste0(prefix, "_logFC") := logFC,
      !!paste0(prefix, "_AveExpr") := AveExpr,
      !!paste0(prefix, "_t") := t,
      !!paste0(prefix, "_P.Value") := P.Value,
      !!paste0(prefix, "_adj.P.Val") := adj.P.Val
    )
}

fraction_string_to_numeric <- function(x) {
  sapply(strsplit(x, "/"), function(parts) {
    if (length(parts) != 2) {
      return(NA_real_)
    }
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
}

make_signature_plot <- function(sig_scores, tissue, focus_only = FALSE) {
  plot_df <- sig_scores
  if (focus_only) {
    plot_df <- plot_df %>%
      dplyr::filter(signature %in% c("SPN_LGE_identity", "Synaptic_maturation", "Immature_guidance", "MGE_interneuron"))
  }
  plot_df <- plot_df %>%
    dplyr::filter(is.finite(score), !is.na(genotype))
  if (nrow(plot_df) == 0) {
    return(
      ggplot() +
        theme_void() +
        annotate("text", x = 0, y = 0, label = paste0("No finite signature scores available for ", tissue, "."), size = 5)
    )
  }
  ggplot(plot_df, aes(x = genotype, y = score, color = genotype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.12, linewidth = 0.5) +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.3, alpha = 0.9) +
    facet_wrap(~signature, scales = "free_y", ncol = ifelse(focus_only, 2, 3)) +
    scale_color_manual(values = DISCOVERY_COLORS) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(tools::toTitleCase(tissue), " bulk RNA: genotype-series signature scores"),
      subtitle = "Scores are mean gene-wise z scores across each signature.",
      x = NULL,
      y = "Signature score",
      color = "Genotype"
    )
}

make_cpcdh_plot <- function(log_expr, sample_info, tissue) {
  cpcdh_genes <- resolve_cpcdh_genes(rownames(log_expr))
  if (length(cpcdh_genes) < 2) {
    return(NULL)
  }
  family_df <- data.frame(
    gene = cpcdh_genes,
    family = dplyr::case_when(
      grepl("^Pcdha", cpcdh_genes) ~ "Pcdha",
      grepl("^Pcdhb", cpcdh_genes) ~ "Pcdhb",
      grepl("^Pcdhg", cpcdh_genes) ~ "Pcdhg",
      TRUE ~ "Other"
    ),
    stringsAsFactors = FALSE
  )
  expr_long <- as.data.frame(log_expr[cpcdh_genes, , drop = FALSE]) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "log_expr") %>%
    dplyr::left_join(family_df, by = "gene") %>%
    dplyr::left_join(sample_info, by = "sample")
  family_summary <- expr_long %>%
    dplyr::group_by(sample, genotype, family) %>%
    dplyr::summarise(mean_log_expr = mean(log_expr, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(is.finite(mean_log_expr))
  if (nrow(family_summary) == 0) {
    return(NULL)
  }
  ggplot(family_summary, aes(x = genotype, y = mean_log_expr, color = genotype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.12, linewidth = 0.5) +
    geom_point(position = position_jitter(width = 0.10, height = 0), size = 2.0) +
    facet_wrap(~family, scales = "free_y") +
    scale_color_manual(values = DISCOVERY_COLORS) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(tools::toTitleCase(tissue), " cPcdh family expression"),
      subtitle = "Family scores use the mean log2(normalized count + 1) across detected cPcdh genes.",
      x = NULL,
      y = "Mean family log expression"
    )
}

make_concordance_plot <- function(merged, tissue) {
  highlight_genes <- unique(c(
    "Ebf1", "Meis2", "Bcl11b", "Grin2a", "Sparcl1", "Gria1",
    head(resolve_cpcdh_genes(merged$gene), 8)
  ))
  merged$label <- ifelse(merged$gene %in% highlight_genes, merged$gene, NA_character_)
  merged$signature_group <- ifelse(is.na(merged$signature_group), "background", merged$signature_group)
  corr_val <- safe_cor(merged$sc_avg_log2FC, merged$HOMO_vs_WT_logFC, method = "spearman")
  ggplot(merged, aes(x = sc_avg_log2FC, y = HOMO_vs_WT_logFC, color = signature_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(alpha = 0.65, size = 1.5) +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = 30, box.padding = 0.25, min.segment.length = 0) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(tools::toTitleCase(tissue), " bulk vs E18.5 snRNA concordance"),
      subtitle = paste0("Comparison uses HOMO vs WT bulk logFC and mut vs wt single-cell global logFC. Spearman rho = ", sprintf("%.2f", corr_val)),
      x = "E18.5 snRNA global log2FC (mut vs wt)",
      y = "Bulk RNA log2FC (HOMO vs WT)",
      color = "Signature"
    )
}

analyze_tissue <- function(tissue) {
  discovery_msg("Running discovery bulk analysis for ", tissue, ".")
  counts_path <- file.path(BULK_RAW, paste0("GSE214689_", tissue, ".nor_counts.tsv.gz"))
  expr <- read_bulk_matrix(counts_path)
  sample_info <- parse_bulk_sample_info(colnames(expr), default_tissue = tissue)
  sample_info$tissue <- tissue
  sample_info <- sample_info %>%
    dplyr::left_join(infer_sex(expr), by = c("sample"))

  write_csv_discovery("bulk", sample_info, paste0("13_", tissue, "_sample_metadata.csv"))

  log_expr <- log2(expr + 1)
  var_genes <- names(sort(matrixStats::rowVars(log_expr), decreasing = TRUE))[seq_len(min(2000, nrow(log_expr)))]
  pca <- stats::prcomp(t(log_expr[var_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
  pca_df <- cbind(sample_info, as.data.frame(pca$x[, 1:4, drop = FALSE]))
  explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = genotype, shape = sex_call, label = sample)) +
    geom_point(size = 3.2) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, box.padding = 0.2, min.segment.length = 0) +
    scale_color_manual(values = DISCOVERY_COLORS) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0(tools::toTitleCase(tissue), " bulk RNA genotype-series PCA"),
      subtitle = "PCA uses the top 2,000 variable genes after log2(normalized count + 1) transform.",
      x = paste0("PC1 (", explained[1], "%)"),
      y = paste0("PC2 (", explained[2], "%)"),
      color = "Genotype",
      shape = "Sex call"
    )
  save_plot_discovery("bulk", p_pca, paste0("13_", tissue, "_pca"), w = 8.6, h = 6.0)

  annotation_col <- data.frame(
    Genotype = sample_info$genotype,
    Sex = sample_info$sex_call,
    row.names = sample_info$sample,
    stringsAsFactors = FALSE
  )
  sample_dist <- 1 - stats::cor(log_expr[var_genes, , drop = FALSE], method = "spearman")
  save_pheatmap_discovery(
    "bulk",
    paste0("13_", tissue, "_sample_distance"),
    sample_dist,
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    main = paste0(tools::toTitleCase(tissue), " sample distance (1 - Spearman)")
  )

  design_obj <- build_bulk_design(sample_info)
  discovery_msg("Using design for ", tissue, ": ", design_obj$design_label)
  fit <- limma::lmFit(log_expr, design_obj$design)
  contrast_matrix <- limma::makeContrasts(
    HET_vs_WT = genotypeHET - genotypeWT,
    HOMO_vs_WT = genotypeHOMO - genotypeWT,
    HOMO_vs_HET = genotypeHOMO - genotypeHET,
    levels = design_obj$design
  )
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE)

  de_tables <- lapply(colnames(contrast_matrix), function(contrast_name) {
    tt <- limma::topTable(fit2, coef = contrast_name, number = Inf, sort.by = "P")
    tt$gene <- rownames(tt)
    tt$contrast <- contrast_name
    tt$tissue <- tissue
    tt
  })
  names(de_tables) <- colnames(contrast_matrix)

  invisible(lapply(names(de_tables), function(contrast_name) {
    write_csv_discovery("bulk", de_tables[[contrast_name]], paste0("13_", tissue, "_", contrast_name, ".csv"))
  }))

  contrast_summary <- dplyr::bind_rows(lapply(names(de_tables), function(contrast_name) {
    tbl <- de_tables[[contrast_name]]
    data.frame(
      tissue = tissue,
      contrast = contrast_name,
      design = design_obj$design_label,
      n_sig = sum(tbl$adj.P.Val < 0.05 & abs(tbl$logFC) >= 0.5, na.rm = TRUE),
      n_up = sum(tbl$adj.P.Val < 0.05 & tbl$logFC > 0.5, na.rm = TRUE),
      n_down = sum(tbl$adj.P.Val < 0.05 & tbl$logFC < -0.5, na.rm = TRUE),
      median_abs_logFC = median(abs(tbl$logFC), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write_csv_discovery("bulk", contrast_summary, paste0("13_", tissue, "_contrast_summary.csv"))

  group_means <- sapply(levels(design_obj$sample_info$genotype), function(gt) {
    rowMeans(log_expr[, design_obj$sample_info$genotype == gt, drop = FALSE], na.rm = TRUE)
  })
  trend_df <- data.frame(
    gene = rownames(log_expr),
    mean_WT = group_means[, "WT"],
    mean_HET = group_means[, "HET"],
    mean_HOMO = group_means[, "HOMO"],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      delta_HET_vs_WT = mean_HET - mean_WT,
      delta_HOMO_vs_WT = mean_HOMO - mean_WT,
      delta_HOMO_vs_HET = mean_HOMO - mean_HET,
      trend_class = mapply(classify_dosage_pattern, mean_WT, mean_HET, mean_HOMO)
    ) %>%
    dplyr::left_join(rename_contrast_cols(de_tables$HET_vs_WT, "HET_vs_WT"), by = "gene") %>%
    dplyr::left_join(rename_contrast_cols(de_tables$HOMO_vs_WT, "HOMO_vs_WT"), by = "gene") %>%
    dplyr::left_join(rename_contrast_cols(de_tables$HOMO_vs_HET, "HOMO_vs_HET"), by = "gene")

  write_csv_discovery("bulk", trend_df, paste0("13_", tissue, "_dosage_patterns.csv"))

  candidate_genes <- unique(c(DISCOVERY_CANDIDATE_GENES, resolve_cpcdh_genes(rownames(log_expr))))
  candidate_summary <- trend_df %>%
    dplyr::filter(gene %in% candidate_genes) %>%
    dplyr::arrange(dplyr::desc(abs(HOMO_vs_WT_logFC)), HOMO_vs_WT_adj.P.Val)
  write_csv_discovery("bulk", candidate_summary, paste0("13_", tissue, "_candidate_gene_dosage_summary.csv"))

  signature_scores <- dplyr::bind_rows(lapply(names(DISCOVERY_SIGNATURES), function(sig_name) {
    sig_genes <- DISCOVERY_SIGNATURES[[sig_name]]
    if (sig_name == "cPcdh_core") {
      sig_genes <- resolve_cpcdh_genes(rownames(log_expr))
    }
    present <- intersect(sig_genes, rownames(log_expr))
    data.frame(
      sample = colnames(log_expr),
      signature = sig_name,
      score = score_signature(log_expr, present),
      n_genes = length(present),
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::left_join(sample_info, by = "sample") %>%
    dplyr::mutate(tissue = tissue)
  write_csv_discovery("bulk", signature_scores, paste0("13_", tissue, "_signature_scores.csv"))

  save_plot_discovery("bulk", make_signature_plot(signature_scores, tissue), paste0("13_", tissue, "_signature_scores"), w = 11.0, h = 8.0)

  cpcdh_plot <- make_cpcdh_plot(log_expr, sample_info, tissue)
  if (!is.null(cpcdh_plot)) {
    save_plot_discovery("bulk", cpcdh_plot, paste0("13_", tissue, "_cpcdh_family_expression"), w = 8.5, h = 4.8)
  }

  merged <- trend_df %>%
    dplyr::select(gene, HOMO_vs_WT_logFC, HOMO_vs_WT_adj.P.Val, trend_class) %>%
    dplyr::inner_join(sc_global, by = "gene")

  signature_map_rows <- lapply(names(DISCOVERY_SIGNATURES), function(sig_name) {
    sig_genes <- DISCOVERY_SIGNATURES[[sig_name]]
    if (sig_name == "cPcdh_core") {
      sig_genes <- resolve_cpcdh_genes(merged$gene)
    }
    present <- intersect(sig_genes, merged$gene)
    if (length(present) == 0) {
      return(data.frame())
    }
    data.frame(gene = present, signature_group = rep(sig_name, length(present)), stringsAsFactors = FALSE)
  })
  signature_map_rows <- signature_map_rows[vapply(signature_map_rows, nrow, integer(1)) > 0]
  signature_map <- if (length(signature_map_rows) == 0) {
    data.frame(gene = character(), signature_group = character(), stringsAsFactors = FALSE)
  } else {
    dplyr::bind_rows(signature_map_rows) %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(signature_group = paste(unique(signature_group), collapse = ";"), .groups = "drop")
  }

  merged <- merged %>% dplyr::left_join(signature_map, by = "gene")
  write_csv_discovery("bulk", merged, paste0("13_", tissue, "_bulk_vs_scrna_concordance.csv"))
  save_plot_discovery("bulk", make_concordance_plot(merged, tissue), paste0("13_", tissue, "_bulk_vs_scrna_concordance"), w = 9.2, h = 7.0)

  signature_summary <- dplyr::bind_rows(lapply(names(DISCOVERY_SIGNATURES), function(sig_name) {
    sig_genes <- DISCOVERY_SIGNATURES[[sig_name]]
    if (sig_name == "cPcdh_core") {
      sig_genes <- resolve_cpcdh_genes(merged$gene)
    }
    sig_genes <- intersect(sig_genes, merged$gene)
    if (length(sig_genes) == 0) {
      return(data.frame())
    }
    sig_tbl <- merged %>% dplyr::filter(gene %in% sig_genes)
    data.frame(
      tissue = tissue,
      signature = sig_name,
      n_genes = nrow(sig_tbl),
      bulk_mean_logFC = mean(sig_tbl$HOMO_vs_WT_logFC, na.rm = TRUE),
      scrna_mean_logFC = mean(sig_tbl$sc_avg_log2FC, na.rm = TRUE),
      bulk_vs_scrna_cor = safe_cor(sig_tbl$HOMO_vs_WT_logFC, sig_tbl$sc_avg_log2FC, method = "spearman"),
      bulk_sig_fraction = mean(sig_tbl$HOMO_vs_WT_adj.P.Val < 0.05 & abs(sig_tbl$HOMO_vs_WT_logFC) >= 0.5, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write_csv_discovery("bulk", signature_summary, paste0("13_", tissue, "_signature_direction_summary.csv"))

  go_sets <- list(
    HOMO_vs_WT_UP = de_tables$HOMO_vs_WT %>% dplyr::filter(adj.P.Val < 0.05, logFC > 0.5) %>% dplyr::pull(gene),
    HOMO_vs_WT_DOWN = de_tables$HOMO_vs_WT %>% dplyr::filter(adj.P.Val < 0.05, logFC < -0.5) %>% dplyr::pull(gene),
    HET_vs_WT_UP = de_tables$HET_vs_WT %>% dplyr::filter(adj.P.Val < 0.05, logFC > 0.5) %>% dplyr::pull(gene),
    HET_vs_WT_DOWN = de_tables$HET_vs_WT %>% dplyr::filter(adj.P.Val < 0.05, logFC < -0.5) %>% dplyr::pull(gene),
    monotonic_up = trend_df %>% dplyr::filter(trend_class == "monotonic_up") %>% dplyr::pull(gene),
    homo_threshold_down = trend_df %>% dplyr::filter(trend_class == "homo_threshold_down") %>% dplyr::pull(gene)
  )
  go_results <- dplyr::bind_rows(lapply(names(go_sets), function(set_name) {
    ego <- run_go_enrichment(go_sets[[set_name]], rownames(log_expr))
    if (nrow(ego) == 0) {
      return(data.frame())
    }
    ego$set <- set_name
    ego$tissue <- tissue
    ego
  }))
  if (nrow(go_results) > 0) {
    write_csv_discovery("bulk", go_results, paste0("13_", tissue, "_GO_enrichment.csv"))
    top_terms <- go_results %>%
      dplyr::group_by(set) %>%
      dplyr::slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(GeneRatioNum = fraction_string_to_numeric(GeneRatio))
    p_go <- ggplot(top_terms, aes(x = GeneRatioNum, y = reorder(Description, GeneRatioNum), size = Count, color = p.adjust)) +
      geom_point() +
      facet_wrap(~set, scales = "free_y") +
      scale_color_viridis_c(option = "magma", direction = -1) +
      theme_classic(base_size = 11) +
      labs(
        title = paste0(tools::toTitleCase(tissue), " GO enrichment across key genotype-series gene sets"),
        x = "Gene ratio",
        y = NULL,
        color = "BH-adjusted p",
        size = "Count"
      )
    save_plot_discovery("bulk", p_go, paste0("13_", tissue, "_go_summary"), w = 12.0, h = 8.2)
  }

  top_heat_genes <- trend_df %>%
    dplyr::filter(HOMO_vs_WT_adj.P.Val < 0.05) %>%
    dplyr::arrange(dplyr::desc(abs(HOMO_vs_WT_logFC))) %>%
    dplyr::slice_head(n = 50) %>%
    dplyr::pull(gene)
  if (length(top_heat_genes) >= 2) {
    heat_mat <- zscore_rows(log_expr[top_heat_genes, sample_info$sample, drop = FALSE])
    save_pheatmap_discovery(
      "bulk",
      paste0("13_", tissue, "_top_dosage_heatmap"),
      heat_mat,
      annotation_col = annotation_col[sample_info$sample, , drop = FALSE],
      show_rownames = TRUE,
      main = paste0(tools::toTitleCase(tissue), " top HOMO vs WT genotype-series genes")
    )
  }

  if (tissue == "brain") {
    brain_focus_plot <- make_signature_plot(signature_scores, tissue, focus_only = TRUE) + p_pca + patchwork::plot_layout(widths = c(1.5, 1))
    save_plot_discovery("bulk", brain_focus_plot, "13_brain_genotype_series_biology_figure", w = 13.2, h = 7.0)
  }

  trend_counts <- trend_df %>%
    dplyr::count(trend_class, sort = TRUE) %>%
    dplyr::mutate(tissue = tissue)
  write_csv_discovery("bulk", trend_counts, paste0("13_", tissue, "_trend_class_counts.csv"))

  summary_lines <- c(
    paste0("# ", tools::toTitleCase(tissue), " bulk RNA discovery summary"),
    "",
    paste0("- Design used: ", design_obj$design_label),
    paste0("- Significant HOMO vs WT genes (|logFC| >= 0.5, BH < 0.05): ", contrast_summary$n_sig[contrast_summary$contrast == "HOMO_vs_WT"]),
    paste0("- Significant HET vs WT genes (|logFC| >= 0.5, BH < 0.05): ", contrast_summary$n_sig[contrast_summary$contrast == "HET_vs_WT"]),
    paste0("- Top dosage pattern class: ", trend_counts$trend_class[1], " (n = ", trend_counts$n[1], ")"),
    paste0("- Bulk vs snRNA HOMO-vs-WT concordance (Spearman): ", sprintf("%.2f", safe_cor(merged$HOMO_vs_WT_logFC, merged$sc_avg_log2FC)))
  )
  write_status_note("bulk", paste0("13_", tissue, "_summary.md"), summary_lines)

  list(
    tissue = tissue,
    contrast_summary = contrast_summary,
    signature_summary = signature_summary,
    trend_counts = trend_counts,
    candidate_summary = candidate_summary
  )
}

tissues <- c("brain", "neuron", "organoid")
results <- lapply(tissues, analyze_tissue)

contrast_summary_all <- dplyr::bind_rows(lapply(results, `[[`, "contrast_summary"))
signature_summary_all <- dplyr::bind_rows(lapply(results, `[[`, "signature_summary"))
trend_counts_all <- dplyr::bind_rows(lapply(results, `[[`, "trend_counts"))
candidate_summary_all <- dplyr::bind_rows(lapply(results, function(x) dplyr::mutate(x$candidate_summary, tissue = x$tissue)))

write_csv_discovery("bulk", contrast_summary_all, "13_all_tissues_contrast_summary.csv")
write_csv_discovery("bulk", signature_summary_all, "13_all_tissues_signature_direction_summary.csv")
write_csv_discovery("bulk", trend_counts_all, "13_all_tissues_trend_class_counts.csv")
write_csv_discovery("bulk", candidate_summary_all, "13_all_tissues_candidate_gene_dosage_summary.csv")

discovery_msg("Bulk genotype-series discovery analysis complete.")
