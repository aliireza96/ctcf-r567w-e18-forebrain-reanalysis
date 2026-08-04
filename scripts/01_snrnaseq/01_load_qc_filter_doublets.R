# ============================================================
# 01_load_qc_filter_doublets.R — load, QC and filtering
# ============================================================

source("scripts/00_setup.R")

data_dirs <- list(
  wt  = file.path(DATA_ROOT, "wt_10x"),
  mut = file.path(DATA_ROOT, "mut_10x")
)

stopifnot(all(sapply(data_dirs, dir.exists)))

objs <- list()

for (nm in names(data_dirs)) {
  message("Loading: ", nm)
  counts <- Read10X(data.dir = data_dirs[[nm]])
  
  so <- CreateSeuratObject(counts = counts, project = nm, min.cells = 3, min.features = 200)
  so$condition <- nm
  so$library_id <- if (nm == "wt") "ZJ6" else "ZJ9"
  so$pipeline <- "BD Rhapsody WTA v1.12"
  
  objs[[nm]] <- so
}

seu <- merge(objs$wt, y = objs$mut, add.cell.ids = names(objs), project = "E18p5_Cortex")
rm(objs); gc()

# ---- QC metrics ----
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^mt-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^(Rpl|Rps)")
# optional:
# seu[["percent.hb"]] <- PercentageFeatureSet(seu, pattern = "^Hb[ab]-")

# ---- QC plots (pre-filter) ----
p_vln <- VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
                 group.by = "condition", pt.size = 0.1, ncol = 4) +
  plot_annotation(title = "QC metrics before filtering")
save_plot(p_vln, "qc_plots/QC_violin_before", w=14, h=5)

p_sc1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("nCount vs %mt (before)")
p_sc2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("nCount vs nFeature (before)")
save_plot(p_sc1 + p_sc2, "qc_plots/QC_scatter_before", w=12, h=5)

# ---- Filtering thresholds (start conservative; adjust after looking at plots) ----
# NOTE: since WT and MUT are different pooled libraries, keep consistent thresholds.
min_features <- 200
max_features <- 6000
min_counts   <- 2000
max_mt       <- 20

seu_f <- subset(seu, subset =
                  nFeature_RNA > min_features &
                  nFeature_RNA < max_features &
                  nCount_RNA   > min_counts   &
                  percent.mt   < max_mt
)

# ---- Summary table ----
qc_summary <- data.frame(
  stage = c("before","after"),
  n_cells = c(ncol(seu), ncol(seu_f)),
  wt_cells = c(sum(seu$condition=="wt"), sum(seu_f$condition=="wt")),
  mut_cells= c(sum(seu$condition=="mut"), sum(seu_f$condition=="mut"))
)
write.csv(qc_summary, file.path(OUT$tables, "QC_cell_counts_before_after.csv"), row.names = FALSE)

# ---- QC plots (post-filter) ----
p_vln2 <- VlnPlot(seu_f, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"),
                  group.by = "condition", pt.size = 0.1, ncol = 4) +
  plot_annotation(title = "QC metrics after filtering")
save_plot(p_vln2, "qc_plots/QC_violin_after", w=14, h=5)

# ---- Sex marker QC (infer sample sex from expression) ----
# Use direct normalized expression summaries so the axes remain non-negative.
sex_female <- c("Xist")
sex_male   <- c("Ddx3y", "Eif2s3y", "Kdm5d", "Uty", "Rps4y1", "Zfy", "Sry")

match_genes <- function(genes, pool) {
  pool_lower <- tolower(pool)
  out <- c()
  for (g in genes) {
    idx <- which(pool_lower == tolower(g))
    if (length(idx) > 0) out <- c(out, pool[idx[1]])
  }
  unique(out)
}

sex_female <- match_genes(sex_female, rownames(seu_f))
sex_male   <- match_genes(sex_male, rownames(seu_f))

if (length(sex_female) + length(sex_male) >= 2) {
  DefaultAssay(seu_f) <- "RNA"
  seu_f <- NormalizeData(seu_f, verbose = FALSE)

  expr_cols <- unique(c(sex_female, sex_male))
  expr_df <- FetchData(seu_f, vars = expr_cols, assay = "RNA", layer = "data")
  expr_df <- as.data.frame(expr_df)

  if (length(sex_female) >= 1) {
    seu_f$sex_female_expr <- if (length(sex_female) == 1) {
      as.numeric(expr_df[[sex_female[1]]])
    } else {
      rowMeans(expr_df[, sex_female, drop = FALSE], na.rm = TRUE)
    }
  }
  if (length(sex_male) >= 1) {
    seu_f$sex_male_expr <- if (length(sex_male) == 1) {
      as.numeric(expr_df[[sex_male[1]]])
    } else {
      rowMeans(expr_df[, sex_male, drop = FALSE], na.rm = TRUE)
    }
  }

  # DotPlot of sex markers by condition and by library_id
  sex_genes <- c(sex_female, sex_male)
  p_sex_dot_cond <- DotPlot(seu_f, features = sex_genes, group.by = "condition", scale = FALSE) +
    RotatedAxis() + ggtitle("Sex marker expression (by condition)")
  save_plot(p_sex_dot_cond, "qc_plots/SexMarkers_Dotplot_by_condition", w=10, h=5)
  
  p_sex_dot_lib <- DotPlot(seu_f, features = sex_genes, group.by = "library_id", scale = FALSE) +
    RotatedAxis() + ggtitle("Sex marker expression (by library)")
  save_plot(p_sex_dot_lib, "qc_plots/SexMarkers_Dotplot_by_library", w=10, h=5)
  
  # Violin plots of normalized expression summaries
  score_feats <- c()
  if ("sex_female_expr" %in% colnames(seu_f@meta.data)) score_feats <- c(score_feats, "sex_female_expr")
  if ("sex_male_expr" %in% colnames(seu_f@meta.data))   score_feats <- c(score_feats, "sex_male_expr")

  if (length(score_feats) > 0) {
    p_sex_vln <- VlnPlot(seu_f, features = score_feats, group.by = "condition", pt.size = 0.1) +
      plot_annotation(title = "Sex marker expression summaries (by condition)")
    save_plot(p_sex_vln, "qc_plots/SexScores_Violin_by_condition", w=8, h=4)
  }

  # Scatter of Y-linked expression vs Xist expression using non-negative normalized values.
  if (all(c("sex_female_expr", "sex_male_expr", "condition") %in% colnames(seu_f@meta.data))) {
    sex_scatter_df <- data.frame(
      sex_male_expr = pmax(seu_f$sex_male_expr, 0),
      sex_female_expr = pmax(seu_f$sex_female_expr, 0),
      condition = seu_f$condition
    )
    p_sex_scatter <- ggplot(sex_scatter_df, aes(x = sex_male_expr, y = sex_female_expr, color = condition)) +
      geom_point(size = 0.3, alpha = 0.5) +
      facet_wrap(~condition) +
      scale_x_continuous(limits = c(0, NA)) +
      scale_y_continuous(limits = c(0, NA)) +
      labs(
        title = "Sex marker expression: Y-linked mean vs Xist",
        x = "Mean Y-linked normalized expression",
        y = "Xist normalized expression"
      ) +
      theme_bw(base_size = 11)
    save_plot(p_sex_scatter, "qc_plots/SexScores_Scatter_male_vs_female", w=7, h=5)
  }
} else {
  message("Sex marker QC skipped: too few sex-linked genes found in this dataset.")
}

# ---- Doublet-detection availability audit ----
# No algorithmic doublet detection or removal is applied in this workflow.
if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
  message("DoubletFinder is unavailable; no algorithmic doublet detection was applied.")
} else {
  message("DoubletFinder is available, but no algorithmic doublet detection was applied.")
}

# Save object
saveRDS(seu_f, file.path(OUT$objects, "01_seu_QCfiltered.rds"))
message("Saved: 01_seu_QCfiltered.rds")
