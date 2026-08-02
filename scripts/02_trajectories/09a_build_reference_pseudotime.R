#!/usr/bin/env Rscript

# ============================================================
# 09a_build_reference_pseudotime.R
# Build reference-only lineage pseudotime/branch labels
# using slingshot (ventral + dorsal references separately).
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(slingshot)
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

stopifnot(file.exists(CFG$paths$ref_ventral_rds))
stopifnot(file.exists(CFG$paths$ref_dorsal_rds))

out_refs <- file.path(CFG$paths$out_dir, "refs")
out_obj <- file.path(out_refs, "objects")
out_tbl <- file.path(out_refs, "tables")
out_plt <- file.path(out_refs, "plots")
dir.create(out_obj, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plt, recursive = TRUE, showWarnings = FALSE)

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
  FetchData(seu, vars = vars, assay = assay, slot = "data")
}

ensure_reductions <- function(seu, pca_npcs = 50L, umap_dims = 1:30) {
  if (!"pca" %in% names(seu@reductions)) {
    if ("SCT" %in% names(seu@assays)) {
      DefaultAssay(seu) <- "SCT"
    } else {
      DefaultAssay(seu) <- "RNA"
      seu <- NormalizeData(seu, verbose = FALSE)
      seu <- FindVariableFeatures(seu, verbose = FALSE)
      seu <- ScaleData(seu, verbose = FALSE)
    }
    seu <- RunPCA(seu, npcs = pca_npcs, verbose = FALSE)
  }
  use_dims <- seq_len(min(max(umap_dims), ncol(Embeddings(seu, "pca"))))
  rerun_umap <- !("umap" %in% names(seu@reductions))
  if (!rerun_umap) {
    um <- seu@reductions$umap
    has_model <- tryCatch(!is.null(um@misc$model), error = function(e) FALSE)
    rerun_umap <- !has_model
  }
  if (rerun_umap) {
    seu <- RunUMAP(seu, reduction = "pca", dims = use_dims, return.model = TRUE, verbose = FALSE)
  }
  seu
}

run_lineage_reference <- function(ref_obj, lineage_name, spec) {
  label_col <- spec$ref_label_col
  ref_name <- spec$reference_name %||% NA_character_

  if (!label_col %in% colnames(ref_obj@meta.data)) {
    if (identical(ref_name, "ventral") && "region" %in% colnames(ref_obj@meta.data)) {
      msg("[", lineage_name, "] ref label column '", label_col, "' missing; using 'region' instead.")
      label_col <- "region"
    } else if (identical(ref_name, "dorsal")) {
      # Dorsal fallback should remain dorsal-like, never region-based.
      dorsal_candidates <- c(
        label_col,
        CFG$cols$dorsal_label_col %||% NA_character_,
        "celltype_broad_dorsal",
        "Gral_cellType",
        "New_cellType"
      )
      dorsal_candidates <- unique(dorsal_candidates[!is.na(dorsal_candidates)])
      hit <- dorsal_candidates[dorsal_candidates %in% colnames(ref_obj@meta.data)]
      if (length(hit) > 0) {
        label_col <- hit[1]
        msg("[", lineage_name, "] ref label column '", spec$ref_label_col, "' missing; using '", label_col, "' for dorsal reference.")
      }
    }
  }

  label_map <- spec$ref_label_map %||% NULL
  if (!is.null(label_map) && label_col %in% colnames(ref_obj@meta.data)) {
    raw_lab <- as.character(ref_obj@meta.data[[label_col]])
    mapped <- unname(label_map[raw_lab])
    mapped[is.na(mapped)] <- raw_lab[is.na(mapped)]
    work_col <- paste0(".__", lineage_name, "_label_working")
    ref_obj@meta.data[[work_col]] <- mapped
    label_col <- work_col
    msg("[", lineage_name, "] using mapped reference labels from ", spec$ref_label_col, ".")
  }

  keep_labels <- spec$ref_keep_labels
  root_label <- spec$root_label
  root_cluster <- spec$root_cluster %||% NA_character_

  if (!label_col %in% colnames(ref_obj@meta.data)) {
    stop("[", lineage_name, "] ref label column missing: ", label_col)
  }

  # If dorsal fallback used Gral/New labels, map to coarse dorsal classes.
  if (identical(ref_name, "dorsal") && label_col %in% c("Gral_cellType", "New_cellType")) {
    raw_lab <- as.character(ref_obj@meta.data[[label_col]])
    mapped <- ifelse(
      grepl("Apical|Radial|RG|progenitor", raw_lab, ignore.case = TRUE), "RG",
      ifelse(
        grepl("Intermediate|IPC|Eomes|Tbr2", raw_lab, ignore.case = TRUE), "IPC",
        ifelse(
          grepl("Excit|Neuron", raw_lab, ignore.case = TRUE), "Excitatory",
          raw_lab
        )
      )
    )
    ref_obj$.__dorsal_label_working <- mapped
    label_col <- ".__dorsal_label_working"
    msg("[", lineage_name, "] mapped ", spec$ref_label_col, " to coarse dorsal classes (RG/IPC/Excitatory).")
  }

  # Map lineage labels to available labels when config uses descriptive aliases.
  keep_labels <- as.character(keep_labels)
  root_label <- as.character(root_label)
  available_labels <- unique(as.character(ref_obj@meta.data[[label_col]]))
  if (length(intersect(keep_labels, available_labels)) == 0) {
    mapped <- character(0)
    for (k in keep_labels) {
      if (grepl("LGE", k, ignore.case = TRUE) && "LGE" %in% available_labels) mapped <- c(mapped, "LGE")
      if (grepl("MGE", k, ignore.case = TRUE) && "MGE" %in% available_labels) mapped <- c(mapped, "MGE")
      if (grepl("CGE", k, ignore.case = TRUE) && "CGE" %in% available_labels) mapped <- c(mapped, "CGE")
    }
    mapped <- unique(mapped)
    if (length(mapped) > 0) {
      msg("[", lineage_name, "] remapped keep_labels -> ", paste(mapped, collapse = ", "))
      keep_labels <- mapped
    }
  }
  if (!root_label %in% available_labels) {
    if (grepl("LGE", root_label, ignore.case = TRUE) && "LGE" %in% available_labels) root_label <- "LGE"
    if (grepl("MGE", root_label, ignore.case = TRUE) && "MGE" %in% available_labels) root_label <- "MGE"
    if (grepl("CGE", root_label, ignore.case = TRUE) && "CGE" %in% available_labels) root_label <- "CGE"
  }

  keep_cells <- rownames(ref_obj@meta.data)[as.character(ref_obj@meta.data[[label_col]]) %in% keep_labels]
  if (length(keep_cells) < 200) {
    stop("[", lineage_name, "] too few cells after ref_keep_labels filtering: ", length(keep_cells))
  }

  seu <- subset(ref_obj, cells = keep_cells)
  seu <- ensure_reductions(seu)

  if (!root_label %in% unique(as.character(seu@meta.data[[label_col]]))) {
    stop("[", lineage_name, "] root_label not present after filtering: ", root_label)
  }

  pca <- Embeddings(seu, "pca")
  use_dims <- seq_len(min(50L, ncol(pca)))
  rd <- pca[, use_dims, drop = FALSE]
  # Use a richer cluster label for slingshot when the lineage label has <2 groups.
  cluster_col <- spec$trajectory_cluster_col %||% label_col
  if (!cluster_col %in% colnames(seu@meta.data) || uniqueN(seu@meta.data[[cluster_col]]) < 2) {
    if ("seurat_clusters" %in% colnames(seu@meta.data) && uniqueN(seu@meta.data$seurat_clusters) >= 2) {
      cluster_col <- "seurat_clusters"
      msg("[", lineage_name, "] using seurat_clusters for slingshot clusters.")
    } else {
      # Fallback deterministic kmeans over PCA for datasets lacking cluster metadata.
      km_k <- max(2L, min(6L, floor(sqrt(nrow(rd) / 200))))
      km <- kmeans(rd[, seq_len(min(10L, ncol(rd))), drop = FALSE], centers = km_k, nstart = 10)
      seu$kmeans_cluster_for_slingshot <- as.factor(km$cluster)
      cluster_col <- "kmeans_cluster_for_slingshot"
      msg("[", lineage_name, "] using kmeans clusters for slingshot (k=", km_k, ").")
    }
  }

  cl <- as.factor(as.character(seu@meta.data[[cluster_col]]))
  names(cl) <- colnames(seu)

  root_cells <- rownames(seu@meta.data)[as.character(seu@meta.data[[label_col]]) == root_label]
  start_level <- root_label
  if (!is.na(root_cluster) && nzchar(root_cluster)) {
    start_level <- as.character(root_cluster)
    if (!start_level %in% levels(cl)) {
      stop("[", lineage_name, "] requested root_cluster not present in ", cluster_col, ": ", start_level)
    }
  } else if (!start_level %in% levels(cl)) {
    # Pick starting slingshot cluster among root-labeled cells by lowest median PC1 (immature end heuristic).
    root_dt <- data.table(
      cell_id = rownames(seu@meta.data),
      sl_cluster = as.character(seu@meta.data[[cluster_col]]),
      is_root = rownames(seu@meta.data) %in% root_cells,
      PC_1 = rd[, 1]
    )
    cand <- root_dt[is_root == TRUE, .(med_pc1 = median(PC_1, na.rm = TRUE)), by = sl_cluster]
    if (nrow(cand) == 0) {
      cand <- root_dt[, .(med_pc1 = median(PC_1, na.rm = TRUE)), by = sl_cluster]
    }
    start_level <- cand[order(med_pc1)][1, sl_cluster]
  }

  msg("[", lineage_name, "] slingshot start = ", start_level, " (root_label=", root_label, "), cells = ", ncol(seu))
  sds <- slingshot(rd, clusterLabels = cl, start.clus = start_level)

  pt <- slingPseudotime(sds)
  cw <- slingCurveWeights(sds)
  if (is.null(dim(pt))) {
    pt <- matrix(pt, ncol = 1)
    colnames(pt) <- "Lineage1"
  }
  if (is.null(dim(cw))) {
    cw <- matrix(cw, ncol = 1)
    colnames(cw) <- "Lineage1"
  }

  branch_id <- apply(cw, 1, function(x) {
    if (all(is.na(x))) return(NA_integer_)
    which.max(x)
  })
  branch_id[is.na(branch_id)] <- 1L

  pt_assigned <- pt[cbind(seq_len(nrow(pt)), pmax(1L, branch_id))]
  pt_assigned <- as.numeric(pt_assigned)

  seu$pseudotime_ref <- pt_assigned
  seu$branch_ref <- paste0("branch_", branch_id)
  seu$ref_lineage_label <- as.character(seu@meta.data[[label_col]])

  # Save reference table for transfer
  tbl <- data.table(
    cell_id = colnames(seu),
    lineage_name = lineage_name,
    ref_label = seu$ref_lineage_label,
    pseudotime_ref = seu$pseudotime_ref,
    branch_ref = seu$branch_ref
  )
  fwrite(tbl, file.path(out_tbl, paste0("09a_", lineage_name, "_reference_pseudotime_table.csv")))

  # UMAP plots
  um <- Embeddings(seu, "umap")
  dfp <- data.table(
    cell_id = rownames(um),
    UMAP_1 = um[, 1],
    UMAP_2 = um[, 2],
    pseudotime_ref = seu$pseudotime_ref,
    branch_ref = seu$branch_ref,
    ref_label = seu$ref_lineage_label
  )

  cent <- dfp[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = ref_label]
  root_xy <- cent[ref_label == root_label]

  p_pt <- ggplot(dfp, aes(UMAP_1, UMAP_2, color = pseudotime_ref)) +
    geom_point(size = 0.25, alpha = 0.8) +
    scale_color_viridis_c(option = "plasma", na.value = "grey85") +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Reference pseudotime: ", lineage_name),
      subtitle = paste0("root = ", root_label),
      color = "Pseudotime"
    )
  if (nrow(root_xy) > 0) {
    p_pt <- p_pt +
      geom_label(
        data = root_xy,
        aes(UMAP_1, UMAP_2, label = paste0("Root: ", ref_label)),
        inherit.aes = FALSE,
        size = 3.2,
        fill = "white",
        alpha = 0.8
      )
  }
  save_pdf(p_pt, file.path(out_plt, paste0("09a_", lineage_name, "_reference_umap_pseudotime")))

  p_branch <- ggplot(dfp, aes(UMAP_1, UMAP_2, color = branch_ref)) +
    geom_point(size = 0.25, alpha = 0.8) +
    theme_classic(base_size = 12) +
    labs(title = paste0("Reference branches: ", lineage_name), color = "Branch")
  save_pdf(p_branch, file.path(out_plt, paste0("09a_", lineage_name, "_reference_umap_branch")))

  p_label <- ggplot(dfp, aes(UMAP_1, UMAP_2, color = ref_label)) +
    geom_point(size = 0.25, alpha = 0.8) +
    theme_classic(base_size = 12) +
    labs(title = paste0("Reference labels: ", lineage_name), color = "Ref label")
  save_pdf(p_label, file.path(out_plt, paste0("09a_", lineage_name, "_reference_umap_labels")))

  # Marker trends in reference
  marker_genes <- unique(spec$marker_genes %||% character(0))
  marker_genes <- marker_genes[marker_genes %in% rownames(seu)]
  if (length(marker_genes) > 0) {
    marker_cells <- dfp$cell_id
    max_marker_cells <- as.integer(spec$marker_trend_max_cells %||% 6000L)
    if (!is.na(max_marker_cells) && max_marker_cells > 0L && length(marker_cells) > max_marker_cells) {
      set.seed(1L)
      marker_cells <- sort(sample(marker_cells, max_marker_cells))
      msg(
        "[", lineage_name, "] marker-trend diagnostic downsampled to ",
        length(marker_cells), " / ", nrow(dfp), " cells."
      )
    }
    df_marker <- dfp[cell_id %in% marker_cells, .(cell_id, pseudotime_ref)]
    expr <- fetch_expr_data(seu[, marker_cells], vars = marker_genes)
    expr <- as.data.table(expr, keep.rownames = "cell_id")
    long <- melt(
      merge(df_marker, expr, by = "cell_id", all.x = TRUE),
      id.vars = c("cell_id", "pseudotime_ref"),
      variable.name = "gene",
      value.name = "expr"
    )
    p_tr <- ggplot(long, aes(x = pseudotime_ref, y = expr)) +
      geom_point(size = 0.1, alpha = 0.2) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#C44E52", linewidth = 0.8) +
      facet_wrap(~ gene, scales = "free_y", ncol = 4) +
      theme_classic(base_size = 11) +
      labs(title = paste0("Reference marker trends: ", lineage_name), x = "Pseudotime", y = "Expression")
    save_pdf(p_tr, file.path(out_plt, paste0("09a_", lineage_name, "_reference_marker_trends")), w = 12, h = 8)
  }

  saveRDS(seu, file.path(out_obj, paste0("09a_", lineage_name, "_reference_lineage.rds")))

  data.table(
    lineage_name = lineage_name,
    n_cells = ncol(seu),
    n_branches = ncol(pt),
    root_label = root_label,
    root_cluster = root_cluster,
    start_cluster = start_level,
    filter_label_col = label_col,
    slingshot_cluster_col = cluster_col,
    pseudotime_min = min(seu$pseudotime_ref, na.rm = TRUE),
    pseudotime_max = max(seu$pseudotime_ref, na.rm = TRUE)
  )
}

msg("Loading references.")
ref_ventral <- readRDS(CFG$paths$ref_ventral_rds)
ref_dorsal <- readRDS(CFG$paths$ref_dorsal_rds)

summ <- vector("list", length(CFG$lineages))
i <- 0L
for (nm in names(CFG$lineages)) {
  i <- i + 1L
  spec <- CFG$lineages[[nm]]
  ref_obj <- if (!is.null(spec$ref_rds)) {
    if (!file.exists(spec$ref_rds)) stop("[", nm, "] missing lineage-specific reference: ", spec$ref_rds)
    readRDS(spec$ref_rds)
  } else if (identical(spec$reference_name, "ventral")) {
    ref_ventral
  } else {
    ref_dorsal
  }
  summ[[i]] <- run_lineage_reference(ref_obj, nm, spec)
}

sum_tbl <- rbindlist(summ, fill = TRUE)
fwrite(sum_tbl, file.path(out_tbl, "09a_reference_lineage_summary.csv"))

msg("09a complete. Outputs under: ", out_refs)
