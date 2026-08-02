#!/usr/bin/env Rscript

# ============================================================
# 09b_map_query_to_reference.R
# Map E18 query subsets to lineage-specific reference trajectories
# and transfer pseudotime/branch labels.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(Matrix)
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

safe_save_pdf <- function(plot_obj, path_no_ext, w = 8.5, h = 6.0, lineage_name = NULL) {
  tryCatch(
    save_pdf(plot_obj, path_no_ext, w = w, h = h),
    error = function(e) {
      msg("[", lineage_name %||% "plot", "] skipped plot ", basename(path_no_ext), ": ", conditionMessage(e))
      NULL
    }
  )
}

save_mapped_query_rds <- function(qry, out_path, lineage_name) {
  ok <- tryCatch({
    gc()
    saveRDS(qry, out_path, compress = FALSE)
    TRUE
  }, error = function(e) {
    msg("[", lineage_name, "] saveRDS failed on full object: ", conditionMessage(e), " ; retrying with DietSeurat.")
    FALSE
  })
  if (ok) return(invisible(TRUE))

  qry_small <- tryCatch({
    assays_keep <- intersect(c("RNA", "SCT"), names(qry@assays))
    if (length(assays_keep) == 0) assays_keep <- names(qry@assays)[1]
    DefaultAssay(qry) <- assays_keep[1]
    dr_keep <- intersect(c("pca", "umap", "ref.umap"), names(qry@reductions))
    DietSeurat(
      object = qry,
      assays = assays_keep,
      dimreducs = dr_keep,
      graphs = NULL
    )
  }, error = function(e) {
    msg("[", lineage_name, "] DietSeurat fallback failed: ", conditionMessage(e), " ; saving minimal assay-only object.")
    assays_keep <- intersect(c("RNA", "SCT"), names(qry@assays))
    if (length(assays_keep) == 0) assays_keep <- names(qry@assays)[1]
    DefaultAssay(qry) <- assays_keep[1]
    qry@reductions <- list()
    qry@graphs <- list()
    qry@commands <- list()
    qry
  })

  gc()
  saveRDS(qry_small, out_path, compress = FALSE)
  invisible(TRUE)
}

stopifnot(file.exists(CFG$paths$query_e18_rds))

out_map <- file.path(CFG$paths$out_dir, "mapping")
out_obj <- file.path(out_map, "objects")
out_tbl <- file.path(out_map, "tables")
out_plt <- file.path(out_map, "plots")
dir.create(out_obj, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_plt, recursive = TRUE, showWarnings = FALSE)

prepare_query_subset <- function(seu_sub, normalization_method = "SCT") {
  if (ncol(seu_sub) < 100) return(seu_sub)
  if (toupper(normalization_method) == "SCT") {
    if (!"SCT" %in% names(seu_sub@assays)) {
      seu_sub <- SCTransform(seu_sub, assay = "RNA", conserve.memory = TRUE, verbose = FALSE)
    }
    DefaultAssay(seu_sub) <- "SCT"
  } else {
    DefaultAssay(seu_sub) <- "RNA"
    seu_sub <- NormalizeData(seu_sub, verbose = FALSE)
    seu_sub <- FindVariableFeatures(seu_sub, verbose = FALSE)
  }
  seu_sub
}

reference_supports_sct <- function(ref_obj) {
  if (!"SCT" %in% names(ref_obj@assays)) return(FALSE)
  sct_assay <- ref_obj[["SCT"]]
  if (inherits(sct_assay, "SCTAssay")) return(TRUE)
  # Seurat v5 objects may store assays differently; this catches legacy slot availability.
  if ("SCTModel.list" %in% slotNames(sct_assay)) return(TRUE)
  FALSE
}

extract_transfer_numeric <- function(pred_obj, n_cells, lineage_name, field_name = "numeric") {
  # Direct numeric vector
  if (is.numeric(pred_obj) && length(pred_obj) == n_cells) {
    return(list(values = as.numeric(pred_obj), score = rep(NA_real_, n_cells)))
  }
  # Common list/data.frame output from label transfer
  if (is.list(pred_obj) && "predicted.id" %in% names(pred_obj)) {
    vals <- suppressWarnings(as.numeric(pred_obj$predicted.id))
    if (length(vals) == n_cells) {
      sc <- if ("prediction.score.max" %in% names(pred_obj)) as.numeric(pred_obj$prediction.score.max) else rep(NA_real_, n_cells)
      return(list(values = vals, score = sc))
    }
  }
  # Sparse matrix / Matrix class from feature transfer
  if (inherits(pred_obj, "Matrix")) {
    rn <- rownames(pred_obj) %||% character()
    use_i <- if ("pseudotime-ref" %in% rn) match("pseudotime-ref", rn) else if ("pseudotime_ref" %in% rn) match("pseudotime_ref", rn) else 1L
    vals <- suppressWarnings(as.numeric(pred_obj[use_i, ]))
    if (length(vals) == n_cells) {
      return(list(values = vals, score = rep(NA_real_, n_cells)))
    }
  }
  # Seurat Assay/Assay5 return from feature transfer.
  if (isS4(pred_obj) && any(class(pred_obj) %in% c("Assay", "Assay5", "SCTAssay"))) {
    mat <- NULL
    # Assay (v4/v5-compat): prefer @data then @counts.
    if ("data" %in% slotNames(pred_obj)) {
      mat <- tryCatch(methods::slot(pred_obj, "data"), error = function(e) NULL)
    }
    if ((is.null(mat) || length(mat) == 0) && "counts" %in% slotNames(pred_obj)) {
      mat <- tryCatch(methods::slot(pred_obj, "counts"), error = function(e) NULL)
    }
    # Assay5: try LayerData.
    if (is.null(mat) || length(mat) == 0) {
      mat <- tryCatch(LayerData(pred_obj, layer = "data"), error = function(e) NULL)
    }
    if ((is.null(mat) || length(mat) == 0)) {
      mat <- tryCatch(LayerData(pred_obj, layer = "counts"), error = function(e) NULL)
    }
    if (!is.null(mat) && length(mat) > 0) {
      rn <- rownames(mat) %||% character()
      use_i <- if ("pseudotime-ref" %in% rn) match("pseudotime-ref", rn) else if ("pseudotime_ref" %in% rn) match("pseudotime_ref", rn) else 1L
      vals <- suppressWarnings(as.numeric(mat[use_i, ]))
      if (length(vals) == n_cells) {
        return(list(values = vals, score = rep(NA_real_, n_cells)))
      }
    }
  }
  # Dense matrix / data.frame
  if (is.matrix(pred_obj) || is.data.frame(pred_obj)) {
    cn <- colnames(pred_obj) %||% character()
    rn <- rownames(pred_obj) %||% character()
    wanted_cols <- c("pseudotime_ref", "predicted.pseudotime_ref", "pseudotime-ref", "predicted.pseudotime-ref")
    wanted_rows <- c("pseudotime_ref", "pseudotime-ref")
    hit_col <- intersect(wanted_cols, cn)
    if (length(hit_col) > 0) {
      vals <- suppressWarnings(as.numeric(pred_obj[, hit_col[1]]))
      if (length(vals) == n_cells) return(list(values = vals, score = rep(NA_real_, n_cells)))
    }
    hit_row <- intersect(wanted_rows, rn)
    if (length(hit_row) > 0) {
      vals <- suppressWarnings(as.numeric(pred_obj[hit_row[1], ]))
      if (length(vals) == n_cells) return(list(values = vals, score = rep(NA_real_, n_cells)))
    }
    if (nrow(pred_obj) == n_cells && ncol(pred_obj) >= 1L) {
      non_score_cols <- setdiff(
        cn[which(vapply(as.data.frame(pred_obj), is.numeric, logical(1)))],
        grep("^prediction\\.score", cn, value = TRUE)
      )
      use_col <- if (length(non_score_cols) > 0) non_score_cols[1] else cn[1]
      vals <- suppressWarnings(as.numeric(pred_obj[, use_col]))
      if (length(vals) == n_cells) return(list(values = vals, score = rep(NA_real_, n_cells)))
    }
    if (ncol(pred_obj) == n_cells && nrow(pred_obj) >= 1L) {
      vals <- suppressWarnings(as.numeric(pred_obj[1, ]))
      if (length(vals) == n_cells) return(list(values = vals, score = rep(NA_real_, n_cells)))
    }
  }
  stop("[", lineage_name, "] could not parse ", field_name, " transfer output. class=", paste(class(pred_obj), collapse = "/"))
}

extract_transfer_label <- function(pred_obj, n_cells, lineage_name, field_name = "label") {
  # Standard TransferData label output
  if (is.list(pred_obj) && "predicted.id" %in% names(pred_obj)) {
    ids <- as.character(pred_obj$predicted.id)
    if (length(ids) == n_cells) {
      sc <- if ("prediction.score.max" %in% names(pred_obj)) as.numeric(pred_obj$prediction.score.max) else rep(NA_real_, n_cells)
      return(list(ids = ids, score = sc))
    }
  }
  if (is.data.frame(pred_obj) || is.matrix(pred_obj)) {
    cn <- colnames(pred_obj) %||% character()
    if ("predicted.id" %in% cn) {
      ids <- as.character(pred_obj[, "predicted.id"])
      sc <- if ("prediction.score.max" %in% cn) as.numeric(pred_obj[, "prediction.score.max"]) else rep(NA_real_, length(ids))
      if (length(ids) == n_cells) return(list(ids = ids, score = sc))
    }
    # If only one column, treat it as predicted id.
    if (nrow(pred_obj) == n_cells && ncol(pred_obj) >= 1L) {
      ids <- as.character(pred_obj[, 1])
      if (length(ids) == n_cells) return(list(ids = ids, score = rep(NA_real_, n_cells)))
    }
  }
  stop("[", lineage_name, "] could not parse ", field_name, " transfer output. class=", paste(class(pred_obj), collapse = "/"))
}

compute_pt_confidence_from_weights <- function(w_obj, ref_pt, query_cells, lineage_name) {
  if (!inherits(w_obj, "Matrix")) {
    stop("[", lineage_name, "] weight object is not Matrix: ", paste(class(w_obj), collapse = "/"))
  }
  w <- w_obj
  qn <- query_cells
  rn <- rownames(w) %||% character()
  cn <- colnames(w) %||% character()

  if (length(rn) > 0 && sum(qn %in% rn) > 0) {
    # rows are query cells
  } else if (length(cn) > 0 && sum(qn %in% cn) > 0) {
    w <- Matrix::t(w)
    rn <- rownames(w) %||% character()
    cn <- colnames(w) %||% character()
  } else if (nrow(w) != length(qn) && ncol(w) == length(qn)) {
    w <- Matrix::t(w)
    rn <- rownames(w) %||% character()
    cn <- colnames(w) %||% character()
  }

  # Reorder rows to query cell order when possible.
  if (length(rn) > 0 && all(qn %in% rn)) {
    w <- w[qn, , drop = FALSE]
  } else if (nrow(w) != length(qn)) {
    stop("[", lineage_name, "] weights matrix rows do not match query cells.")
  }

  # Align reference pseudotime to weight columns.
  ref_cells <- colnames(w) %||% character()
  if (length(ref_cells) > 0 && !is.null(names(ref_pt)) && all(ref_cells %in% names(ref_pt))) {
    ref_vec <- as.numeric(ref_pt[ref_cells])
  } else if (ncol(w) == length(ref_pt)) {
    ref_vec <- as.numeric(ref_pt)
  } else {
    stop("[", lineage_name, "] cannot align weight columns to reference pseudotime.")
  }

  n <- nrow(w)
  score <- rep(NA_real_, n)
  maxw <- rep(NA_real_, n)
  wsd <- rep(NA_real_, n)
  sdat <- Matrix::summary(w)
  if (nrow(sdat) > 0) {
    dtw <- as.data.table(sdat)
    # Matrix::summary names: i, j, x (1-based indices)
    dtw[, pt := ref_vec[j]]
    dtw[, wnorm := x / sum(x), by = i]
    dtw[, m := sum(wnorm * pt), by = i]
    res <- dtw[, .(
      wsd = sqrt(sum(wnorm * (pt - m[1])^2)),
      maxw = max(wnorm)
    ), by = i]
    idx <- as.integer(res$i)
    wsd[idx] <- res$wsd
    maxw[idx] <- res$maxw
    score[idx] <- (1 / (1 + res$wsd)) * res$maxw
  }
  list(score = score, maxw = maxw, wsd = wsd)
}

ensure_sct_reference <- function(ref_obj, lineage_name, spec, dims_use) {
  need_sct <- isTRUE(CFG$mapping$force_sct_reference_model) || toupper(spec$normalization_method %||% CFG$mapping$normalization_method) == "SCT"
  if (!need_sct) return(ref_obj)
  if (!reference_supports_sct(ref_obj)) {
    msg("[", lineage_name, "] building SCT model on reference lineage object.")
    nfeat <- as.integer(spec$anchor_features %||% CFG$mapping$anchor_features %||% 1500L)
    nfeat <- max(500L, nfeat)
    ref_obj <- SCTransform(
      ref_obj,
      assay = "RNA",
      new.assay.name = "SCT",
      variable.features.n = nfeat,
      conserve.memory = TRUE,
      verbose = FALSE
    )
  }
  DefaultAssay(ref_obj) <- "SCT"
  # Ensure PCA exists for the requested dims in SCT space.
  need_pca <- TRUE
  if ("pca" %in% names(ref_obj@reductions)) {
    need_pca <- ncol(Embeddings(ref_obj, "pca")) < max(dims_use)
  }
  if (need_pca) {
    ref_obj <- RunPCA(ref_obj, npcs = max(dims_use), verbose = FALSE)
  }
  ref_obj
}

run_map_lineage <- function(query_obj, lineage_name, spec) {
  ref_path <- file.path(CFG$paths$out_dir, "refs", "objects", paste0("09a_", lineage_name, "_reference_lineage.rds"))
  if (!file.exists(ref_path)) {
    stop("[", lineage_name, "] missing reference lineage object from 09a: ", ref_path)
  }
  ref <- readRDS(ref_path)
  run_meta <- data.table(
    lineage_name = lineage_name,
    filter_label_col = NA_character_,
    slingshot_cluster_col = NA_character_,
    root_label = NA_character_,
    start_cluster = NA_character_
  )
  sum_path <- file.path(CFG$paths$out_dir, "refs", "tables", "09a_reference_lineage_summary.csv")
  if (file.exists(sum_path)) {
    sm <- tryCatch(fread(sum_path), error = function(e) NULL)
    if (!is.null(sm) && "lineage_name" %in% colnames(sm)) {
      lineage_key <- lineage_name
      hit <- sm[sm$lineage_name == lineage_key, ]
      if (nrow(hit) > 0) {
        run_meta <- as.data.table(hit[1])
      }
    }
  }
  if (!all(c("pseudotime_ref", "branch_ref") %in% colnames(ref@meta.data))) {
    # Fallback: recover from 09a table if metadata columns are missing in object.
    tab_path <- file.path(CFG$paths$out_dir, "refs", "tables", paste0("09a_", lineage_name, "_reference_pseudotime_table.csv"))
    if (!file.exists(tab_path)) {
      stop("[", lineage_name, "] reference lineage object lacks pseudotime_ref/branch_ref and fallback table not found.")
    }
    tab <- fread(tab_path)
    if (!all(c("cell_id", "pseudotime_ref", "branch_ref") %in% colnames(tab))) {
      stop("[", lineage_name, "] fallback table missing required columns.")
    }
    idx <- match(colnames(ref), tab$cell_id)
    ref$pseudotime_ref <- tab$pseudotime_ref[idx]
    ref$branch_ref <- tab$branch_ref[idx]
    if ("ref_label" %in% colnames(tab)) {
      ref$ref_lineage_label <- tab$ref_label[idx]
    }
  }
  ref_cells <- colnames(ref)
  ref_pt <- ref@meta.data[ref_cells, "pseudotime_ref", drop = TRUE]
  ref_branch <- ref@meta.data[ref_cells, "branch_ref", drop = TRUE]
  ref_label <- if ("ref_lineage_label" %in% colnames(ref@meta.data)) {
    ref@meta.data[ref_cells, "ref_lineage_label", drop = TRUE]
  } else {
    ref@meta.data[ref_cells, 1, drop = TRUE]
  }
  names(ref_pt) <- ref_cells
  names(ref_branch) <- ref_cells
  names(ref_label) <- ref_cells
  # Ensure transfer vectors are valid (named and non-missing).
  if (all(is.na(ref_pt))) stop("[", lineage_name, "] reference pseudotime_ref is all NA.")
  if (anyNA(ref_pt)) ref_pt[is.na(ref_pt)] <- median(ref_pt, na.rm = TRUE)
  ref_branch <- as.character(ref_branch)
  names(ref_branch) <- ref_cells
  ref_branch[is.na(ref_branch)] <- "branch_1"
  ref_label <- as.character(ref_label)
  names(ref_label) <- ref_cells
  ref_label[is.na(ref_label)] <- "unknown"
  # Seurat v5 TransferData treats numeric refdata as feature-transfer input (matrix),
  # not as label-transfer input (vector). Keep pseudotime as 1xN matrix.
  ref_pt_mat <- matrix(
    as.numeric(ref_pt),
    nrow = 1L,
    dimnames = list("pseudotime-ref", ref_cells)
  )
  ref_branch_fac <- factor(ref_branch)
  names(ref_branch_fac) <- ref_cells
  ref_label_fac <- factor(ref_label)
  names(ref_label_fac) <- ref_cells

  # Force query broad labels from celltype_broad (do not use celltype_label).
  broad_col <- "celltype_broad"
  if (!broad_col %in% colnames(query_obj@meta.data)) {
    stop("Query broad_col missing: ", broad_col)
  }

  keep_broad <- spec$query_keep_broad
  q_cells <- rownames(query_obj@meta.data)[query_obj@meta.data[[broad_col]] %in% keep_broad]
  if (length(q_cells) < 100) {
    msg("[", lineage_name, "] skip: too few query cells after broad filter (", length(q_cells), ").")
    return(data.table(lineage_name = lineage_name, n_query_cells = length(q_cells), mapped = FALSE))
  }

  qry <- subset(query_obj, cells = q_cells)

  requested_norm <- toupper(spec$normalization_method %||% CFG$mapping$normalization_method)
  use_norm <- if (requested_norm == "SCT") "SCT" else "LogNormalize"

  if (use_norm != "SCT") {
    DefaultAssay(ref) <- "RNA"
    ref <- NormalizeData(ref, verbose = FALSE)
    ref <- FindVariableFeatures(ref, verbose = FALSE)
  }

  qry <- prepare_query_subset(qry, normalization_method = use_norm)

  dims_cfg <- spec$dims %||% CFG$mapping$dims
  dims_cfg <- as.integer(dims_cfg)
  dims_cfg <- dims_cfg[is.finite(dims_cfg) & dims_cfg > 0]
  if (length(dims_cfg) == 0) dims_cfg <- 1:20
  dims_use <- seq_len(max(dims_cfg))
  k_weight_use <- as.integer(spec$k_weight %||% CFG$mapping$k_weight)
  if (!is.finite(k_weight_use) || k_weight_use < 5) k_weight_use <- 30L

  # Build an SCT model on reference lineage objects if needed.
  if (use_norm == "SCT") {
    ref <- ensure_sct_reference(ref, lineage_name, spec, dims_use = dims_use)
  } else {
    if (!"pca" %in% names(ref@reductions)) {
      ref <- RunPCA(ref, npcs = max(dims_use), verbose = FALSE)
    }
    if (ncol(Embeddings(ref, "pca")) < max(dims_use)) {
      ref <- RunPCA(ref, npcs = max(dims_use), verbose = FALSE)
    }
  }
  if (!"pca" %in% names(ref@reductions)) {
    stop("[", lineage_name, "] reference PCA reduction missing after preparation.")
  }
  dims_use <- seq_len(min(max(dims_cfg), ncol(Embeddings(ref, "pca"))))

  # Use explicit shared features to avoid "No features to use" failures.
  anchor_n <- as.integer(spec$anchor_features %||% CFG$mapping$anchor_features %||% 1500L)
  anchor_n <- max(500L, anchor_n)
  ref_vf <- VariableFeatures(ref)
  qry_vf <- VariableFeatures(qry)
  shared_features <- intersect(ref_vf, qry_vf)
  if (length(shared_features) < 200) {
    shared_features <- intersect(rownames(ref), rownames(qry))
  }
  if (length(shared_features) == 0) {
    stop("[", lineage_name, "] no shared features between reference and query for anchor finding.")
  }
  features_use <- head(shared_features, min(anchor_n, length(shared_features)))

  msg("[", lineage_name, "] FindTransferAnchors (", use_norm, ", dims=", max(dims_use), ", features=", length(features_use), ").")
  anchors <- tryCatch(
    FindTransferAnchors(
      reference = ref,
      query = qry,
      normalization.method = use_norm,
      reference.reduction = "pca",
      dims = dims_use,
      features = features_use
    ),
    error = function(e) {
      if (identical(use_norm, "SCT") && isTRUE(CFG$mapping$allow_lognorm_fallback)) {
        msg("[", lineage_name, "] SCT anchors failed (", conditionMessage(e), "); retrying with LogNormalize because allow_lognorm_fallback=TRUE.")
        use_norm <<- "LogNormalize"
        DefaultAssay(ref) <- "RNA"
        ref <<- NormalizeData(ref, verbose = FALSE)
        ref <<- FindVariableFeatures(ref, nfeatures = 2000, verbose = FALSE)
        if (!"pca" %in% names(ref@reductions) || ncol(Embeddings(ref, "pca")) < max(dims_use)) {
          ref <<- RunPCA(ref, npcs = max(dims_use), verbose = FALSE)
        }
        DefaultAssay(qry) <- "RNA"
        qry <<- NormalizeData(qry, verbose = FALSE)
        qry <<- FindVariableFeatures(qry, nfeatures = 2000, verbose = FALSE)
        shared2 <- intersect(VariableFeatures(ref), VariableFeatures(qry))
        if (length(shared2) < 200) shared2 <- intersect(rownames(ref), rownames(qry))
        if (length(shared2) == 0) stop("[", lineage_name, "] no shared features after LogNormalize fallback.")
        features2 <- head(shared2, min(anchor_n, length(shared2)))
        return(FindTransferAnchors(
          reference = ref,
          query = qry,
          normalization.method = use_norm,
          reference.reduction = "pca",
          dims = dims_use,
          features = features2
        ))
      }
      stop("[", lineage_name, "] FindTransferAnchors failed: ", conditionMessage(e))
    }
  )

  msg("[", lineage_name, "] Transfer pseudotime/branch/ref_label.")
  pred_pt <- tryCatch(
    TransferData(
      anchorset = anchors,
      refdata = ref_pt_mat,
      dims = dims_use,
      k.weight = k_weight_use
    ),
    error = function(e) stop("[", lineage_name, "] TransferData pseudotime failed: ", conditionMessage(e))
  )
  pred_branch <- tryCatch(
    TransferData(anchorset = anchors, refdata = ref_branch_fac, dims = dims_use, k.weight = k_weight_use),
    error = function(e) stop("[", lineage_name, "] TransferData branch failed: ", conditionMessage(e))
  )
  pred_label <- tryCatch(
    TransferData(anchorset = anchors, refdata = ref_label_fac, dims = dims_use, k.weight = k_weight_use),
    error = function(e) stop("[", lineage_name, "] TransferData label failed: ", conditionMessage(e))
  )

  pt_parsed <- extract_transfer_numeric(pred_pt, n_cells = ncol(qry), lineage_name = lineage_name, field_name = "pseudotime")
  br_parsed <- extract_transfer_label(pred_branch, n_cells = ncol(qry), lineage_name = lineage_name, field_name = "branch")
  lb_parsed <- extract_transfer_label(pred_label, n_cells = ncol(qry), lineage_name = lineage_name, field_name = "ref_label")
  qry$pseudotime_transfer <- pt_parsed$values
  # Build a pseudotime confidence from transfer weights (anchor/weight based),
  # because feature transfer does not reliably provide prediction.score.max.
  pt_conf <- tryCatch(
    {
      w <- TransferData(
        anchorset = anchors,
        refdata = ref_pt_mat,
        dims = dims_use,
        k.weight = k_weight_use,
        only.weights = TRUE
      )
      compute_pt_confidence_from_weights(
        w_obj = w,
        ref_pt = ref_pt,
        query_cells = colnames(qry),
        lineage_name = lineage_name
      )
    },
    error = function(e) {
      msg("[", lineage_name, "] pseudotime confidence fallback (no weights): ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(pt_conf)) {
    qry$pseudotime_transfer_score <- pt_conf$score
    qry$pseudotime_transfer_maxweight <- pt_conf$maxw
    qry$pseudotime_transfer_wsd <- pt_conf$wsd
  } else {
    qry$pseudotime_transfer_score <- pt_parsed$score
  }
  qry$branch_transfer <- br_parsed$ids
  qry$branch_transfer_score <- br_parsed$score
  qry$ref_label_transfer <- lb_parsed$ids
  qry$ref_label_transfer_score <- lb_parsed$score
  used_manual <- FALSE

  if ("filter_label_col" %in% colnames(run_meta)) qry$ref_filter_label_col <- as.character(run_meta$filter_label_col[1])
  if ("slingshot_cluster_col" %in% colnames(run_meta)) qry$ref_slingshot_cluster_col <- as.character(run_meta$slingshot_cluster_col[1])
  if ("root_label" %in% colnames(run_meta)) qry$ref_root_label <- as.character(run_meta$root_label[1])
  if ("start_cluster" %in% colnames(run_meta)) qry$ref_start_cluster <- as.character(run_meta$start_cluster[1])

  # Optional projection to reference UMAP model
  projection_ok <- FALSE
  do_projection <- isTRUE(CFG$mapping$run_mapquery_projection) && !used_manual
  if (do_projection && "umap" %in% names(ref@reductions)) {
    mq <- tryCatch(
      suppressWarnings(MapQuery(
        anchorset = anchors,
        query = qry,
        reference = ref,
        refdata = NULL,
        reference.reduction = "pca",
        reduction.model = "umap",
        verbose = FALSE
      )),
      error = function(e) {
        msg("[", lineage_name, "] MapQuery skipped: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(mq)) {
      qry <- mq
      projection_ok <- TRUE
    }
  } else if (isTRUE(CFG$mapping$run_mapquery_projection) && used_manual) {
    msg("[", lineage_name, "] MapQuery projection skipped because manual transfer was used.")
  }

  # Save mapped object + table
  save_mapped_query_rds(
    qry = qry,
    out_path = file.path(out_obj, paste0("09b_", lineage_name, "_mapped_query.rds")),
    lineage_name = lineage_name
  )

  keep_cols <- unique(c(
    "cell_id",
    CFG$cols$genotype_col,
    CFG$cols$broad_col,
    CFG$cols$query_sample_col,
    "pseudotime_transfer",
    "pseudotime_transfer_score",
    "pseudotime_transfer_maxweight",
    "pseudotime_transfer_wsd",
    "branch_transfer",
    "branch_transfer_score",
    "ref_label_transfer",
    "ref_label_transfer_score"
  ))
  keep_cols <- keep_cols[keep_cols == "cell_id" | keep_cols %in% colnames(qry@meta.data)]
  tbl <- as.data.table(qry@meta.data, keep.rownames = "cell_id")[, ..keep_cols]
  tbl[, lineage_name := lineage_name]
  if ("filter_label_col" %in% colnames(run_meta)) tbl[, ref_filter_label_col := as.character(run_meta$filter_label_col[1])]
  if ("slingshot_cluster_col" %in% colnames(run_meta)) tbl[, ref_slingshot_cluster_col := as.character(run_meta$slingshot_cluster_col[1])]
  if ("root_label" %in% colnames(run_meta)) tbl[, ref_root_label := as.character(run_meta$root_label[1])]
  if ("start_cluster" %in% colnames(run_meta)) tbl[, ref_start_cluster := as.character(run_meta$start_cluster[1])]
  fwrite(tbl, file.path(out_tbl, paste0("09b_", lineage_name, "_mapped_table.csv")))

  # Plot transferred pseudotime
  red_name <- if ("ref.umap" %in% names(qry@reductions)) "ref.umap" else if ("umap" %in% names(qry@reductions)) "umap" else NA_character_
  if (!is.na(red_name)) {
    em <- Embeddings(qry, red_name)
    dplot <- data.table(
      UMAP_1 = em[, 1],
      UMAP_2 = em[, 2],
      pseudotime_transfer = qry$pseudotime_transfer,
      branch_transfer = qry$branch_transfer
    )
    dplot <- dplot[is.finite(UMAP_1) & is.finite(UMAP_2)]
    dplot[, pseudotime_transfer := suppressWarnings(as.numeric(pseudotime_transfer))]
    dplot[!is.finite(pseudotime_transfer), pseudotime_transfer := NA_real_]
    dplot[, branch_transfer := as.character(branch_transfer)]
    dplot[is.na(branch_transfer) | branch_transfer == "", branch_transfer := "unknown"]
    dplot[, branch_transfer := factor(branch_transfer)]
    if (CFG$cols$genotype_col %in% colnames(qry@meta.data)) {
      dplot[, genotype := as.character(qry@meta.data[[CFG$cols$genotype_col]])]
    } else {
      dplot[, genotype := "NA"]
    }

    if (nrow(dplot) > 0) {
      p1 <- ggplot(dplot, aes(UMAP_1, UMAP_2, color = pseudotime_transfer)) +
        geom_point(size = 0.25, alpha = 0.8) +
        scale_color_viridis_c(option = "plasma", na.value = "grey85") +
        theme_classic(base_size = 12) +
        labs(title = paste0("Mapped query pseudotime: ", lineage_name), color = "Transferred PT")
      safe_save_pdf(
        p1,
        file.path(out_plt, paste0("09b_", lineage_name, "_mapped_umap_pseudotime")),
        lineage_name = lineage_name
      )

      p2 <- ggplot(dplot, aes(UMAP_1, UMAP_2, color = branch_transfer)) +
        geom_point(size = 0.25, alpha = 0.8) +
        theme_classic(base_size = 12) +
        labs(title = paste0("Mapped query branches: ", lineage_name), color = "Transferred branch")
      safe_save_pdf(
        p2,
        file.path(out_plt, paste0("09b_", lineage_name, "_mapped_umap_branch")),
        lineage_name = lineage_name
      )
    }

    if (length(unique(dplot$genotype)) > 1) {
      p3 <- ggplot(dplot, aes(UMAP_1, UMAP_2, color = pseudotime_transfer)) +
        geom_point(size = 0.2, alpha = 0.8) +
        scale_color_viridis_c(option = "plasma", na.value = "grey85") +
        facet_wrap(~ genotype) +
        theme_classic(base_size = 12) +
        labs(title = paste0("Mapped query pseudotime split by genotype: ", lineage_name))
      safe_save_pdf(
        p3,
        file.path(out_plt, paste0("09b_", lineage_name, "_mapped_umap_pseudotime_split_genotype")),
        w = 11, h = 5, lineage_name = lineage_name
      )
    }
  }

  med_pt <- if (all(is.na(qry$pseudotime_transfer_score))) NA_real_ else median(qry$pseudotime_transfer_score, na.rm = TRUE)
  med_branch <- if (all(is.na(qry$branch_transfer_score))) NA_real_ else median(qry$branch_transfer_score, na.rm = TRUE)
  data.table(
    lineage_name = lineage_name,
    n_query_cells = ncol(qry),
    mapped = TRUE,
    used_manual_transfer = used_manual,
    projection_ok = projection_ok,
    ref_filter_label_col = if ("filter_label_col" %in% colnames(run_meta)) as.character(run_meta$filter_label_col[1]) else NA_character_,
    ref_slingshot_cluster_col = if ("slingshot_cluster_col" %in% colnames(run_meta)) as.character(run_meta$slingshot_cluster_col[1]) else NA_character_,
    median_transfer_score_pt = med_pt,
    median_transfer_score_branch = med_branch
  )
}

msg("Loading query object.")
query <- readRDS(CFG$paths$query_e18_rds)

sum_rows <- vector("list", length(CFG$lineages))
i <- 0L
for (nm in names(CFG$lineages)) {
  i <- i + 1L
  sum_rows[[i]] <- run_map_lineage(query, nm, CFG$lineages[[nm]])
}

sum_tbl <- rbindlist(sum_rows, fill = TRUE)
fwrite(sum_tbl, file.path(out_tbl, "09b_mapping_summary.csv"))

msg("09b complete. Outputs under: ", out_map)
