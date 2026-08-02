#!/usr/bin/env Rscript

# ============================================================
# 08_build_dorsal_reference_dibella.R
# Build dorsal embryonic reference from Di Bella et al. processed matrix
# Files:
#   - gene_sorted-matrix.mtx.gz
#   - genes.tsv
#   - barcodes.tsv
#   - metaData_scDevSC.txt
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(ggplot2)
})

# ------------------------------
# Parameters
# ------------------------------
# GSE153164 (Di Bella et al. 2021). Set REF_DIBELLA to the download directory.
base_dir <- Sys.getenv("REF_DIBELLA", unset = "data/references/DiBella_GSE153164")
mtx_file <- file.path(base_dir, "gene_sorted-matrix.mtx.gz")
genes_file <- file.path(base_dir, "genes.tsv")
barcodes_file <- file.path(base_dir, "barcodes.tsv")
metadata_file <- file.path(base_dir, "metaData_scDevSC.txt")

out_dir <- file.path(ROOT, "results/dibella_dorsal_reference")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
plots_dir <- file.path(out_dir, "plots")
tables_dir <- file.path(out_dir, "tables")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Filters
keep_embryonic_only <- TRUE      # drop P1/P4
keep_biosamples <- NULL          # e.g. c("E12","E13","E14","E15","E16","E17","E18_S1","E18_S3")
keep_gral_celltypes <- c("Apical progenitors", "Intermediate progenitors", "Excitatory neurons", "Glia")
exclude_null_gral <- TRUE
exclude_scrublet_doublets <- TRUE

# QC
min_features <- 200L
max_percent_mt <- 20
use_ncount_cap <- TRUE
max_ncount <- 50000L

# Batch correction / integration
do_sct_integration <- TRUE
batch_correction_method <- "harmony" # "harmony" (memory-safe) or "sct_rpca"
integration_batch_col <- "biosample_id"  # "biosample_id" or "sample_id" (if present)
sct_integration_nfeatures <- 3000L
sct_anchor_dims <- 1:50
sct_umap_dims <- 1:30
sct_ncells <- 5000L
force_sequential_future <- TRUE
future_globals_maxsize <- Inf

# Non-integrated fallback dimensionality reduction
n_variable_features <- 3000L
npcs <- 50L
umap_dims <- 1:30
pt_size <- 0.25
plot_w <- 9
plot_h <- 7

# Outputs
out_rds <- file.path(out_dir, "DiBella_dorsal_reference_embryonic_filtered.rds")
out_rds_integrated <- file.path(out_dir, "DiBella_dorsal_reference_SCTintegrated.rds")
out_cells_before <- file.path(tables_dir, "cells_before_filtering.csv")
out_cells_after_meta <- file.path(tables_dir, "cells_after_metadata_filters.csv")
out_qc_sample <- file.path(tables_dir, "qc_before_after_by_biosample.csv")
out_qc_stage <- file.path(tables_dir, "qc_before_after_by_stage.csv")
out_kept_cell_ids <- file.path(tables_dir, "kept_cell_ids.txt")

out_umap_stage <- file.path(plots_dir, "08_UMAP_by_stage.pdf")
out_umap_biosample <- file.path(plots_dir, "08_UMAP_by_biosample.pdf")
out_umap_gral <- file.path(plots_dir, "08_UMAP_by_Gral_cellType.pdf")
out_umap_broad <- file.path(plots_dir, "08_UMAP_by_celltype_broad_dorsal.pdf")
out_umap_biosample_split_stage <- file.path(plots_dir, "08_UMAP_by_biosample_split_by_stage.pdf")

# ------------------------------
# Helpers
# ------------------------------
msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

parse_stage <- function(biosample_id) {
  # E10, E18_S1, E18_S3 -> E10.5, E18.5
  x <- sub("_S[0-9]+$", "", biosample_id)
  n <- suppressWarnings(as.integer(sub("^E", "", x)))
  out <- ifelse(!is.na(n), sprintf("E%d.5", n), NA_character_)
  out
}

assert_exists <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
}

# ------------------------------
# Validate inputs
# ------------------------------
assert_exists(mtx_file)
assert_exists(genes_file)
assert_exists(barcodes_file)
assert_exists(metadata_file)

# ------------------------------
# Read metadata and filter target cells
# ------------------------------
msg("Reading metadata: ", metadata_file)
meta <- fread(metadata_file, na.strings = c("NA", ""))

# Remove metadata type row
if ("NAME" %in% names(meta)) {
  meta <- meta[NAME != "TYPE"]
}

# Parse fields
meta[, biosample_id := as.character(biosample_id)]
meta[, stage := parse_stage(biosample_id)]
meta[, sample_id := if ("donor_id" %in% names(meta)) as.character(donor_id) else as.character(biosample_id)]
meta[, dataset := "DiBella_Nature2021"]
meta[, modality := "scRNA"]
meta[, scrublet_doublet := tolower(as.character(scrublet_doublet)) %in% c("true", "t", "1")]

before_counts <- meta[, .(n_cells = .N), by = .(biosample_id, stage, Gral_cellType)][order(stage, biosample_id, Gral_cellType)]
fwrite(before_counts, out_cells_before)

msg("Applying metadata filters.")
meta_keep <- copy(meta)

if (isTRUE(keep_embryonic_only)) {
  meta_keep <- meta_keep[grepl("^E", biosample_id)]
}

if (!is.null(keep_biosamples)) {
  meta_keep <- meta_keep[biosample_id %in% keep_biosamples]
}

if (isTRUE(exclude_null_gral)) {
  meta_keep <- meta_keep[!is.na(Gral_cellType) & Gral_cellType != "Null"]
}

if (!is.null(keep_gral_celltypes) && length(keep_gral_celltypes) > 0) {
  meta_keep <- meta_keep[Gral_cellType %in% keep_gral_celltypes]
}

if (isTRUE(exclude_scrublet_doublets) && "scrublet_doublet" %in% names(meta_keep)) {
  meta_keep <- meta_keep[is.na(scrublet_doublet) | scrublet_doublet == FALSE]
}

if (nrow(meta_keep) == 0) {
  stop("No cells left after metadata filtering. Check keep_biosamples / keep_gral_celltypes.")
}

after_meta_counts <- meta_keep[, .(n_cells = .N), by = .(biosample_id, stage, Gral_cellType)][order(stage, biosample_id, Gral_cellType)]
fwrite(after_meta_counts, out_cells_after_meta)

target_cells <- unique(meta_keep$NAME)
msg("Cells retained by metadata filters: ", length(target_cells), " / ", nrow(meta))

# ------------------------------
# Read genes + barcodes
# ------------------------------
msg("Reading genes/barcodes.")
genes_dt <- fread(genes_file, header = FALSE)
barcodes_dt <- fread(barcodes_file, header = FALSE)

gene_symbols <- if (ncol(genes_dt) >= 2) genes_dt[[2]] else genes_dt[[1]]
gene_symbols <- make.unique(as.character(gene_symbols))
barcodes <- as.character(barcodes_dt[[1]])

# Keep only selected columns from the matrix
keep_idx <- which(barcodes %in% target_cells)
if (length(keep_idx) == 0) {
  stop("No overlap between metadata NAME and barcodes.tsv.")
}
kept_barcodes <- barcodes[keep_idx]

meta_keep <- meta_keep[match(kept_barcodes, NAME)]
if (any(is.na(meta_keep$NAME))) {
  stop("Failed to align metadata rows to kept barcodes.")
}
writeLines(kept_barcodes, out_kept_cell_ids)

# ------------------------------
# Read matrix and subset to target cells
# ------------------------------
msg("Reading sparse matrix (this can take a while): ", mtx_file)
mat <- tryCatch(
  readMM(mtx_file),
  error = function(e1) {
    msg("Direct readMM(path) failed, retrying via gzfile connection.")
    con <- gzfile(mtx_file, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    tryCatch(
      readMM(con),
      error = function(e2) {
        stop(
          "Failed to read matrix market file.\n",
          "Path mode error: ", conditionMessage(e1), "\n",
          "Connection mode error: ", conditionMessage(e2)
        )
      }
    )
  }
)

if (nrow(mat) != length(gene_symbols)) {
  stop("Gene dimension mismatch: matrix rows=", nrow(mat), " genes.tsv rows=", length(gene_symbols))
}
if (ncol(mat) != length(barcodes)) {
  stop("Cell dimension mismatch: matrix cols=", ncol(mat), " barcodes.tsv rows=", length(barcodes))
}

rownames(mat) <- gene_symbols
colnames(mat) <- barcodes

msg("Subsetting matrix to filtered cells.")
mat <- mat[, keep_idx, drop = FALSE]
mat <- as(mat, "dgCMatrix")
gc()

# ------------------------------
# Build Seurat object + metadata
# ------------------------------
msg("Building Seurat object.")
seu <- CreateSeuratObject(
  counts = mat,
  project = "DiBella_dorsal_reference",
  min.cells = 0,
  min.features = 0
)

meta_df <- as.data.frame(meta_keep)
rownames(meta_df) <- meta_df$NAME
common_cols <- setdiff(colnames(meta_df), colnames(seu@meta.data))
seu <- AddMetaData(seu, metadata = meta_df[colnames(seu), common_cols, drop = FALSE])

# Broad dorsal label from provided annotation
seu$celltype_broad_dorsal <- "Other"
seu$celltype_broad_dorsal[seu$Gral_cellType == "Apical progenitors"] <- "RG"
seu$celltype_broad_dorsal[seu$Gral_cellType == "Intermediate progenitors"] <- "IPC"
seu$celltype_broad_dorsal[seu$Gral_cellType == "Excitatory neurons"] <- "Excitatory"
seu$celltype_broad_dorsal[seu$Gral_cellType == "Glia"] <- "Glia"

# ------------------------------
# QC on filtered set
# ------------------------------
msg("Computing QC metrics.")
mt_genes <- grep("^mt-", rownames(seu), value = TRUE)
if (length(mt_genes) == 0) {
  mt_genes <- grep("^MT-", rownames(seu), value = TRUE)
}
if (length(mt_genes) > 0) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mt_genes)
} else {
  seu$percent.mt <- 0
}

qc_meta <- as.data.table(seu@meta.data, keep.rownames = "cell_id")
qc_meta[, pass_min_features := nFeature_RNA >= min_features]
qc_meta[, pass_mt := percent.mt <= max_percent_mt]
if (isTRUE(use_ncount_cap)) {
  qc_meta[, pass_ncount_cap := nCount_RNA <= max_ncount]
} else {
  qc_meta[, pass_ncount_cap := TRUE]
}
qc_meta[, pass_qc := pass_min_features & pass_mt & pass_ncount_cap]

qc_before_sample <- qc_meta[, .(n_before = .N), by = .(biosample_id, stage)]
qc_after_sample <- qc_meta[pass_qc == TRUE, .(n_after = .N), by = .(biosample_id, stage)]
qc_sum_sample <- merge(qc_before_sample, qc_after_sample, by = c("biosample_id", "stage"), all.x = TRUE)
qc_sum_sample[is.na(n_after), n_after := 0L]
qc_sum_sample[, frac_retained := n_after / n_before]
setorder(qc_sum_sample, stage, biosample_id)
fwrite(qc_sum_sample, out_qc_sample)

qc_before_stage <- qc_meta[, .(n_before = .N), by = .(stage)]
qc_after_stage <- qc_meta[pass_qc == TRUE, .(n_after = .N), by = .(stage)]
qc_sum_stage <- merge(qc_before_stage, qc_after_stage, by = "stage", all.x = TRUE)
qc_sum_stage[is.na(n_after), n_after := 0L]
qc_sum_stage[, frac_retained := n_after / n_before]
setorder(qc_sum_stage, stage)
fwrite(qc_sum_stage, out_qc_stage)

cells_keep_qc <- qc_meta[pass_qc == TRUE, cell_id]
msg("Cells after QC: ", length(cells_keep_qc), " / ", nrow(qc_meta))
if (length(cells_keep_qc) < 1000) {
  warning("Low number of cells after QC: ", length(cells_keep_qc))
}
seu <- subset(seu, cells = cells_keep_qc)

# ------------------------------
# Batch correction + dimensional reduction
# ------------------------------
if (isTRUE(do_sct_integration)) {
  if (requireNamespace("future", quietly = TRUE)) {
    if (is.infinite(future_globals_maxsize)) {
      options(future.globals.maxSize = +Inf)
    } else {
      options(future.globals.maxSize = as.numeric(future_globals_maxsize))
    }
    if (isTRUE(force_sequential_future)) {
      future::plan("sequential")
    }
  }

  if (!integration_batch_col %in% colnames(seu@meta.data)) {
    stop("integration_batch_col not found in metadata: ", integration_batch_col)
  }

  msg("Running batch correction (", batch_correction_method, ") by ", integration_batch_col, ".")

  cc <- Seurat::cc.genes.updated.2019

  if (tolower(batch_correction_method) == "harmony") {
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Package 'harmony' is required for batch_correction_method='harmony'. Please install it first.")
    }

    DefaultAssay(seu) <- "RNA"
    s_features <- intersect(cc$s.genes, rownames(seu))
    g2m_features <- intersect(cc$g2m.genes, rownames(seu))

    msg("Preparing global SCT model for Harmony.")
    seu <- NormalizeData(seu, verbose = FALSE)
    regress_vars <- c("percent.mt")
    if (length(s_features) >= 10 && length(g2m_features) >= 10) {
      seu <- CellCycleScoring(
        seu,
        s.features = s_features,
        g2m.features = g2m_features,
        set.ident = FALSE
      )
      regress_vars <- c(regress_vars, "S.Score", "G2M.Score")
    } else {
      warning("Skipping cell-cycle regression (insufficient S/G2M genes).")
    }

    seu <- SCTransform(
      seu,
      assay = "RNA",
      vars.to.regress = regress_vars,
      ncells = sct_ncells,
      conserve.memory = TRUE,
      verbose = FALSE
    )
    DefaultAssay(seu) <- "SCT"
    seu <- RunPCA(seu, npcs = max(sct_anchor_dims), verbose = FALSE)
    seu <- harmony::RunHarmony(
      object = seu,
      group.by.vars = integration_batch_col,
      reduction.use = "pca",
      dims.use = sct_anchor_dims,
      assay.use = "SCT",
      reduction.save = "harmony",
      verbose = FALSE
    )
    seu <- RunUMAP(seu, reduction = "harmony", dims = sct_umap_dims, verbose = FALSE)
  } else if (tolower(batch_correction_method) == "sct_rpca") {
    # Avoid future global-size crashes in SCT integration.
    if (requireNamespace("future", quietly = TRUE)) {
      options(future.globals.maxSize = +Inf)
      future::plan("sequential")
    }

    obj_list <- SplitObject(seu, split.by = integration_batch_col)
    obj_list <- obj_list[sapply(obj_list, ncol) > 50]
    if (length(obj_list) < 2) {
      stop("Need at least two batches with >50 cells for SCT integration.")
    }

    for (nm in names(obj_list)) {
      msg("SCT prep: ", nm, " (", ncol(obj_list[[nm]]), " cells)")
      DefaultAssay(obj_list[[nm]]) <- "RNA"

      s_features <- intersect(cc$s.genes, rownames(obj_list[[nm]]))
      g2m_features <- intersect(cc$g2m.genes, rownames(obj_list[[nm]]))

      obj_list[[nm]] <- NormalizeData(obj_list[[nm]], verbose = FALSE)
      regress_vars <- c("percent.mt")

      if (length(s_features) >= 10 && length(g2m_features) >= 10) {
        obj_list[[nm]] <- CellCycleScoring(
          obj_list[[nm]],
          s.features = s_features,
          g2m.features = g2m_features,
          set.ident = FALSE
        )
        regress_vars <- c(regress_vars, "S.Score", "G2M.Score")
      } else {
        warning("Skipping cell-cycle regression for batch ", nm, " (insufficient S/G2M genes).")
      }

      obj_list[[nm]] <- SCTransform(
        obj_list[[nm]],
        assay = "RNA",
        vars.to.regress = regress_vars,
        ncells = sct_ncells,
        conserve.memory = TRUE,
        verbose = FALSE
      )
    }

    msg("Selecting integration features.")
    integ_features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = sct_integration_nfeatures)
    obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = integ_features, verbose = FALSE)

    # Required for FindIntegrationAnchors(reduction = "rpca"):
    # each split object must already contain a PCA reduction.
    msg("Running PCA per batch for RPCA anchors.")
    obj_list <- lapply(obj_list, function(x) {
      DefaultAssay(x) <- "SCT"
      RunPCA(
        object = x,
        features = integ_features,
        npcs = max(sct_anchor_dims),
        verbose = FALSE
      )
    })

    msg("Finding SCT anchors with RPCA reduction.")
    anchors <- FindIntegrationAnchors(
      object.list = obj_list,
      normalization.method = "SCT",
      anchor.features = integ_features,
      reduction = "rpca",
      dims = sct_anchor_dims,
      verbose = FALSE
    )

    msg("Integrating data (SCT).")
    seu <- IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      dims = sct_anchor_dims,
      verbose = FALSE
    )

    DefaultAssay(seu) <- "integrated"
    seu <- RunPCA(seu, npcs = max(sct_anchor_dims), verbose = FALSE)
    seu <- RunUMAP(seu, dims = sct_umap_dims, reduction = "pca", verbose = FALSE)
  } else {
    stop("Unknown batch_correction_method: ", batch_correction_method, ". Use 'harmony' or 'sct_rpca'.")
  }
} else {
  msg("Running non-integrated NormalizeData/FindVariableFeatures/ScaleData/PCA/UMAP.")
  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_variable_features, verbose = FALSE)
  seu <- ScaleData(seu, vars.to.regress = "percent.mt", verbose = FALSE)
  seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  seu <- RunUMAP(seu, dims = umap_dims, reduction = "pca", verbose = FALSE)
}

# ------------------------------
# Plots
# ------------------------------
p_stage <- DimPlot(seu, reduction = "umap", group.by = "stage", pt.size = pt_size) +
  ggtitle("Di Bella processed dorsal reference: stage") +
  theme_classic(base_size = 12)

p_biosample <- DimPlot(seu, reduction = "umap", group.by = "biosample_id", pt.size = pt_size) +
  ggtitle("Di Bella processed dorsal reference: biosample_id") +
  theme_classic(base_size = 12)

p_biosample_split <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "biosample_id",
  split.by = "stage",
  pt.size = pt_size,
  ncol = 4
) +
  ggtitle("Di Bella processed dorsal reference: biosample_id split by stage") +
  theme_classic(base_size = 12)

p_gral <- DimPlot(seu, reduction = "umap", group.by = "Gral_cellType", pt.size = pt_size) +
  ggtitle("Di Bella processed dorsal reference: Gral_cellType") +
  theme_classic(base_size = 12)

p_broad <- DimPlot(seu, reduction = "umap", group.by = "celltype_broad_dorsal", pt.size = pt_size) +
  ggtitle("Di Bella processed dorsal reference: broad dorsal labels") +
  theme_classic(base_size = 12)

ggsave(out_umap_stage, p_stage, width = plot_w, height = plot_h, useDingbats = FALSE)
ggsave(out_umap_biosample, p_biosample, width = plot_w, height = plot_h, useDingbats = FALSE)
ggsave(out_umap_biosample_split_stage, p_biosample_split, width = plot_w * 1.6, height = plot_h * 1.2, useDingbats = FALSE)
ggsave(out_umap_gral, p_gral, width = plot_w, height = plot_h, useDingbats = FALSE)
ggsave(out_umap_broad, p_broad, width = plot_w, height = plot_h, useDingbats = FALSE)

# ------------------------------
# Save object
# ------------------------------
save_target <- if (isTRUE(do_sct_integration)) out_rds_integrated else out_rds
saveRDS(seu, save_target)
msg("Done.")
msg("Saved object: ", save_target)
msg("Saved tables in: ", tables_dir)
msg("Saved plots in: ", plots_dir)
