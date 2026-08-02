#!/usr/bin/env Rscript

# ============================================================
# 07_build_ventral_reference_mayer.R
# Build ventral telencephalon developmental reference from
# Mayer et al. GSE103983 (Drop-seq GE dataset)
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# ------------------------------
# Parameters (edit here)
# ------------------------------
# GSE103983 (Mayer et al. 2018). Set REF_MAYER to the downloaded file.
input_file <- Sys.getenv("REF_MAYER", unset = "data/references/GSE103983_dropseq.csv.gz")
output_dir <- file.path(ROOT, "results/mayer_ventral_reference")

# Read mode:
# - "full" : read all cells at once
# - "sequential_by_region" : read MGE/CGE/LGE columns separately and cbind
read_mode <- "sequential_by_region"
regions_to_keep <- c("MGE", "CGE", "LGE")

# Core QC
min_features <- 700L

# Optional QC toggles
use_log10_outlier_filter <- TRUE
log10_sd_threshold <- 3

use_loess_residual_filter <- TRUE
loess_sd_threshold <- 3
loess_span <- 0.75
min_cells_for_sample_outlier_model <- 100L

# Optional contaminant removal
remove_contaminants <- TRUE
contam_genes <- c("Neurod6", "Igfbp7")

# Normalize / dimensional reduction
n_variable_features <- 3000L
n_pcs <- 50L
umap_dims <- 1:30

# Plot settings
pt_size <- 0.25
plot_width <- 9
plot_height <- 7

# ------------------------------
# Paths
# ------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(output_dir, "objects", "mayer_GSE103983_ventral_reference_clean.rds")
out_qc <- file.path(output_dir, "tables", "mayer_GSE103983_qc_summary_by_sample.csv")
out_qc_cell <- file.path(output_dir, "tables", "mayer_GSE103983_qc_cell_flags.csv.gz")
out_meta <- file.path(output_dir, "tables", "mayer_GSE103983_metadata_allcells.csv.gz")
out_plot_region <- file.path(output_dir, "plots", "UMAP_mayer_reference_by_region.pdf")
out_plot_sample <- file.path(output_dir, "plots", "UMAP_mayer_reference_by_sample.pdf")

# ------------------------------
# Helpers
# ------------------------------
message2 <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")

safe_log10 <- function(x) log10(pmax(x, 1))

parse_cell_metadata <- function(cell_ids) {
  sample <- sub("_[^_]+$", "", cell_ids)
  barcode <- sub("^.*_", "", cell_ids)
  region <- toupper(sub("_.*$", "", cell_ids))
  stage <- ifelse(region == "MGE", "E13.5",
                  ifelse(region %in% c("CGE", "LGE"), "E14.5", NA_character_))
  data.frame(
    cell_id = cell_ids,
    sample = sample,
    barcode = barcode,
    region = region,
    stage = stage,
    stringsAsFactors = FALSE
  )
}

find_gene_name <- function(gene, gene_universe) {
  if (gene %in% gene_universe) return(gene)
  hit <- gene_universe[toupper(gene_universe) == toupper(gene)]
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

get_counts_layer <- function(seu_obj) {
  mat <- tryCatch(GetAssayData(seu_obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  mat <- tryCatch(GetAssayData(seu_obj, assay = "RNA", slot = "counts"), error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  mat <- tryCatch(LayerData(seu_obj, assay = "RNA", layer = "counts"), error = function(e) NULL)
  if (!is.null(mat)) return(mat)
  stop("Could not access RNA counts layer.")
}

filter_log10_outliers_by_sample <- function(meta_df, sd_thresh = 3, min_cells = 100) {
  keep <- rep(TRUE, nrow(meta_df))
  names(keep) <- rownames(meta_df)
  for (s in unique(meta_df$sample)) {
    idx <- which(meta_df$sample == s)
    if (length(idx) < min_cells) next
    x <- safe_log10(meta_df$nCount_RNA[idx])
    y <- safe_log10(meta_df$nFeature_RNA[idx])
    x_ok <- abs(x - mean(x, na.rm = TRUE)) <= sd_thresh * sd(x, na.rm = TRUE)
    y_ok <- abs(y - mean(y, na.rm = TRUE)) <= sd_thresh * sd(y, na.rm = TRUE)
    keep[idx] <- x_ok & y_ok
  }
  keep
}

filter_loess_residual_outliers_by_sample <- function(meta_df, sd_thresh = 3, span = 0.75, min_cells = 100) {
  keep <- rep(TRUE, nrow(meta_df))
  names(keep) <- rownames(meta_df)
  for (s in unique(meta_df$sample)) {
    idx <- which(meta_df$sample == s)
    if (length(idx) < min_cells) next
    x <- safe_log10(meta_df$nCount_RNA[idx])
    y <- safe_log10(meta_df$nFeature_RNA[idx])
    fit <- tryCatch(loess(y ~ x, span = span), error = function(e) NULL)
    if (is.null(fit)) next
    resid <- y - stats::predict(fit, x)
    rmu <- mean(resid, na.rm = TRUE)
    rsd <- sd(resid, na.rm = TRUE)
    if (!is.finite(rsd) || rsd == 0) next
    keep[idx] <- abs(resid - rmu) <= sd_thresh * rsd
  }
  keep
}

to_sparse_counts <- function(dt, gene_col) {
  genes <- as.character(dt[[gene_col]])
  expr_cols <- setdiff(colnames(dt), gene_col)
  m <- as.matrix(dt[, expr_cols, with = FALSE])
  storage.mode(m) <- "numeric"
  sp <- Matrix(m, sparse = TRUE)
  rownames(sp) <- genes
  colnames(sp) <- expr_cols
  sp
}

read_header_names <- function(file) {
  hdr <- fread(
    cmd = sprintf("zcat %s", shQuote(file)),
    nrows = 0,
    sep = ",",
    check.names = FALSE,
    data.table = TRUE
  )
  colnames(hdr)
}

read_subset_sparse <- function(file, gene_col, selected_cells) {
  if (length(selected_cells) == 0) return(NULL)
  dt <- fread(
    cmd = sprintf("zcat %s", shQuote(file)),
    sep = ",",
    select = c(gene_col, selected_cells),
    check.names = FALSE,
    data.table = TRUE,
    showProgress = TRUE
  )
  to_sparse_counts(dt, gene_col)
}

message2("Reading input matrix: ", input_file)
stopifnot(file.exists(input_file))

# ------------------------------
# Read matrix (memory-aware)
# ------------------------------
col_names <- read_header_names(input_file)
if (length(col_names) < 2) stop("Input appears malformed: expected gene column + cell columns.")

gene_col <- col_names[1]
all_cells <- col_names[-1]
all_meta <- parse_cell_metadata(all_cells)

# Keep target regions only (MGE/CGE/LGE by default)
target_cells <- all_meta$cell_id[all_meta$region %in% regions_to_keep]
if (length(target_cells) == 0) stop("No cells matched requested regions: ", paste(regions_to_keep, collapse = ", "))

counts <- NULL
if (read_mode == "full") {
  message2("Read mode = full")
  dt <- fread(
    cmd = sprintf("zcat %s", shQuote(input_file)),
    sep = ",",
    select = c(gene_col, target_cells),
    check.names = FALSE,
    data.table = TRUE,
    showProgress = TRUE
  )
  counts <- to_sparse_counts(dt, gene_col)
} else if (read_mode == "sequential_by_region") {
  message2("Read mode = sequential_by_region")
  region_mats <- list()
  for (rg in regions_to_keep) {
    rg_cells <- all_meta$cell_id[all_meta$region == rg]
    if (length(rg_cells) == 0) next
    message2("Reading region ", rg, " (", length(rg_cells), " cells)")
    sp <- read_subset_sparse(input_file, gene_col, rg_cells)
    if (is.null(sp)) next
    region_mats[[rg]] <- sp
  }
  if (length(region_mats) == 0) stop("No region matrices were loaded.")
  common_genes <- Reduce(intersect, lapply(region_mats, rownames))
  if (length(common_genes) == 0) stop("No common genes across region-wise reads.")
  region_mats <- lapply(region_mats, function(m) m[common_genes, , drop = FALSE])
  # Matrix::cBind is deprecated/defunct in newer Matrix versions; use Reduce(cbind2).
  counts <- Reduce(Matrix::cbind2, region_mats)
  # Keep only requested cell order
  target_cells <- target_cells[target_cells %in% colnames(counts)]
  counts <- counts[, target_cells, drop = FALSE]
} else {
  stop("Unsupported read_mode: ", read_mode, ". Use 'full' or 'sequential_by_region'.")
}

# Ensure gene names are unique for Seurat compatibility
if (anyDuplicated(rownames(counts)) > 0) {
  message2("Duplicated gene symbols found: making unique names.")
  rownames(counts) <- make.unique(rownames(counts))
}

message2("Counts matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

# ------------------------------
# Build Seurat object + metadata
# ------------------------------
seu <- CreateSeuratObject(counts = counts, project = "Mayer_GSE103983_GE", min.cells = 0, min.features = 0)
meta <- parse_cell_metadata(colnames(seu))
rownames(meta) <- meta$cell_id
seu$cell_id <- meta$cell_id
seu$sample <- meta$sample
seu$barcode <- meta$barcode
seu$region <- meta$region
seu$stage <- meta$stage
seu$dataset <- "Mayer_GSE103983"

write_csv(meta, out_meta)

# ------------------------------
# QC filters
# ------------------------------
meta0 <- seu@meta.data
qc_flags <- data.frame(
  cell_id = rownames(meta0),
  sample = meta0$sample,
  region = meta0$region,
  stage = meta0$stage,
  nCount_RNA = meta0$nCount_RNA,
  nFeature_RNA = meta0$nFeature_RNA,
  pass_min_features = meta0$nFeature_RNA >= min_features,
  pass_log10_outlier = TRUE,
  pass_loess_residual = TRUE,
  pass_qc = TRUE,
  stringsAsFactors = FALSE
)
rownames(qc_flags) <- qc_flags$cell_id

if (use_log10_outlier_filter) {
  keep_log <- filter_log10_outliers_by_sample(meta0, sd_thresh = log10_sd_threshold, min_cells = min_cells_for_sample_outlier_model)
  qc_flags[names(keep_log), "pass_log10_outlier"] <- keep_log
}

if (use_loess_residual_filter) {
  keep_loess <- filter_loess_residual_outliers_by_sample(
    meta0,
    sd_thresh = loess_sd_threshold,
    span = loess_span,
    min_cells = min_cells_for_sample_outlier_model
  )
  qc_flags[names(keep_loess), "pass_loess_residual"] <- keep_loess
}

qc_flags$pass_qc <- qc_flags$pass_min_features & qc_flags$pass_log10_outlier & qc_flags$pass_loess_residual
cells_after_qc <- rownames(qc_flags)[qc_flags$pass_qc]
message2("Cells after QC: ", length(cells_after_qc), " / ", ncol(seu))
seu_qc <- subset(seu, cells = cells_after_qc)

# ------------------------------
# Remove contaminants (optional)
# ------------------------------
qc_flags$pass_contaminant_filter <- TRUE
if (remove_contaminants) {
  counts_qc <- get_counts_layer(seu_qc)
  row_universe <- rownames(counts_qc)
  neurod6_name <- find_gene_name(contam_genes[1], row_universe)
  igfbp7_name <- find_gene_name(contam_genes[2], row_universe)

  neurod6_pos <- rep(FALSE, ncol(seu_qc))
  igfbp7_pos <- rep(FALSE, ncol(seu_qc))
  names(neurod6_pos) <- colnames(seu_qc)
  names(igfbp7_pos) <- colnames(seu_qc)

  if (!is.na(neurod6_name)) neurod6_pos <- as.vector(counts_qc[neurod6_name, ] > 0)
  if (!is.na(igfbp7_name)) igfbp7_pos <- as.vector(counts_qc[igfbp7_name, ] > 0)

  contam <- neurod6_pos | igfbp7_pos
  keep_cells <- colnames(seu_qc)[!contam]
  seu_clean <- subset(seu_qc, cells = keep_cells)

  qc_flags[colnames(seu_qc), "pass_contaminant_filter"] <- !contam
  message2("Contaminant-filtered cells: ", ncol(seu_clean), " / ", ncol(seu_qc))
} else {
  seu_clean <- seu_qc
}

qc_flags$pass_final <- qc_flags$pass_qc & qc_flags$pass_contaminant_filter

# ------------------------------
# QC summary table
# ------------------------------
qc_summary <- qc_flags %>%
  group_by(sample, region, stage) %>%
  summarise(
    n_before = n(),
    n_after_min_features = sum(pass_min_features),
    n_after_log10_filter = sum(pass_min_features & pass_log10_outlier),
    n_after_loess_filter = sum(pass_qc),
    n_after_contaminant_filter = sum(pass_final),
    frac_retained = n_after_contaminant_filter / n_before,
    .groups = "drop"
  ) %>%
  arrange(region, sample)

write_csv(qc_summary, out_qc)
write_csv(qc_flags, out_qc_cell)

# ------------------------------
# Normalize + regress cell cycle + PCA/UMAP
# ------------------------------
DefaultAssay(seu_clean) <- "RNA"
seu_clean <- NormalizeData(seu_clean, verbose = FALSE)
seu_clean <- FindVariableFeatures(seu_clean, selection.method = "vst", nfeatures = n_variable_features, verbose = FALSE)

cc <- Seurat::cc.genes.updated.2019
s_genes <- intersect(cc$s.genes, rownames(seu_clean))
g2m_genes <- intersect(cc$g2m.genes, rownames(seu_clean))

if (length(s_genes) >= 10 && length(g2m_genes) >= 10) {
  seu_clean <- CellCycleScoring(
    seu_clean,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = FALSE
  )
  seu_clean <- ScaleData(seu_clean, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
} else {
  warning("Not enough cell-cycle genes found in dataset. Skipping CellCycleScoring/CC regression.")
  seu_clean <- ScaleData(seu_clean, verbose = FALSE)
}

seu_clean <- RunPCA(seu_clean, npcs = n_pcs, verbose = FALSE)
seu_clean <- RunUMAP(seu_clean, dims = umap_dims, reduction = "pca", verbose = FALSE)

# ------------------------------
# Save UMAP plots
# ------------------------------
p_region <- DimPlot(seu_clean, reduction = "umap", group.by = "region", pt.size = pt_size) +
  ggtitle("Mayer GSE103983 reference: UMAP by region") +
  theme_classic(base_size = 12)

p_sample <- DimPlot(seu_clean, reduction = "umap", group.by = "sample", pt.size = pt_size) +
  ggtitle("Mayer GSE103983 reference: UMAP by sample") +
  theme_classic(base_size = 12)

ggsave(out_plot_region, plot = p_region, width = plot_width, height = plot_height, useDingbats = FALSE)
ggsave(out_plot_sample, plot = p_sample, width = plot_width, height = plot_height, useDingbats = FALSE)

# ------------------------------
# Save object
# ------------------------------
saveRDS(seu_clean, out_rds)

message2("Done.")
message2("Saved RDS: ", out_rds)
message2("Saved QC summary: ", out_qc)
message2("Saved UMAPs: ", out_plot_region, " ; ", out_plot_sample)
