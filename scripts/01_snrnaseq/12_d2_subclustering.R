#!/usr/bin/env Rscript
# D2 spiny projection neuron subclustering.
#
# Splits the D2 class into a differentiating and a precursor-like state, then tests
# each state for a genotype difference in abundance and in transferred pseudotime.
# This separates two effects that are confounded at class level: a change in the
# mix of states, and a change in maturation within a state.
#
# Writes Table_S12_D2_subclustering.csv plus the per-cell and marker tables that
# the plotting script reads.
#
# Usage:  Rscript 12_d2_subclustering.R [out_dir]
#
# Requires PROJECT_ROOT to point at the analysis root (see 00_setup.R).

suppressMessages({
  library(Seurat)
  library(Matrix)
})

args     <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[1] else "."
root     <- Sys.getenv("PROJECT_ROOT", unset = NA)
if (is.na(root)) stop("set PROJECT_ROOT to the analysis root directory")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

SEED       <- 42
N_PCS      <- 20    # RunPCA npcs; FindNeighbors uses the first N_DIMS of these
N_DIMS     <- 15
RESOLUTION <- 0.2
D2_LABEL   <- "Striatal SPN D2 (indirect pathway)"

# Analysed scope: nuclei retained for the composition and maturation analyses.
# Abundance is expressed against this denominator, not against the D2 class,
# so a change in one state cannot be masked by a change in the other.
SCOPE <- c(WT = 8444, MUT = 10852)

# Variable features whose detection rate differs by more than this factor between
# genotypes are excluded from the PCA input. The filter is derived from the data
# rather than from a fixed gene list, so that any transcript behaving like ambient
# contamination is caught regardless of its identity. Genes are excluded from
# feature selection only; counts are untouched.
ASYMMETRY_FACTOR <- 8
PSEUDOCOUNT      <- 0.002

# Markers used to name the two states. Terminal SPN identity versus retained
# ganglionic-eminence precursor identity.
MARKERS <- list(
  "SPN differentiation" = c("Rarb", "Bcl11b", "Foxp1", "Adora2a", "Gpr88", "Ppp1r1b"),
  "LGE precursor"       = c("Six3", "Dlx1", "Tshz1", "Sox4")
)

message("loading annotated object")
seu <- readRDS(file.path(root, "E18p5_clean", "results", "objects", "03_seu_annotated.rds"))

# The saved object defaults to the SCT assay, which was used for the main clustering.
# Subclustering runs on RNA counts: SCT residuals are fitted across all cell types and
# are not the right variance model for a within-class split. Set this explicitly,
# because the result depends on it.
DefaultAssay(seu) <- "RNA"

md <- seu@meta.data

stopifnot("celltype_fine" %in% colnames(md))
d2_cells <- colnames(seu)[md$celltype_fine == D2_LABEL]
message(sprintf("D2 nuclei: %d", length(d2_cells)))

set.seed(SEED)
sub <- subset(seu, cells = d2_cells)
sub <- NormalizeData(sub, verbose = FALSE)
sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)

# Derive the exclusion set: variable features whose detection rate is strongly
# genotype-asymmetric are technical rather than biological, and would otherwise
# split the subclustering by library.
gw    <- colnames(sub)[md[colnames(sub), "condition"] %in% c("wt", "WT")]
gm    <- setdiff(colnames(sub), gw)
expr0 <- GetAssayData(sub, layer = "data")
hvg   <- VariableFeatures(sub)
pw    <- Matrix::rowMeans(expr0[hvg, gw, drop = FALSE] > 0)
pm    <- Matrix::rowMeans(expr0[hvg, gm, drop = FALSE] > 0)
ratio <- (pm + PSEUDOCOUNT) / (pw + PSEUDOCOUNT)
amb   <- hvg[ratio > ASYMMETRY_FACTOR | ratio < 1 / ASYMMETRY_FACTOR]

message(sprintf("variable features with >%dx genotype detection asymmetry (excluded): %s",
                ASYMMETRY_FACTOR,
                if (length(amb)) paste(amb, collapse = ", ") else "none"))
message(sprintf("  %d of %d excluded", length(amb), length(hvg)))

VariableFeatures(sub) <- setdiff(hvg, amb)

sub <- ScaleData(sub, verbose = FALSE)
sub <- RunPCA(sub, npcs = N_PCS, verbose = FALSE)
sub <- FindNeighbors(sub, dims = seq_len(N_DIMS), verbose = FALSE)
sub <- FindClusters(sub, resolution = RESOLUTION, verbose = FALSE)

# Take the assignment from the resolution-specific column rather than Idents(),
# which can carry a stale factor inherited from the parent object.
res_col <- sprintf("%s_snn_res.%s", DefaultAssay(sub), RESOLUTION)
if (!res_col %in% colnames(sub@meta.data))
  stop(sprintf("expected cluster column %s not found", res_col))
clus <- as.character(sub@meta.data[[res_col]])
if (length(unique(clus)) != 2)
  warning(sprintf("expected 2 subclusters at resolution %.2f, got %d",
                  RESOLUTION, length(unique(clus))))

# Name the states from marker expression rather than from cluster number, which is
# not stable across Seurat versions. The state with lower terminal-SPN detection is
# the precursor-like one.
expr     <- GetAssayData(sub, layer = "data")
spn_gset <- intersect(MARKERS[["SPN differentiation"]], rownames(expr))
spn_frac <- sapply(unique(clus), function(k)
  mean(colMeans(expr[spn_gset, clus == k, drop = FALSE] > 0)))
prec_id  <- names(which.min(spn_frac))

state <- ifelse(clus == prec_id, "precursor-like", "differentiating")
geno  <- ifelse(md[colnames(sub), "condition"] %in% c("wt", "WT"), "WT", "MUT")

message(sprintf("terminal-SPN detection by cluster: %s",
                paste(sprintf("%s=%.3f", names(spn_frac), spn_frac), collapse = ", ")))

# Transferred pseudotime comes from the ventral LGE/SPN projection output, joined by
# cell id. It is not a column on the annotated object.
map_f <- file.path(root, "E18p5_clean", "results", "projection_09_rooted_v2",
                   "mapping", "tables", "09b_ventral_lge_spn_mapped_table.csv")
if (!file.exists(map_f))
  stop(sprintf("mapped table not found: %s\n  run the projection pipeline first", map_f))
mp <- read.csv(map_f, stringsAsFactors = FALSE)
pt <- mp$pseudotime_transfer[match(colnames(sub), mp$cell_id)]
message(sprintf("pseudotime joined for %d of %d nuclei", sum(!is.na(pt)), length(pt)))

# Per-cell table for the plotting script.
umap <- Embeddings(seu, "umap")[colnames(sub), ]
write.csv(data.frame(cell = colnames(sub), UMAP1 = umap[, 1], UMAP2 = umap[, 2],
                     subcluster = state, geno = geno, pt = pt),
          file.path(out_dir, "d2_subcluster_cells.csv"), row.names = FALSE)

# Background nuclei for the UMAP context layer.
set.seed(SEED)
bg <- Embeddings(seu, "umap")
bg <- bg[sample(nrow(bg), min(12000, nrow(bg))), ]
write.csv(data.frame(UMAP1 = bg[, 1], UMAP2 = bg[, 2]),
          file.path(out_dir, "d2_umap_background.csv"), row.names = FALSE)

# Marker detection rates per state.
mrows <- NULL
for (grp in names(MARKERS)) {
  for (g in intersect(MARKERS[[grp]], rownames(expr))) {
    mrows <- rbind(mrows, data.frame(
      group = grp, gene = g,
      diff  = 100 * mean(expr[g, state == "differentiating"] > 0),
      prec  = 100 * mean(expr[g, state == "precursor-like"]  > 0)))
  }
}
write.csv(mrows, file.path(out_dir, "d2_subcluster_markers.csv"), row.names = FALSE)

# Abundance and maturation tests per state.
out <- NULL
for (s in c("precursor-like", "differentiating")) {
  nw <- sum(state == s & geno == "WT")
  nm <- sum(state == s & geno == "MUT")
  pw <- 100 * nw / SCOPE[["WT"]]
  pm <- 100 * nm / SCOPE[["MUT"]]
  ab <- fisher.test(matrix(c(nm, SCOPE[["MUT"]] - nm,
                             nw, SCOPE[["WT"]]  - nw), nrow = 2))
  vw <- pt[state == s & geno == "WT"];  vw <- vw[!is.na(vw)]
  vm <- pt[state == s & geno == "MUT"]; vm <- vm[!is.na(vm)]
  mt <- wilcox.test(vm, vw)
  out <- rbind(out, data.frame(
    subcluster = s, n_wt = nw, n_mut = nm,
    pct_wt = round(pw, 3), pct_mut = round(pm, 3),
    abundance_ratio = round(pm / pw, 2), abundance_p = signif(ab$p.value, 3),
    pt_median_wt = round(median(vw), 2), pt_median_mut = round(median(vm), 2),
    pt_delta = round(median(vm) - median(vw), 2), pt_p = signif(mt$p.value, 3)))
}

hdr <- c(
  sprintf("# D2 SPN subclustering. resolution %.2f on the RNA assay, %d PCs, seed %d.",
          RESOLUTION, N_PCS, SEED),
  sprintf("# Abundance normalised to analysed-scope nuclei: WT %d, MUT %d.",
          SCOPE[["WT"]], SCOPE[["MUT"]]),
  "# abundance_p: Fisher exact test on nucleus counts against the analysed scope.",
  "# pt_p: Wilcoxon rank-sum test on per-cell transferred pseudotime, within state.",
  "# State names are assigned from terminal-SPN marker detection, not cluster number.")
f <- file.path(out_dir, "Table_S12_D2_subclustering.csv")
writeLines(hdr, f)
suppressWarnings(write.table(out, f, sep = ",", row.names = FALSE, append = TRUE))

print(out, row.names = FALSE)
message(sprintf("\nwrote %s", f))
