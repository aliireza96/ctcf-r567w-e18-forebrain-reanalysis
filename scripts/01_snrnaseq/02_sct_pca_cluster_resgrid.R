# ============================================================
# 02_sct_pca_cluster_resgrid.R — SCT, PCA, clustering resolution grid
# ============================================================

source("scripts/00_setup.R")
seu <- readRDS(file.path(OUT$objects, "01_seu_QCfiltered.rds"))

# ---- future / parallel settings for SCT ----
suppressPackageStartupMessages(library(future))

# Use sequential to avoid huge object export to workers (most robust)
plan(sequential)

# Also raise the export limit just in case (2-8 GB are common for scRNA objects)
options(future.globals.maxSize = 8 * 1024^3)  # 8 GB

# ---- SCTransform ----
seu <- SCTransform(seu, verbose = FALSE, method = "glmGamPoi", conserve.memory = TRUE)

# ---- PCA ----
seu <- RunPCA(seu, verbose = FALSE)
p_elbow <- ElbowPlot(seu, ndims = 60) + ggtitle("PCA elbow (SCT)")
save_plot(p_elbow, "qc_plots/PCA_elbow", w=6, h=4)

dims_use <- 1:30
seu <- FindNeighbors(seu, dims = dims_use, verbose = FALSE)

# ---- Resolution grid ----
res_grid <- c(0.2, 0.4, 0.5, 0.6, 0.8)
umaps <- list()

for (r in res_grid) {
  seu <- FindClusters(seu, resolution = r, verbose = FALSE)
  seu <- RunUMAP(seu, dims = dims_use, verbose = FALSE)
  p <- DimPlot(seu, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle(paste0("Resolution = ", r)) + NoLegend()
  umaps[[as.character(r)]] <- p
}

save_plot(wrap_plots(umaps, ncol = 2), "umaps/UMAP_resolution_grid", w=12, h=10)

# ---- Choose final resolution (you prefer 0.5) ----
final_res <- 0.5
seu <- FindClusters(seu, resolution = final_res, verbose = FALSE)
seu <- RunUMAP(seu, dims = dims_use, verbose = FALSE)

p1 <- DimPlot(seu, label = TRUE, repel = TRUE) + ggtitle("Clusters (final)") + NoLegend()
p2 <- DimPlot(seu, group.by = "condition") + ggtitle("Condition")
save_plot(p1 + p2, "umaps/UMAP_clusters_condition_final", w=12, h=5)

p_split <- DimPlot(seu, group.by="seurat_clusters", split.by="condition", label=TRUE, repel=TRUE) + NoLegend()
save_plot(p_split, "umaps/UMAP_split_by_condition_clusters", w=14, h=6)

saveRDS(seu, file.path(OUT$objects, "02_seu_SCT_clustered_res0.5.rds"))
message("Saved: 02_seu_SCT_clustered_res0.5.rds")
