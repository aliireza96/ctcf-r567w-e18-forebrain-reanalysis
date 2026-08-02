# ============================================================
# 03_markers_annotation.R — markers + manual annotation scaffolding
# ============================================================

source("scripts/00_setup.R")
seu <- readRDS(file.path(OUT$objects, "02_seu_SCT_clustered_res0.5.rds"))

# Prepare SCT markers
seu <- PrepSCTFindMarkers(seu)

# Find markers
markers <- FindAllMarkers(
  seu, assay = "SCT", only.pos = TRUE,
  min.pct = 0.2, logfc.threshold = 0.25
)

write_csv(markers, file.path(OUT$tables, "03_all_markers.csv"))

# Top markers heatmap (cap features for readability)
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

p_hm <- DoHeatmap(seu, features = unique(top10$gene), assay = "SCT") + NoLegend() +
  ggtitle("Top markers per cluster")
save_plot(p_hm, "dotplots/Heatmap_top10_markers", w=20, h=18)

# Canonical marker dotplot (edit list as needed)
canonical <- c("Slc17a7","Gad1","Gad2","Foxp2","Tbr1","Satb2","Eomes",
               "Sox2","Mki67","Aqp4","Pdgfra","Olig1","Olig2","Cx3cr1","P2ry12",
               "Pecam1","Rgs5")
canonical <- canonical[canonical %in% rownames(seu)]

p_dot <- DotPlot(seu, features = canonical, scale = FALSE) + RotatedAxis() + ggtitle("Canonical markers")
save_plot(p_dot, "dotplots/Dotplot_canonical_markers", w=12, h=6)

# ---- Save top 50 markers per cluster (annotation-friendly) ----
# Rationale for pct.1 filter:
# pct.1 is the fraction of cells in the cluster expressing the gene.
# Requiring pct.1 >= 0.25 prioritizes markers that are more consistently expressed
# across the cluster (reduces "spiky"/rare markers that are harder to interpret for annotation).

top_n <- 50

top50_markers <- markers %>%
  dplyr::filter(p_val_adj < 0.05, pct.1 >= 0.25) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(cluster, dplyr::desc(avg_log2FC))

readr::write_csv(
  top50_markers,
  file.path(OUT$tables, paste0("03_top", top_n, "_markers_per_cluster_padj0.05_pct1ge0.25.csv"))
)

# --------------------------
# Manual annotation section
# --------------------------
# Create a template mapping table to edit
cluster_levels <- levels(seu$seurat_clusters)
annot_template <- data.frame(
  seurat_cluster = cluster_levels,
  celltype_fine = NA_character_,
  celltype_broad = NA_character_
)

write_csv(annot_template, file.path(OUT$tables, "03_annotation_template_fillme.csv"))

message("Fill in 03_annotation_template_fillme.csv, then rerun annotation step below.")

# ---- If you have your mapping already, load it ----
# Example: use your prior mapping (replace with your final mapping)
# mapping <- read_csv(file.path(OUT$tables, "03_annotation_template_fillme.csv"))
# seu$celltype_fine  <- plyr::mapvalues(as.character(seu$seurat_clusters), mapping$seurat_cluster, mapping$celltype_fine)
# seu$celltype_broad <- plyr::mapvalues(as.character(seu$seurat_clusters), mapping$seurat_cluster, mapping$celltype_broad)

saveRDS(seu, file.path(OUT$objects, "03_seu_markers_ready_for_annotation.rds"))
