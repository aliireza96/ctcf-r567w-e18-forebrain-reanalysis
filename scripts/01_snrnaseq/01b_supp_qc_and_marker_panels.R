# ============================================================
# 01b_supp_qc_and_marker_panels.R
# Supplementary Figure 1 panels that needed regenerating or did not exist.
#
# 1. QC violins AFTER filtering.
#    The existing qc_plots/QC_violin_after.pdf is broken: only nFeature_RNA carries a
#    facet title and the remaining metrics are clipped, so the panel is unreadable.
#    Cause: in 01_load_qc_filter_doublets.R the VlnPlot patchwork is combined with
#    plot_annotation() and saved at w=14,h=5; the same call works for the pre-filter
#    plot but not post-filter. Rebuilt here as an explicit wrap_plots() of four
#    single-feature violins, each with its own title, which is robust to that.
#
# 2. Dorsal / excitatory marker dotplot, as a companion to the existing ventral one.
#    Built with the same DotPlot idiom, grouping and theme as the ventral panel in
#    03b_apply_annotation.R, so the two read as a matched pair.
#
# Outputs:
#   plots/qc_plots/01b_QC_violin_after_filtering.pdf
#   plots/dotplots/01b_Dotplot_dorsal_markers_by_finetype.pdf
# ============================================================

source("scripts/00_setup.R")
suppressPackageStartupMessages(library(patchwork))

# ---- 1. QC violins after filtering ------------------------------------------
# 01 saves the filtered object as part of the QC step; fall back to the annotated
# object if the QC-stage object is unavailable.
qc_obj <- file.path(OUT$objects, "01_seu_QCfiltered.rds")
stopifnot(file.exists(qc_obj))
seu_f <- readRDS(qc_obj)
message("QC object: ", ncol(seu_f), " nuclei; conditions: ",
        paste(names(table(seu_f$condition)), table(seu_f$condition), collapse = " "))

qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo")
qc_feats <- qc_feats[qc_feats %in% colnames(seu_f@meta.data)]
stopifnot(length(qc_feats) == 4)

# Point rasterisation. Each violin draws one mark per nucleus, so four violins over
# 28,802 nuclei produced ~115,000 vector operations and a file that is slow to open in
# a vector editor. ggrastr::rasterise() converts ONLY the jittered-point layer to an
# embedded bitmap at print resolution; the violin outlines, axes, tick labels and titles
# remain live vector text, so the panel is still fully editable.
RASTER_DPI <- 450

vln_one <- function(f) {
  p <- Seurat::VlnPlot(seu_f, features = f, group.by = "condition", pt.size = 0.1)
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    p <- ggrastr::rasterise(p, layers = "Point", dpi = RASTER_DPI, dev = "ragg")
  }
  p +
    ggtitle(f) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          axis.title.x = element_blank())
}
p_qc <- patchwork::wrap_plots(lapply(qc_feats, vln_one), nrow = 1)
save_plot(p_qc, "qc_plots/01b_QC_violin_after_filtering", w = 12, h = 4)

# ---- 2. dorsal / excitatory marker dotplot ----------------------------------
# Same construction as the ventral panel: unscaled DotPlot over celltype_fine with
# RotatedAxis and the minimal theme, so the pair is directly comparable.
kept_obj <- file.path(OUT$objects, "03_seu_annotated_KEPT.rds")
stopifnot(file.exists(kept_obj))
seu_kept <- readRDS(kept_obj)

# Dorsal telencephalic programme: progenitor -> intermediate progenitor ->
# deep-layer -> upper-layer -> glial, plus pan-excitatory identity.
dorsal_markers <- c(
  "Pax6", "Sox2", "Hes5",           # radial glia / apical progenitors
  "Eomes", "Neurod2",               # intermediate progenitors / newborn EN
  "Tbr1", "Bcl11b", "Fezf2",        # deep layer / corticofugal
  "Satb2", "Cux1", "Rorb",          # upper layer / callosal
  "Reln",                           # Cajal-Retzius
  "Aldoc", "Slc1a3",                # astrocyte lineage
  "Olig1", "Pdgfra",                # OPC
  "Slc17a6", "Slc17a7"              # pan-glutamatergic
)
dorsal_markers <- dorsal_markers[dorsal_markers %in% rownames(seu_kept)]
message("dorsal markers found: ", length(dorsal_markers), " of 18 -> ",
        paste(dorsal_markers, collapse = ", "))
stopifnot(length(dorsal_markers) >= 10)

p_dors <- Seurat::DotPlot(seu_kept, features = dorsal_markers,
                          group.by = "celltype_fine", scale = FALSE) +
  ggtitle("Dorsal / excitatory-focused markers (fine types)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
        panel.grid = element_blank()) +
  Seurat::RotatedAxis()
save_plot(p_dors, "dotplots/01b_Dotplot_dorsal_markers_by_finetype", w = 13, h = 6)

# ---- 3. genotype overlay UMAP, rasterised points ----------------------------
# Same DimPlot call as in 03b_apply_annotation.R (line ~400), re-emitted here with the
# point layer rasterised. Regenerating it from 03b would re-run the whole annotation
# step, so it is reproduced in place instead.
seu_full <- readRDS(file.path(OUT$objects, "03_seu_annotated.rds"))
p_cond <- Seurat::DimPlot(seu_full, group.by = "condition",
                          label = FALSE, repel = FALSE, pt.size = 0.25)
if (requireNamespace("ggrastr", quietly = TRUE)) {
  p_cond <- ggrastr::rasterise(p_cond, layers = "Point", dpi = RASTER_DPI, dev = "ragg")
}
p_cond <- p_cond + ggtitle("WT vs MUT overlay (all clusters)")
if (exists("theme_pub")) p_cond <- p_cond + theme_pub()
save_plot(p_cond, "umaps/01b_UMAP_condition_overlay_rasterised", w = 8, h = 7)

message("01b done")
