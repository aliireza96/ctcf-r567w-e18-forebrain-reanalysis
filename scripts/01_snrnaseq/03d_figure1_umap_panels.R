# ============================================================
# 03d_figure1_umap_panels.R
# Figure 1 UMAP panels, same plot idiom as 03c_dissection_bias_and_tel_filter.R.
#
# Changes relative to 03c, requested for the figure:
#   1. grouping: upper-layer EN merged, deep-layer EN merged, SPN-D1 / SPN-D2 kept separate
#      (i.e. celltype_broad, but with the SPN merge undone)
#   2. aspect ratio: near-square per UMAP (03c used w=15,h=8 which reads too wide)
#
# Outputs:
#   plots/umaps/03d_UMAP_figure1_celltype_scope.pdf   (Figure 1C)
#   plots/umaps/03d_UMAP_figure1_region_all.pdf       (Figure 1A, region, square)
#   plots/umaps/03d_UMAP_figure1_celltype_all.pdf     (optional companion, all cells by celltype)
# ============================================================

source("scripts/00_setup.R")

# Load the annotated objects first; the region palette below is derived from the
# object's own factor levels exactly as in 03c.
seu_all <- readRDS(file.path(OUT$objects, "03_seu_annotated.rds"))

# Reuse the exact boxed-label helper from 03c so the panels keep the same look.
# Evaluate only the function definition, not the whole script, to avoid
# re-running the filtering and re-writing the objects.
src_03c <- readLines("scripts/03c_dissection_bias_and_tel_filter.R")
fn_start <- grep("^add_boxed_labels <- function", src_03c)[1]
fn_end   <- grep("^\\}$", src_03c)
fn_end   <- fn_end[fn_end > fn_start][1]
eval(parse(text = paste(src_03c[fn_start:fn_end], collapse = "\n")))
stopifnot(exists("add_boxed_labels"))

# Single shared region palette, same construction as 03c
region_levels <- levels(seu_all$anatomical_region_display)
region_colors <- stats::setNames(
  grDevices::hcl.colors(length(region_levels), palette = "Dark 3"),
  region_levels
)

# Point rasterisation. A UMAP of ~29,000 nuclei is ~29,000 vector marks, which makes the
# PDF slow to open in Illustrator for no visual benefit at print size. rasterise() converts
# ONLY the point layer to an embedded bitmap at print resolution; axes, tick labels, titles
# and the boxed cluster labels stay live vector text and remain editable.
RASTER_DPI <- 450
rast_points <- function(p) {
  if (requireNamespace("ggrastr", quietly = TRUE)) {
    ggrastr::rasterise(p, layers = "Point", dpi = RASTER_DPI, dev = "ragg")
  } else p
}

upper_labels <- c("IT-L2/3", "IT-L2/4")
deep_labels  <- c("Corticofugal", "Deep corticofugal", "Deep cortical")

add_figure_grouping <- function(obj) {
  lab <- as.character(obj$celltype_label)
  obj$celltype_figure <- dplyr::case_when(
    lab %in% upper_labels ~ "Upper layer EN",
    lab %in% deep_labels  ~ "Deep layer EN",
    TRUE ~ lab
  )
  obj
}

# ---- Figure 1C: annotated UMAP within the cortical/striatal analysis scope ----
seu_cortstr <- readRDS(file.path(OUT$objects, "03_seu_annotated_CORTSTR.rds"))
seu_cortstr <- add_figure_grouping(seu_cortstr)
message("Figure grouping levels (scope): ",
        paste(sort(unique(seu_cortstr$celltype_figure)), collapse = ", "))

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
p_fig1c <- DimPlot(
  seu_cortstr, reduction = "umap", group.by = "celltype_figure",
  label = !has_ggrepel, repel = !has_ggrepel, pt.size = 0.25
) +
  ggtitle("Cortical + Striatal analysis scope")
# rasterise the cell points before the boxed labels are added, so the labels stay vector
p_fig1c <- rast_points(p_fig1c)
if (has_ggrepel) {
  p_fig1c <- add_boxed_labels(p_fig1c, seu_cortstr, "celltype_figure", label.size = 3.0)
}
# Labels sit on the clusters (centroid label purity 0.98-1.00 for every group), so the
# colour legend is redundant; dropping it frees the space the wide layout was wasting.
p_fig1c <- p_fig1c + theme(legend.position = "none")
# near-square plotting area; 03c used w=15,h=8
save_plot(p_fig1c, "umaps/03d_UMAP_figure1_celltype_scope", w = 7, h = 6.5)

# ---- Figure 1A: all cells coloured by anatomical region (square) -------------
# The region panel KEEPS its legend. Region centroids are compact enough to label
# (purity 0.96-1.00), but "Cortical + Striatal" has a median centroid spread of 6.5
# UMAP units versus 0.2-1.8 for every other region: a single centroid label would sit
# in the middle of a region spanning most of the plot and would read as pointing at
# one lobe of it. A legend is the honest encoding for this panel.
p_fig1a <- DimPlot(
  seu_all, reduction = "umap", group.by = "anatomical_region_display",
  cols = unname(region_colors[region_levels]), pt.size = 0.25
) +
  ggtitle("All cells colored by anatomical region")
p_fig1a <- rast_points(p_fig1a)
# The 9-entry legend takes roughly 2.9 in of the width, so the plotting panel is much
# narrower than the canvas. UMAP data aspect here is 27.7 / 23.0 = 1.20; at w = 8 the
# panel came out at 0.84 (visibly squished). w = 10 puts the panel at 1.16, matching
# the data. The labelled panels below need no such allowance (no legend).
save_plot(p_fig1a, "umaps/03d_UMAP_figure1_region_all", w = 10, h = 7)

# ---- optional companion: all cells by cell type, figure grouping --------------
seu_all <- add_figure_grouping(seu_all)
p_fig1a2 <- DimPlot(
  seu_all, reduction = "umap", group.by = "celltype_figure",
  label = !has_ggrepel, repel = !has_ggrepel, pt.size = 0.25
) +
  ggtitle("All cells colored by cell type")
p_fig1a2 <- rast_points(p_fig1a2)
if (has_ggrepel) {
  p_fig1a2 <- add_boxed_labels(p_fig1a2, seu_all, "celltype_figure", label.size = 2.8)
}
p_fig1a2 <- p_fig1a2 + theme(legend.position = "none")
save_plot(p_fig1a2, "umaps/03d_UMAP_figure1_celltype_all", w = 7, h = 6.5)

message("03d done")
