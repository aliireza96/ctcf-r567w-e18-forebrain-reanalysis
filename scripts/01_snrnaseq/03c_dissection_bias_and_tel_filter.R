# ============================================================
# 03c_dissection_bias_and_tel_filter.R
# Purpose:
#   - Add/validate region metadata from region-aware master table
#   - Produce dissection-bias figures/tables
#   - Write cortical/striatal-focused object for downstream analyses
#   - Write a sensitivity object that also includes Hypothalamus_glutamatergic
# ============================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

in_obj <- file.path(OUT$objects, "03_seu_annotated.rds")
in_region_map <- file.path(OUT$tables, "03_annotation_master_table_with_regions.csv")
out_obj_cortstr <- file.path(OUT$objects, "03_seu_annotated_CORTSTR.rds")
out_obj_cortstr_hypo <- file.path(OUT$objects, "03_seu_annotated_CORTSTR_plus_HYPOGLUT.rds")
out_obj_tel_legacy <- file.path(OUT$objects, "03_seu_annotated_TEL.rds")

stopifnot(file.exists(in_obj))
stopifnot(file.exists(in_region_map))

seu <- readRDS(in_obj)
region_map <- readr::read_csv(in_region_map, show_col_types = FALSE)

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

req_cols <- c("Cluster Number", "Anatomical Region")
if (!all(req_cols %in% colnames(region_map))) {
  stop("Region master must contain: Cluster Number, Anatomical Region")
}

cluster_num <- suppressWarnings(as.integer(region_map[["Cluster Number"]]))
map_tbl <- data.frame(
  seurat_cluster = as.character(cluster_num),
  anatomical_region = as.character(region_map[["Anatomical Region"]]),
  is_telencephalic = if ("Telencephalic?" %in% colnames(region_map)) as_bool(region_map[["Telencephalic?"]]) else if ("Telencephalic" %in% colnames(region_map)) as_bool(region_map[["Telencephalic"]]) else NA,
  include_telencephalon = if ("Include in telencephalon proportion analysis?" %in% colnames(region_map)) as_bool(region_map[["Include in telencephalon proportion analysis?"]]) else NA,
  stringsAsFactors = FALSE
)
map_tbl <- map_tbl[!is.na(cluster_num), , drop = FALSE]
map_tbl <- map_tbl %>% dplyr::distinct(seurat_cluster, .keep_all = TRUE)

if (all(is.na(map_tbl$include_telencephalon))) {
  map_tbl$include_telencephalon <- map_tbl$is_telencephalic
}
map_tbl$include_telencephalon[is.na(map_tbl$include_telencephalon)] <- FALSE

map_region <- setNames(map_tbl$anatomical_region, map_tbl$seurat_cluster)
map_tel <- setNames(map_tbl$is_telencephalic, map_tbl$seurat_cluster)
map_tel_keep <- setNames(map_tbl$include_telencephalon, map_tbl$seurat_cluster)

seu$anatomical_region <- unname(map_region[as.character(seu$seurat_clusters)])
seu$anatomical_region_display <- dplyr::case_when(
  grepl("^telencephalon$", seu$anatomical_region, ignore.case = TRUE) ~ "Cortical + Striatal",
  TRUE ~ as.character(seu$anatomical_region)
)
seu$is_telencephalic <- unname(map_tel[as.character(seu$seurat_clusters)])
seu$include_telencephalon <- unname(map_tel_keep[as.character(seu$seurat_clusters)])
seu$include_cortical_striatal_scope <- seu$include_telencephalon

hypo_label_hit <- rep(FALSE, ncol(seu))
for (col in c("celltype_label", "celltype_broad")) {
  if (col %in% colnames(seu@meta.data)) {
    hypo_label_hit <- hypo_label_hit | as.character(seu@meta.data[[col]]) == "Hypothalamus_glutamatergic"
  }
}
seu$include_hypothalamus_glutamatergic_sensitivity <- hypo_label_hit
seu$include_cortical_striatal_plus_hypoglut_sensitivity <- seu$include_cortical_striatal_scope | seu$include_hypothalamus_glutamatergic_sensitivity

if (any(!is.na(seu$anatomical_region))) {
  reg_levels <- unique(map_tbl$anatomical_region[!is.na(map_tbl$anatomical_region) & map_tbl$anatomical_region != ""])
  reg_levels <- reg_levels[reg_levels %in% unique(as.character(seu$anatomical_region))]
  if (length(reg_levels) == 0) {
    reg_levels <- sort(unique(as.character(seu$anatomical_region)))
    reg_levels <- reg_levels[!is.na(reg_levels) & reg_levels != ""]
  }
  seu$anatomical_region <- factor(seu$anatomical_region, levels = reg_levels)
  reg_display_levels <- dplyr::case_when(
    grepl("^telencephalon$", reg_levels, ignore.case = TRUE) ~ "Cortical + Striatal",
    TRUE ~ reg_levels
  )
  seu$anatomical_region_display <- factor(seu$anatomical_region_display, levels = unique(reg_display_levels))
}

# Single shared region palette used for BOTH UMAP and barplots
region_levels <- levels(seu$anatomical_region_display)
region_colors <- stats::setNames(
  grDevices::hcl.colors(length(region_levels), palette = "Dark 3"),
  region_levels
)

# Add bold boxed labels at group centroids on UMAP
add_boxed_labels <- function(p, seu_obj, group.by, label.size = 3.2, max.overlaps = 200) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) return(p)

  emb <- as.data.frame(Seurat::Embeddings(seu_obj, reduction = "umap"))
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(emb))) {
    umap_cols <- grep("^UMAP", colnames(emb), value = TRUE, ignore.case = TRUE)
    if (length(umap_cols) >= 2) {
      emb <- emb[, umap_cols[1:2], drop = FALSE]
      colnames(emb) <- c("UMAP_1", "UMAP_2")
    } else if (ncol(emb) >= 2) {
      colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
    } else {
      stop("UMAP embeddings not found or have unexpected shape.")
    }
  }

  emb$grp <- as.character(seu_obj[[group.by, drop = TRUE]])
  centers <- emb %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )

  p + ggrepel::geom_label_repel(
    data = centers,
    aes(x = UMAP_1, y = UMAP_2, label = grp),
    inherit.aes = FALSE,
    size = label.size,
    fontface = "bold",
    label.size = 0.28,
    fill = scales::alpha("white", 0.9),
    color = "black",
    segment.color = "grey30",
    box.padding = 0.35,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = max.overlaps,
    seed = 7
  )
}

saveRDS(seu, in_obj)
message("Updated full annotated object with region metadata: ", in_obj)

# UMAPs on all annotated cells (no exclusion)
if ("umap" %in% names(seu@reductions)) {
  ann_col <- if ("celltype_label" %in% colnames(seu@meta.data)) {
    "celltype_label"
  } else if ("celltype_fine" %in% colnames(seu@meta.data)) {
    "celltype_fine"
  } else {
    "seurat_clusters"
  }
  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  p_umap_annot_all <- DimPlot(
    seu, reduction = "umap", group.by = ann_col,
    label = !has_ggrepel, repel = !has_ggrepel, pt.size = 0.25
  ) +
    ggtitle("All annotated cells (no exclusion)")
  if (has_ggrepel) {
    p_umap_annot_all <- add_boxed_labels(p_umap_annot_all, seu, ann_col, label.size = 3.0)
  }
  save_plot(p_umap_annot_all, "umaps/03c_UMAP_all_cells_annotated_ALL", w = 15, h = 8)

  p_umap_region_all <- DimPlot(
    seu, reduction = "umap", group.by = "anatomical_region_display",
    cols = unname(region_colors[levels(seu$anatomical_region_display)]), pt.size = 0.25
  ) +
    ggtitle("All cells colored by anatomical region")
  save_plot(p_umap_region_all, "umaps/03c_UMAP_anatomical_region_ALL", w = 11, h = 8)

  if ("condition" %in% colnames(seu@meta.data)) {
    p_umap_region_split <- DimPlot(
      seu, reduction = "umap", group.by = "anatomical_region_display",
      split.by = "condition", cols = unname(region_colors[levels(seu$anatomical_region_display)]), pt.size = 0.25, ncol = 2
    ) +
      ggtitle("Anatomical region map split by condition (all cells)")
    save_plot(p_umap_region_split, "umaps/03c_UMAP_anatomical_region_split_by_condition_ALL", w = 13, h = 6)
  }
}

# Region composition tables
sample_col <- c("sample_id", "library_id", "orig.ident", "condition")
sample_col <- sample_col[sample_col %in% colnames(seu@meta.data)][1]

if (is.na(sample_col)) {
  stop("No sample-like column found. Expected one of: sample_id, library_id, orig.ident, condition")
}

df_region <- seu@meta.data %>%
  dplyr::mutate(
    sample_id = as.character(.data[[sample_col]]),
    condition = as.character(condition),
    Region = as.character(anatomical_region_display)
  ) %>%
  dplyr::filter(!is.na(Region), Region != "") %>%
  dplyr::count(sample_id, condition, Region, name = "n_cells") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

readr::write_csv(df_region, file.path(OUT$tables, "03c_region_composition_by_sample.csv"))

df_region_cond <- df_region %>%
  dplyr::group_by(condition, Region) %>%
  dplyr::summarise(n_cells = sum(n_cells), .groups = "drop") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
  dplyr::ungroup()

readr::write_csv(df_region_cond, file.path(OUT$tables, "03c_region_composition_by_condition.csv"))

# Apply the same factor level order as UMAP so colors are identical across figures
df_region$Region <- factor(df_region$Region, levels = region_levels)
df_region_cond$Region <- factor(df_region_cond$Region, levels = region_levels)

# Non-cortical/striatal burden per sample
non_tel_order <- df_region %>%
  dplyr::mutate(is_non_tel = !grepl("^cortical \\+ striatal$", Region, ignore.case = TRUE)) %>%
  dplyr::group_by(sample_id, condition) %>%
  dplyr::summarise(non_tel_frac = sum(frac[is_non_tel], na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(condition, dplyr::desc(non_tel_frac))
readr::write_csv(non_tel_order, file.path(OUT$tables, "03c_non_cortical_striatal_fraction_by_sample.csv"))

df_region$sample_id <- factor(df_region$sample_id, levels = non_tel_order$sample_id)

p_region_by_sample <- ggplot(df_region, aes(x = sample_id, y = frac, fill = Region)) +
  geom_col(width = 0.9) +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = region_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Region composition per sample/library (dissection bias check)",
    x = "Sample/Library",
    y = "Fraction of cells"
  ) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
save_plot(p_region_by_sample, "qc_plots/03c_DissectionBias_region_composition_by_sample", w = 13, h = 6)

p_region_by_condition <- ggplot(df_region_cond, aes(x = condition, y = frac, fill = Region)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = region_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Region composition by condition",
    x = NULL,
    y = "Fraction of cells"
  ) +
  theme_bw(base_size = 12)
save_plot(p_region_by_condition, "qc_plots/03c_DissectionBias_region_composition_by_condition", w = 8, h = 5)

# Alternate version with total n shown above each condition bar
cond_n <- df_region_cond %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(total_n = sum(n_cells), .groups = "drop")

p_region_by_condition_with_n <- p_region_by_condition +
  geom_text(
    data = cond_n,
    aes(x = condition, y = 1.03, label = paste0("n=", scales::comma(total_n))),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4
  ) +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  theme(plot.margin = margin(t = 12, r = 8, b = 8, l = 8))
save_plot(p_region_by_condition_with_n, "qc_plots/03c_DissectionBias_region_composition_by_condition_with_n", w = 8, h = 5)

# Main cortical + striatal scope object
if (any(seu$include_cortical_striatal_scope, na.rm = TRUE)) {
  seu_cortstr <- subset(seu, subset = include_cortical_striatal_scope)
  saveRDS(seu_cortstr, out_obj_cortstr)
  saveRDS(seu_cortstr, out_obj_tel_legacy)
  message("Saved cortical/striatal-focused object: ", out_obj_cortstr)
  message("Updated legacy compatibility object: ", out_obj_tel_legacy)

  if ("umap" %in% names(seu_cortstr@reductions)) {
    ann_col_scope <- if ("celltype_label" %in% colnames(seu_cortstr@meta.data)) {
      "celltype_label"
    } else if ("celltype_fine" %in% colnames(seu_cortstr@meta.data)) {
      "celltype_fine"
    } else {
      "seurat_clusters"
    }
    has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

    p_scope_annot <- DimPlot(
      seu_cortstr, reduction = "umap", group.by = ann_col_scope,
      label = !has_ggrepel, repel = !has_ggrepel, pt.size = 0.25
    ) +
      ggtitle("Cortical + Striatal analysis scope")
    if (has_ggrepel) {
      p_scope_annot <- add_boxed_labels(p_scope_annot, seu_cortstr, ann_col_scope, label.size = 3.0)
    }
    save_plot(p_scope_annot, "umaps/03c_UMAP_cortical_striatal_scope_annotated", w = 15, h = 8)

    p_scope_region <- DimPlot(
      seu_cortstr, reduction = "umap", group.by = "anatomical_region_display",
      cols = unname(region_colors[levels(seu$anatomical_region_display)]), pt.size = 0.25
    ) +
      ggtitle("Cortical + Striatal analysis scope colored by anatomical region")
    save_plot(p_scope_region, "umaps/03c_UMAP_cortical_striatal_scope_region", w = 11, h = 8)
  }
} else {
  stop("No cells marked include_cortical_striatal_scope == TRUE. Check region mapping file.")
}

# Sensitivity object: main scope + Hypothalamus_glutamatergic
if (any(seu$include_cortical_striatal_plus_hypoglut_sensitivity, na.rm = TRUE)) {
  seu_cortstr_hypo <- subset(seu, subset = include_cortical_striatal_plus_hypoglut_sensitivity)
  if ("include_in_analysis" %in% colnames(seu_cortstr_hypo@meta.data)) {
    hyp_idx <- seu_cortstr_hypo$include_hypothalamus_glutamatergic_sensitivity %in% TRUE
    seu_cortstr_hypo$include_in_analysis[hyp_idx] <- TRUE
  }
  saveRDS(seu_cortstr_hypo, out_obj_cortstr_hypo)
  message("Saved cortical/striatal + Hypothalamus_glutamatergic sensitivity object: ", out_obj_cortstr_hypo)

  if ("umap" %in% names(seu_cortstr_hypo@reductions)) {
    ann_col_sens <- if ("celltype_label" %in% colnames(seu_cortstr_hypo@meta.data)) {
      "celltype_label"
    } else if ("celltype_fine" %in% colnames(seu_cortstr_hypo@meta.data)) {
      "celltype_fine"
    } else {
      "seurat_clusters"
    }
    has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

    p_sens_annot <- DimPlot(
      seu_cortstr_hypo, reduction = "umap", group.by = ann_col_sens,
      label = !has_ggrepel, repel = !has_ggrepel, pt.size = 0.25
    ) +
      ggtitle("Cortical + Striatal scope + Hypothalamus_glutamatergic sensitivity set")
    if (has_ggrepel) {
      p_sens_annot <- add_boxed_labels(p_sens_annot, seu_cortstr_hypo, ann_col_sens, label.size = 3.0)
    }
    save_plot(p_sens_annot, "umaps/03c_UMAP_cortical_striatal_plus_hypoglut_sensitivity_annotated", w = 15, h = 8)

    p_sens_region <- DimPlot(
      seu_cortstr_hypo, reduction = "umap", group.by = "anatomical_region_display",
      cols = unname(region_colors[levels(seu$anatomical_region_display)]), pt.size = 0.25
    ) +
      ggtitle("Cortical + Striatal + Hypothalamus_glutamatergic sensitivity set by anatomical region")
    save_plot(p_sens_region, "umaps/03c_UMAP_cortical_striatal_plus_hypoglut_sensitivity_region", w = 11, h = 8)
  }
} else {
  stop("No cells marked include_cortical_striatal_plus_hypoglut_sensitivity == TRUE. Check annotation labels and region mapping.")
}

message("03c complete: dissection-bias outputs + cortical/striatal scope objects written.")
