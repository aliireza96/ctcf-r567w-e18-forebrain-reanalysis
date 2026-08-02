# ============================================================
# 03b_apply_annotation.R — apply manual annotation mapping
# ============================================================

source("scripts/00_setup.R")

in_obj  <- file.path(OUT$objects, "03_seu_markers_ready_for_annotation.rds")
in_map_region   <- file.path(OUT$tables, "03_annotation_master_table_with_regions.csv")
in_map_template <- file.path(OUT$tables, "03_annotation_template_fillme.csv")
in_map <- if (file.exists(in_map_region)) in_map_region else in_map_template
out_obj <- file.path(OUT$objects, "03_seu_annotated.rds")
out_obj_kept <- file.path(OUT$objects, "03_seu_annotated_KEPT.rds")
out_obj_tel <- file.path(OUT$objects, "03_seu_annotated_TEL.rds")

stopifnot(file.exists(in_obj))
stopifnot(file.exists(in_map))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

seu <- readRDS(in_obj)
mapping_raw <- readr::read_csv(in_map, show_col_types = FALSE)

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[y %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

# Build a unified mapping schema from either:
#   A) region-aware master table (preferred), or
#   B) legacy annotation template.
if (all(c("Cluster Number", "Detailed Label", "UMAP Label") %in% colnames(mapping_raw))) {
  cluster_num <- suppressWarnings(as.integer(mapping_raw[["Cluster Number"]]))

  # Allow a separate keep/exclude flag for downstream analysis (independent of telencephalon-only subset).
  include_analysis <- rep(NA, nrow(mapping_raw))
  if ("include_in_analysis" %in% colnames(mapping_raw)) {
    include_analysis <- as_bool(mapping_raw[["include_in_analysis"]])
  } else if ("Include in analysis?" %in% colnames(mapping_raw)) {
    include_analysis <- as_bool(mapping_raw[["Include in analysis?"]])
  } else if ("Recommended Set" %in% colnames(mapping_raw)) {
    rs <- tolower(trimws(as.character(mapping_raw[["Recommended Set"]])))
    include_analysis[grepl("^primary", rs) | grepl("^keep", rs)] <- TRUE
    include_analysis[grepl("^exclude", rs)] <- FALSE
  }

  mapping <- data.frame(
    seurat_cluster = as.character(cluster_num),
    celltype_fine = as.character(mapping_raw[["Detailed Label"]]),
    celltype_label = as.character(mapping_raw[["UMAP Label"]]),
    celltype_broad = if ("Major Class" %in% colnames(mapping_raw)) as.character(mapping_raw[["Major Class"]]) else as.character(mapping_raw[["Detailed Label"]]),
    anatomical_region = if ("Anatomical Region" %in% colnames(mapping_raw)) as.character(mapping_raw[["Anatomical Region"]]) else NA_character_,
    is_telencephalic = if ("Telencephalic?" %in% colnames(mapping_raw)) as_bool(mapping_raw[["Telencephalic?"]]) else if ("Telencephalic" %in% colnames(mapping_raw)) as_bool(mapping_raw[["Telencephalic"]]) else NA,
    include_telencephalon = if ("Include in telencephalon proportion analysis?" %in% colnames(mapping_raw)) as_bool(mapping_raw[["Include in telencephalon proportion analysis?"]]) else NA,
    include_in_analysis = include_analysis,
    stringsAsFactors = FALSE
  )
  mapping <- mapping[!is.na(cluster_num), , drop = FALSE]
  if (all(is.na(mapping$include_telencephalon))) {
    mapping$include_telencephalon <- mapping$is_telencephalic
  }
  mapping$include_in_analysis[is.na(mapping$include_in_analysis)] <- mapping$include_telencephalon[is.na(mapping$include_in_analysis)]
  mapping$include_in_analysis[is.na(mapping$include_in_analysis)] <- FALSE
} else {
  req_cols <- c("seurat_cluster", "celltype_fine", "celltype_broad")
  if (!all(req_cols %in% colnames(mapping_raw))) {
    stop("Mapping CSV must contain either region-master columns (Cluster Number, Detailed Label, UMAP Label) or template columns (seurat_cluster, celltype_fine, celltype_broad).")
  }
  mapping <- data.frame(
    seurat_cluster = as.character(mapping_raw$seurat_cluster),
    celltype_fine = as.character(mapping_raw$celltype_fine),
    celltype_broad = as.character(mapping_raw$celltype_broad),
    celltype_label = if ("celltype_label" %in% colnames(mapping_raw)) as.character(mapping_raw$celltype_label) else as.character(mapping_raw$celltype_fine),
    include_in_analysis = if ("include_in_analysis" %in% colnames(mapping_raw)) as_bool(mapping_raw$include_in_analysis) else TRUE,
    anatomical_region = if ("anatomical_region" %in% colnames(mapping_raw)) as.character(mapping_raw$anatomical_region) else NA_character_,
    is_telencephalic = if ("is_telencephalic" %in% colnames(mapping_raw)) as_bool(mapping_raw$is_telencephalic) else NA,
    include_telencephalon = if ("include_telencephalon" %in% colnames(mapping_raw)) as_bool(mapping_raw$include_telencephalon) else NA,
    stringsAsFactors = FALSE
  )
}

mapping <- mapping %>%
  dplyr::mutate(
    seurat_cluster = as.character(seurat_cluster),
    celltype_fine = as.character(celltype_fine),
    celltype_broad = as.character(celltype_broad),
    celltype_label = as.character(celltype_label),
    anatomical_region = as.character(anatomical_region),
    is_telencephalic = as_bool(is_telencephalic),
    include_telencephalon = as_bool(include_telencephalon),
    include_in_analysis = as_bool(include_in_analysis)
  ) %>%
  dplyr::filter(!is.na(seurat_cluster), seurat_cluster != "") %>%
  dplyr::distinct(seurat_cluster, .keep_all = TRUE)

# Keep microglia in the KEPT analysis set (requested).
microglia_idx <- grepl("microglia", tolower(trimws(mapping$celltype_fine))) |
  grepl("microglia", tolower(trimws(mapping$celltype_label)))
mapping$include_in_analysis[microglia_idx & is.na(mapping$include_in_analysis)] <- TRUE
mapping$include_in_analysis[microglia_idx & mapping$include_in_analysis == FALSE] <- TRUE

# Build broad groups exactly as requested:
# - merge upper-layer excitatory clusters
# - merge deep-layer excitatory clusters
# - merge SPN clusters
# - keep all other groups as their own label
upper_clusters <- c("3", "5")
deep_clusters  <- c("4", "6", "11")
spn_clusters   <- c("2", "12")

upper_labels <- c("it-l2/3", "it-l2/4")
deep_labels  <- c("corticofugal", "deep corticofugal", "deep cortical")
spn_labels   <- c("spn-d1", "spn-d2")

mapping <- mapping %>%
  dplyr::mutate(
    celltype_broad_original = celltype_broad,
    celltype_broad = dplyr::case_when(
      seurat_cluster %in% upper_clusters | tolower(trimws(celltype_label)) %in% upper_labels ~ "Upper layer EN",
      seurat_cluster %in% deep_clusters  | tolower(trimws(celltype_label)) %in% deep_labels  ~ "Deep layer EN",
      seurat_cluster %in% spn_clusters   | tolower(trimws(celltype_label)) %in% spn_labels   ~ "SPNs",
      TRUE ~ celltype_label
    )
  )

# If telencephalon flag is missing, fall back to include_in_analysis.
if (all(is.na(mapping$include_telencephalon))) {
  mapping$include_telencephalon <- mapping$include_in_analysis
}
mapping$include_telencephalon[is.na(mapping$include_telencephalon)] <- FALSE
mapping$include_in_analysis[is.na(mapping$include_in_analysis)] <- mapping$include_telencephalon[is.na(mapping$include_in_analysis)]

# Check all clusters are covered
clusters <- levels(seu$seurat_clusters)
missing <- setdiff(clusters, mapping$seurat_cluster)
if (length(missing) > 0) {
  stop("Mapping is missing clusters: ", paste(missing, collapse = ", "))
}

# Check for NA/blank annotations
if (any(is.na(mapping$celltype_fine)) || any(mapping$celltype_fine == "")) {
  stop("celltype_fine has missing/blank entries. Fill all rows in the mapping CSV.")
}
if (any(is.na(mapping$celltype_broad)) || any(mapping$celltype_broad == "")) {
  stop("celltype_broad has missing/blank entries. Fill all rows in the mapping CSV.")
}
if (any(is.na(mapping$celltype_label)) || any(mapping$celltype_label == "")) {
  stop("celltype_label has missing/blank entries. Fill all rows (or remove the column to default to celltype_fine).")
}
if (any(is.na(mapping$include_in_analysis))) {
  stop("include_in_analysis has missing/blank entries. Use TRUE/FALSE for every cluster.")
}

readr::write_csv(mapping, file.path(OUT$tables, "03b_cluster_mapping_expanded.csv"))

# ---- Apply mapping (no plyr) ----
map_fine  <- setNames(mapping$celltype_fine,  mapping$seurat_cluster)
map_broad <- setNames(mapping$celltype_broad, mapping$seurat_cluster)
map_label <- setNames(mapping$celltype_label, mapping$seurat_cluster)
map_keep  <- setNames(mapping$include_in_analysis, mapping$seurat_cluster)
map_region <- setNames(mapping$anatomical_region, mapping$seurat_cluster)
map_tel <- setNames(mapping$is_telencephalic, mapping$seurat_cluster)
map_tel_keep <- setNames(mapping$include_telencephalon, mapping$seurat_cluster)

seu$celltype_fine  <- unname(map_fine[as.character(seu$seurat_clusters)])
seu$celltype_broad <- unname(map_broad[as.character(seu$seurat_clusters)])
seu$celltype_label <- unname(map_label[as.character(seu$seurat_clusters)])
seu$include_in_analysis <- unname(map_keep[as.character(seu$seurat_clusters)])
seu$anatomical_region <- unname(map_region[as.character(seu$seurat_clusters)])
seu$is_telencephalic <- unname(map_tel[as.character(seu$seurat_clusters)])
seu$include_telencephalon <- unname(map_tel_keep[as.character(seu$seurat_clusters)])

# Factor levels (stable ordering)
fine_levels  <- unique(mapping$celltype_fine)
broad_levels <- unique(mapping$celltype_broad)
label_levels <- unique(mapping$celltype_label)

seu$celltype_fine  <- factor(seu$celltype_fine,  levels = fine_levels)
seu$celltype_broad <- factor(seu$celltype_broad, levels = broad_levels)
seu$celltype_label <- factor(seu$celltype_label, levels = label_levels)
if (any(!is.na(seu$anatomical_region))) {
  reg_levels <- unique(mapping$anatomical_region[!is.na(mapping$anatomical_region) & mapping$anatomical_region != ""])
  seu$anatomical_region <- factor(seu$anatomical_region, levels = reg_levels)
}

# Save annotated object (full) and a kept-only version for downstream analysis
saveRDS(seu, out_obj)
message("Saved annotated object (full): ", out_obj)

# ---- Annotation summary tables ----
ann_summary <- seu@meta.data %>%
  count(condition, celltype_broad, celltype_fine, include_in_analysis, name = "n_cells") %>%
  group_by(condition) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  ungroup()

readr::write_csv(ann_summary, file.path(OUT$tables, "03b_annotation_counts_by_condition.csv"))

# Convenience: summary restricted to clusters marked include_in_analysis == TRUE
ann_summary_kept <- ann_summary %>%
  filter(include_in_analysis) %>%
  group_by(condition) %>%
  mutate(frac = n_cells / sum(n_cells)) %>%
  ungroup()
readr::write_csv(ann_summary_kept, file.path(OUT$tables, "03b_annotation_counts_by_condition_KEPT.csv"))

# ---- Cluster presence by condition (all clusters) ----
cluster_presence <- seu@meta.data %>%
  dplyr::count(seurat_clusters, condition, name = "n_cells") %>%
  tidyr::pivot_wider(names_from = condition, values_from = n_cells, values_fill = 0)
if (!("wt" %in% colnames(cluster_presence))) cluster_presence$wt <- 0L
if (!("mut" %in% colnames(cluster_presence))) cluster_presence$mut <- 0L
cluster_presence <- cluster_presence %>%
  dplyr::mutate(
    total = wt + mut,
    wt_frac = dplyr::if_else(total > 0, wt / total, NA_real_),
    mut_frac = dplyr::if_else(total > 0, mut / total, NA_real_),
    exclusivity = dplyr::case_when(
      wt == 0 & mut > 0 ~ "mut_only",
      mut == 0 & wt > 0 ~ "wt_only",
      TRUE ~ "both"
    )
  ) %>%
  dplyr::arrange(as.numeric(as.character(seurat_clusters)))
readr::write_csv(cluster_presence, file.path(OUT$tables, "03b_cluster_presence_by_condition.csv"))

# ---- Region composition summaries (dissection bias diagnostics; all clusters) ----
if ("anatomical_region" %in% colnames(seu@meta.data)) {
  sample_col <- c("sample_id", "library_id", "orig.ident", "condition")
  sample_col <- sample_col[sample_col %in% colnames(seu@meta.data)][1]
  
  if (!is.na(sample_col)) {
    df_region <- seu@meta.data %>%
      dplyr::mutate(
        sample_id = as.character(.data[[sample_col]]),
        condition = as.character(condition),
        Region = as.character(anatomical_region)
      ) %>%
      dplyr::filter(!is.na(Region), Region != "") %>%
      dplyr::count(sample_id, condition, Region, name = "n_cells") %>%
      dplyr::group_by(sample_id) %>%
      dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
      dplyr::ungroup()
    readr::write_csv(df_region, file.path(OUT$tables, "03b_region_composition_by_sample.csv"))
    
    df_region_cond <- seu@meta.data %>%
      dplyr::mutate(
        condition = as.character(condition),
        Region = as.character(anatomical_region)
      ) %>%
      dplyr::filter(!is.na(Region), Region != "") %>%
      dplyr::count(condition, Region, name = "n_cells") %>%
      dplyr::group_by(condition) %>%
      dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
      dplyr::ungroup()
    readr::write_csv(df_region_cond, file.path(OUT$tables, "03b_region_composition_by_condition.csv"))
  }
}

# ---- Plot helpers (publication-ish) ----
# Distinct qualitative palette (no manual colors needed; derived from colorspace)
get_palette <- function(n) {
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    # Fallback: ggplot hue palette
    return(scales::hue_pal()(n))
  }
  colorspace::qualitative_hcl(n, palette = "Dark 3")
}

# Add boxed labels at cluster centroids (repelled) to a DimPlot
add_boxed_labels <- function(p, seu, group.by, label.size = 3.2, max.overlaps = 200) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) return(p)
  
  emb <- as.data.frame(Seurat::Embeddings(seu, reduction = "umap"))
  # Normalize UMAP column names for robustness across Seurat versions
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(emb))) {
    umap_cols <- grep("^UMAP", colnames(emb), value = TRUE, ignore.case = TRUE)
    if (length(umap_cols) >= 2) {
      emb <- emb[, umap_cols[1:2], drop = FALSE]
      colnames(emb) <- c("UMAP_1", "UMAP_2")
    } else if (ncol(emb) >= 2) {
      colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
    } else {
      stop("UMAP embeddings not found or have unexpected shape. Check that UMAP was run.")
    }
  }
  emb$grp <- seu[[group.by, drop = TRUE]]
  
  centers <- emb %>%
    group_by(grp) %>%
    summarise(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(grp = as.character(grp))
  
  p +
    ggrepel::geom_label_repel(
      data = centers,
      aes(x = UMAP_1, y = UMAP_2, label = grp),
      size = label.size,
      label.size = 0.25,
      fill = "white",
      color = "black",
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0,
      max.overlaps = max.overlaps,
      seed = 1
    )
}

# Add boxed labels on split UMAP panels (e.g., split.by = condition)
add_boxed_labels_split <- function(p, seu, group.by, split.by, label.size = 3.4, max.overlaps = 200) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) return(p)

  emb <- as.data.frame(Seurat::Embeddings(seu, reduction = "umap"))
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(emb))) {
    umap_cols <- grep("^UMAP", colnames(emb), value = TRUE, ignore.case = TRUE)
    if (length(umap_cols) >= 2) {
      emb <- emb[, umap_cols[1:2], drop = FALSE]
      colnames(emb) <- c("UMAP_1", "UMAP_2")
    } else if (ncol(emb) >= 2) {
      colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
    } else {
      stop("UMAP embeddings not found or have unexpected shape. Check that UMAP was run.")
    }
  }

  emb$grp <- as.character(seu[[group.by, drop = TRUE]])
  emb$split_val <- as.character(seu[[split.by, drop = TRUE]])

  centers <- emb %>%
    dplyr::group_by(split_val, grp) %>%
    dplyr::summarise(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )
  centers[[split.by]] <- factor(centers$split_val, levels = levels(factor(seu[[split.by, drop = TRUE]])))
  centers$split_val <- NULL

  p +
    ggrepel::geom_label_repel(
      data = centers,
      aes(x = UMAP_1, y = UMAP_2, label = grp),
      inherit.aes = FALSE,
      size = label.size,
      fontface = "bold",
      label.size = 0.28,
      label.padding = grid::unit(0.15, "lines"),
      label.r = grid::unit(0.08, "lines"),
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

# Theme tuned for figures
theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.text  = element_text(size = 10)
    )
}

# ---- Plot: all clusters for dissection-bias defense ----
p_cluster_num_all <- Seurat::DimPlot(
  seu,
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) +
  ggtitle("All clusters (numbered)") +
  theme_pub() +
  theme(legend.position = "none")
save_plot(p_cluster_num_all, "umaps/UMAP_clusters_numbered_ALL", w = 9, h = 7)

p_cond_overlay_all <- Seurat::DimPlot(
  seu,
  group.by = "condition",
  label = FALSE,
  repel = FALSE,
  pt.size = 0.25
) +
  ggtitle("WT vs MUT overlay (all clusters)") +
  theme_pub()
save_plot(p_cond_overlay_all, "umaps/UMAP_condition_overlay_ALL", w = 9, h = 7)

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

p_split_clusters_all <- Seurat::DimPlot(
  seu,
  group.by = "seurat_clusters",
  split.by = "condition",
  label = !has_ggrepel,
  repel = !has_ggrepel,
  pt.size = 0.2,
  ncol = 2
) +
  ggtitle("All clusters split by condition") +
  theme_pub() +
  theme(legend.position = "none")

if (has_ggrepel) {
  p_split_clusters_all <- add_boxed_labels_split(
    p_split_clusters_all,
    seu = seu,
    group.by = "seurat_clusters",
    split.by = "condition",
    label.size = 3.1
  )
}

save_plot(p_split_clusters_all, "umaps/UMAP_split_by_condition_clusters_ALL", w = 14, h = 6)

save_plot(p_cluster_num_all + p_cond_overlay_all, "umaps/UMAP_cluster_numbers_and_condition_overlay_ALL", w = 14, h = 6)

# ---- Plot: dissection bias by region ----
if ("anatomical_region" %in% colnames(seu@meta.data)) {
  sample_col <- c("sample_id", "library_id", "orig.ident", "condition")
  sample_col <- sample_col[sample_col %in% colnames(seu@meta.data)][1]
  
  if (!is.na(sample_col)) {
    df_region <- seu@meta.data %>%
      dplyr::mutate(
        sample_id = as.character(.data[[sample_col]]),
        condition = as.character(condition),
        Region = as.character(anatomical_region)
      ) %>%
      dplyr::filter(!is.na(Region), Region != "") %>%
      dplyr::count(sample_id, condition, Region, name = "n_cells") %>%
      dplyr::group_by(sample_id) %>%
      dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
      dplyr::ungroup()
    
    non_tel_order <- df_region %>%
      dplyr::mutate(is_non_tel = !grepl("^telencephalon$", Region, ignore.case = TRUE)) %>%
      dplyr::group_by(sample_id, condition) %>%
      dplyr::summarise(non_tel_frac = sum(frac[is_non_tel], na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(condition, dplyr::desc(non_tel_frac))
    
    df_region$sample_id <- factor(df_region$sample_id, levels = non_tel_order$sample_id)
    
    p_region_by_sample <- ggplot(df_region, aes(x = sample_id, y = frac, fill = Region)) +
      geom_col(width = 0.9) +
      facet_grid(. ~ condition, scales = "free_x", space = "free_x") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        title = "Region composition per sample/library (dissection bias check)",
        x = "Sample/Library",
        y = "Fraction of cells"
      ) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    save_plot(p_region_by_sample, "qc_plots/DissectionBias_region_composition_by_sample", w = 13, h = 6)
    
    df_region_cond <- df_region %>%
      dplyr::group_by(condition, Region) %>%
      dplyr::summarise(n_cells = sum(n_cells), .groups = "drop") %>%
      dplyr::group_by(condition) %>%
      dplyr::mutate(frac = n_cells / sum(n_cells)) %>%
      dplyr::ungroup()
    
    p_region_by_condition <- ggplot(df_region_cond, aes(x = condition, y = frac, fill = Region)) +
      geom_col(width = 0.8) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        title = "Region composition by condition",
        x = NULL,
        y = "Fraction of cells"
      ) +
      theme_bw(base_size = 12)
    save_plot(p_region_by_condition, "qc_plots/DissectionBias_region_composition_by_condition", w = 8, h = 5)
  }
}

# Create "kept-only" object for plotting and downstream analysis
seu_kept <- subset(seu, subset = include_in_analysis)
saveRDS(seu_kept, out_obj_kept)
message("Saved kept-only object: ", out_obj_kept)

# Save explicit telencephalon-focused object for downstream analyses.
if ("include_telencephalon" %in% colnames(seu@meta.data) && any(seu$include_telencephalon, na.rm = TRUE)) {
  seu_tel <- subset(seu, subset = include_telencephalon)
  saveRDS(seu_tel, out_obj_tel)
  message("Saved telencephalon object: ", out_obj_tel)
}

# Drop unused factor levels so excluded categories cannot appear in legends/axes.
seu_kept$celltype_fine  <- droplevels(seu_kept$celltype_fine)
seu_kept$celltype_broad <- droplevels(seu_kept$celltype_broad)
seu_kept$celltype_label <- droplevels(seu_kept$celltype_label)

# ---- Plot: UMAP (labels, kept-only, boxed labels) ----
p_fine <- Seurat::DimPlot(
  seu_kept,
  group.by = "celltype_label",
  label = FALSE,
  repel = FALSE,
  pt.size = 0.25
) +
  ggtitle("E18.5 annotated cell types (labels) — kept clusters") +
  scale_color_manual(values = get_palette(length(levels(seu_kept$celltype_label)))) +
  theme_pub()

p_fine <- add_boxed_labels(p_fine, seu_kept, group.by = "celltype_label", label.size = 3.0)

save_plot(p_fine, "umaps/UMAP_celltype_fine_labeled_KEPT", w = 9, h = 7)

# ---- Plot: UMAP (broad, kept-only, boxed labels) ----
p_broad <- Seurat::DimPlot(
  seu_kept,
  group.by = "celltype_broad",
  label = FALSE,
  repel = FALSE,
  pt.size = 0.25
) +
  ggtitle("E18.5 annotated classes (broad) — kept clusters") +
  scale_color_manual(values = get_palette(length(levels(seu_kept$celltype_broad)))) +
  theme_pub()

p_broad <- add_boxed_labels(p_broad, seu_kept, group.by = "celltype_broad", label.size = 3.3)

save_plot(p_broad, "umaps/UMAP_celltype_broad_labeled_KEPT", w = 8, h = 6)

# ---- Plot: Split by condition (fine label, kept-only) ----
# Use short plotting label if provided (celltype_label)
p_split <- Seurat::DimPlot(
  seu_kept,
  group.by = "celltype_label",
  split.by = "condition",
  pt.size = 0.2,
  ncol = 2
) +
  ggtitle("Annotated cell types (kept) split by condition") +
  scale_color_manual(values = get_palette(length(levels(seu_kept$celltype_label)))) +
  theme_pub() +
  theme(legend.position = "none")

save_plot(p_split, "umaps/UMAP_celltype_fine_split_by_condition_KEPT", w = 14, h = 6)

# ---- Plot: Left labeled UMAP, right WT/MUT overlay (kept-only) ----
p_left <- Seurat::DimPlot(
  seu_kept,
  group.by = "celltype_label",
  label = FALSE,
  repel = FALSE,
  pt.size = 0.25
) +
  ggtitle("Celltype labels (kept)") +
  scale_color_manual(values = get_palette(length(levels(seu_kept$celltype_label)))) +
  theme_pub()

p_left <- add_boxed_labels(p_left, seu_kept, group.by = "celltype_label", label.size = 3.0)

p_right <- Seurat::DimPlot(
  seu_kept,
  group.by = "condition",
  label = FALSE,
  repel = FALSE,
  pt.size = 0.25
) +
  ggtitle("WT vs MUT (overlay)") +
  theme_pub()

save_plot(p_left + p_right, "umaps/UMAP_labels_vs_condition_overlay_KEPT", w = 14, h = 6)

# ---- Improved marker DotPlots (less broad; interpretable) ----
# Instead of one giant "broad" list, use a compact panel of canonical markers per broad class.
marker_panels <- list(
  "Cycling Radial Glia" = c("Mki67", "Top2a", "Pax6", "Fabp7", "Slc1a3"),
  "Excitatory"          = c("Satb2", "Neurod2", "Neurod6", "Tbr1", "Fezf2", "Tle4", "Bcl11b"),

  # Ventral GABA subtypes
  "MGE-derived Cells"    = c("Gad1", "Gad2", "Lhx6", "Pvalb", "Sst"),
  "CGE-derived Cells"    = c("Gad1", "Gad2", "Adarb2", "Nr2f2", "Reln", "Vip", "Lhx6"),
  "LGE-derived Cells"    = c("Gad1", "Gad2", "Isl1", "Ebf1", "Gpr88", "Oprm1", "Drd2", "Penk", "Adora2a", "Tshz1", "Dlx1", "Sp9"),
  
  # Cluster 0 (Hypothalamic glutamatergic, Trhde+, imprinted-gene enriched)
  "Cl0_Hypothalamic glutamatergic" = c(
    "Slc17a6", "Slc6a7", "Trhde",
    "Magel2", "Peg3", "Nnat", "Meg3", "Rian",
    "Scg2", "Vgf", "Chga", "Pcsk1n"
  ),
  "Cl0_RuleOut_Thal_GABA_Ctx" = c(
    "Tcf7l2", "Gbx2",          # thalamus anchors
    "Gad1", "Gad2", "Slc6a1",  # GABAergic
    "Slc17a7"                  # cortical VGLUT1 excitatory
  ),
  
  # Cluster 9 (Cerebellar granule lineage, Barhl1+, likely differentiating/postmitotic)
  "Cl9_Cb_GranuleLineage" = c(
    "Barhl1", "Pax6",
    "Zic1", "Zic2", "Zic4", "Zic5",
    "Cbln1",
    "Neurod1", "Neurod2",
    "Dcc", "Unc5c", "Dcx"
  ),
  "Cl9_StageCheck_GCP_vs_Postmitotic" = c(
    "Atoh1",                   # proliferative GCP driver
    "Mki67", "Top2a", "Pcna",  # cycling
    "Mcm2", "Mcm5"
  ),
  "Cl9_RuleOut_CR_Meninges_Thal" = c(
    "Trp73", "Lhx1", "Lhx5", "Calb2", "Ebf2",  # Cajal-Retzius
    "Foxc1", "Col1a1", "Dcn",                  # meninges
    "Tcf7l2", "Gbx2"                           # thalamus/diencephalon
  ),
  
  "OPC"           = c("Pdgfra", "Cspg4", "Olig2", "Sox10", "Gpr17"),
  "Immature Astrocytes"  = c("Apoe", "Aldh1l1", "Slc1a3", "Nfia", "Sox9"),
  "Microglia"           = c("Cx3cr1", "C1qc", "Trem2", "P2ry13", "Fcrls"),
  "Endothelial cells"       = c("Cldn5", "Robo4", "Pcdh12", "Cd93", "Aplnr")
)

# Helper: save one dotplot per panel (readable in a paper)
plot_dot_panel <- function(seu_obj, features, title, group.by = "celltype_broad") {
  feats <- unique(features)
  feats <- feats[feats %in% rownames(seu_obj)]
  if (length(feats) < 4) return(NULL)
  
  p <- Seurat::DotPlot(seu_obj, features = feats, group.by = group.by, scale = FALSE) +
    ggtitle(title) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      axis.title = element_blank(),
      panel.grid = element_blank()
    ) +
    Seurat::RotatedAxis()
  
  return(p)
}

for (nm in names(marker_panels)) {
  p_dot <- plot_dot_panel(
    seu_obj = seu_kept,
    features = marker_panels[[nm]],
    title = paste0("Marker sanity-check: ", nm),
    group.by = "celltype_broad"
  )
  if (!is.null(p_dot)) {
    save_plot(p_dot, paste0("dotplots/Dotplot_markers_", gsub("[^A-Za-z0-9]+", "_", nm)), w = 12, h = 5.5)
  } else {
    message("DotPlot skipped for panel '", nm, "': too few markers found in dataset.")
  }
}

# Optional: a single "focused ventral interneuron" dotplot (useful for your CGE/LGE boundary question)
ventral_markers <- c("Nr2f2", "Isl1", "Six3", "Ecel1", "Reln", "Vip", "Lhx6", "Sst", "Pvalb", "Tac1", "Gad1", "Gad2")
ventral_markers <- ventral_markers[ventral_markers %in% rownames(seu_kept)]
if (length(ventral_markers) >= 6) {
  p_vent <- Seurat::DotPlot(seu_kept, features = ventral_markers, group.by = "celltype_fine", scale = FALSE) +
    ggtitle("Ventral / interneuron-focused markers (fine types)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
          panel.grid = element_blank()) +
    Seurat::RotatedAxis()
  save_plot(p_vent, "dotplots/Dotplot_ventral_markers_by_finetype", w = 13, h = 6)
}

message("03b complete: annotation applied + publication-quality plots/tables saved.")
