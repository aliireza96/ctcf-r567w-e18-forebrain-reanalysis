#!/usr/bin/env Rscript

# ============================================================
# 09f_combined_lineage_umaps.R
# Requested visualization panels:
# 1) E18 query-only UMAP with combined mapped pseudotime (ventral+dorsal)
#    + split WT/MUT panels
# 2) Per-lineage two-panel figure:
#    left = reference/root UMAP colored by pseudotime_ref (root highlighted)
#    right = query UMAP split by genotype colored by mapped pseudotime
# 3) Export exact roots used for each lineage
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(grid)
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

msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

save_pdf <- function(plot_obj, path_no_ext, w = 8.5, h = 6.0) {
  ggsave(paste0(path_no_ext, ".pdf"), plot_obj, width = w, height = h, useDingbats = FALSE)
}

save_side_by_side <- function(p_left, p_right, path_no_ext, w = 13, h = 5.8) {
  grDevices::pdf(paste0(path_no_ext, ".pdf"), width = w, height = h, useDingbats = FALSE)
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 1, ncol = 2)
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

canon_genotype <- function(x, wt_values, mut_values) {
  x <- as.character(x)
  x_low <- tolower(trimws(x))
  wt_low <- tolower(trimws(as.character(wt_values)))
  mut_low <- tolower(trimws(as.character(mut_values)))
  out <- rep(NA_character_, length(x))
  out[x_low %in% wt_low] <- "WT"
  out[x_low %in% mut_low] <- "MUT"
  out[is.na(out) & grepl("(^wt$|wild|control|ctrl)", x_low)] <- "WT"
  out[is.na(out) & grepl("(mut|r567w|homo|hom|ko|case)", x_low)] <- "MUT"
  out
}

resolve_genotype_col <- function(md, preferred = NULL) {
  cand <- unique(c(preferred, "genotype", "genotype2", "condition"))
  cand <- cand[!is.na(cand) & nzchar(cand)]
  hit <- cand[cand %in% colnames(md)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

extract_umap_dt <- function(seu, preferred = c("umap", "ref.umap"), tag = "object") {
  red_hit <- preferred[preferred %in% names(seu@reductions)]
  if (length(red_hit) == 0) {
    stop("[", tag, "] none of reductions found: ", paste(preferred, collapse = ", "))
  }
  red <- red_hit[1]
  um <- Embeddings(seu, red)
  data.table(cell_id = rownames(um), UMAP_1 = um[, 1], UMAP_2 = um[, 2], reduction_used = red)
}

scale01 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[2] <= rng[1]) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

read_mapping_tables <- function(base_dir) {
  map_dir <- file.path(base_dir, "mapping", "tables")
  files <- list.files(map_dir, pattern = "^09b_.*_mapped_table\\.csv$", full.names = TRUE)
  if (length(files) == 0) stop("No 09b mapped tables found under: ", map_dir)

  out <- vector("list", length(files))
  for (i in seq_along(files)) {
    f <- files[i]
    d <- fread(f)
    if (!"lineage_name" %in% colnames(d)) {
      nm <- sub("^09b_(.*)_mapped_table\\.csv$", "\\1", basename(f))
      d[, lineage_name := nm]
    }
    req <- c("cell_id", "lineage_name", "pseudotime_transfer")
    miss <- req[!req %in% colnames(d)]
    if (length(miss) > 0) {
      stop("Mapped table missing required columns in ", basename(f), ": ", paste(miss, collapse = ", "))
    }
    if (!"pseudotime_transfer_score" %in% colnames(d)) d[, pseudotime_transfer_score := NA_real_]
    d[, pseudotime_transfer := suppressWarnings(as.numeric(pseudotime_transfer))]
    d[, pseudotime_transfer_score := suppressWarnings(as.numeric(pseudotime_transfer_score))]
    d <- d[!is.na(pseudotime_transfer)]
    d[, pseudotime_scaled_lineage := scale01(pseudotime_transfer), by = lineage_name]
    out[[i]] <- d[, .(cell_id, lineage_name, pseudotime_transfer, pseudotime_transfer_score, pseudotime_scaled_lineage)]
  }

  map_all <- rbindlist(out, fill = TRUE)
  map_best <- copy(map_all)
  map_best[is.na(pseudotime_transfer_score), pseudotime_transfer_score := -Inf]
  setorder(map_best, cell_id, -pseudotime_transfer_score)
  map_best <- map_best[, .SD[1], by = cell_id]
  map_best[pseudotime_transfer_score == -Inf, pseudotime_transfer_score := NA_real_]

  list(all = map_all, best = map_best)
}

read_roots_summary <- function(base_dir) {
  f <- file.path(base_dir, "refs", "tables", "09a_reference_lineage_summary.csv")
  if (!file.exists(f)) stop("Missing root summary table: ", f)
  rs <- fread(f)
  req <- c("lineage_name", "root_label", "start_cluster", "filter_label_col")
  miss <- req[!req %in% colnames(rs)]
  if (length(miss) > 0) stop("Root summary missing columns: ", paste(miss, collapse = ", "))
  rs
}

make_query_combined_df <- function(query_obj, map_best) {
  um <- extract_umap_dt(query_obj, preferred = c("umap"), tag = "query_e18")
  md0 <- query_obj@meta.data
  # Guard against pre-existing cell_id column in meta.data
  if ("cell_id" %in% colnames(md0)) md0$cell_id <- NULL
  md <- as.data.table(md0, keep.rownames = "cell_id")
  gcol <- resolve_genotype_col(md, preferred = CFG$cols$genotype_col)
  if (is.null(gcol)) md[, genotype2 := NA_character_] else md[, genotype2 := get(gcol)]
  md[, genotype2 := canon_genotype(genotype2, CFG$genotype$wt_values, CFG$genotype$mut_values)]

  out <- merge(um, md, by = "cell_id", all.x = TRUE, all.y = FALSE)
  out <- merge(out, map_best, by = "cell_id", all.x = TRUE, all.y = FALSE)
  out
}

plot_query_combined <- function(dt, out_plots) {
  p_all <- ggplot(dt, aes(UMAP_1, UMAP_2, color = pseudotime_scaled_lineage)) +
    geom_point(size = 0.22, alpha = 0.9) +
    scale_color_viridis_c(option = "plasma", na.value = "grey85", limits = c(0, 1), name = "Scaled PT") +
    theme_classic(base_size = 12) +
    labs(
      title = "E18.5 query UMAP: combined mapped pseudotime (ventral + dorsal)",
      subtitle = "All query cells shown; unmapped cells in grey",
      x = "UMAP_1", y = "UMAP_2"
    )
  save_pdf(p_all, file.path(out_plots, "09f_queryE18_combined_pseudotime_allcells"), w = 8.2, h = 6.2)

  q2 <- dt[!is.na(genotype2)]
  p_split <- NULL
  if (nrow(q2) > 0 && uniqueN(q2$genotype2) >= 2) {
    p_split <- ggplot(q2, aes(UMAP_1, UMAP_2, color = pseudotime_scaled_lineage)) +
      geom_point(size = 0.22, alpha = 0.9) +
      scale_color_viridis_c(option = "plasma", na.value = "grey85", limits = c(0, 1), name = "Scaled PT") +
      facet_wrap(~ genotype2) +
      theme_classic(base_size = 12) +
      labs(
        title = "E18.5 query UMAP: combined mapped pseudotime by genotype",
        subtitle = "WT vs MUT on the same E18 manifold",
        x = "UMAP_1", y = "UMAP_2"
      )
    save_pdf(p_split, file.path(out_plots, "09f_queryE18_combined_pseudotime_split_genotype"), w = 10.8, h = 5.4)
  }

  # Add original query cell-type panel and save combined figure requested by user.
  ct_col <- if (CFG$cols$broad_col %in% colnames(dt)) CFG$cols$broad_col else if ("celltype_label" %in% colnames(dt)) "celltype_label" else NA_character_
  if (!is.na(ct_col)) {
    p_ct <- ggplot(dt, aes(UMAP_1, UMAP_2, color = .data[[ct_col]])) +
      geom_point(size = 0.22, alpha = 0.9) +
      theme_classic(base_size = 12) +
      labs(
        title = "E18.5 query UMAP: original cell types",
        subtitle = paste0("Color by ", ct_col),
        x = "UMAP_1", y = "UMAP_2", color = "Cell type"
      )
    save_pdf(p_ct, file.path(out_plots, "09f_queryE18_combined_celltype"), w = 8.8, h = 6.2)
    if (!is.null(p_split)) {
      save_side_by_side(
        p_left = p_split + theme(legend.position = "right"),
        p_right = p_ct + theme(legend.position = "right"),
        path_no_ext = file.path(out_plots, "09f_queryE18_combined_pseudotime_split_genotype_celltype"),
        w = 15.8, h = 5.6
      )
    }
  }
}

make_reference_root_df <- function(base_dir, lineage_name, root_label, filter_label_col) {
  obj_path <- file.path(base_dir, "refs", "objects", paste0("09a_", lineage_name, "_reference_lineage.rds"))
  if (!file.exists(obj_path)) stop("Missing reference lineage object: ", obj_path)
  ref <- readRDS(obj_path)

  um <- extract_umap_dt(ref, preferred = c("umap", "ref.umap"), tag = paste0("ref_", lineage_name))
  md0 <- ref@meta.data
  # Guard against pre-existing cell_id column in meta.data
  if ("cell_id" %in% colnames(md0)) md0$cell_id <- NULL
  md <- as.data.table(md0, keep.rownames = "cell_id")
  dt <- merge(um, md, by = "cell_id", all.x = TRUE, all.y = FALSE)

  if (!"pseudotime_ref" %in% colnames(dt)) stop("[", lineage_name, "] missing pseudotime_ref in reference object")
  dt[, pseudotime_ref := suppressWarnings(as.numeric(pseudotime_ref))]
  dt[, pseudotime_ref_scaled := scale01(pseudotime_ref)]

  dt[, root_flag := FALSE]
  if (!is.na(filter_label_col) && nzchar(filter_label_col) && filter_label_col %in% colnames(dt)) {
    dt[, root_flag := as.character(get(filter_label_col)) == as.character(root_label)]
  }

  dt
}

make_query_lineage_df <- function(query_dt, map_all, lineage_name) {
  lin <- lineage_name
  m <- map_all[lineage_name == lin, .(cell_id, pseudotime_transfer, pseudotime_transfer_score)]
  if (nrow(m) == 0) return(NULL)
  m[, pseudotime_transfer_scaled := scale01(pseudotime_transfer)]

  q <- merge(
    query_dt[, .(cell_id, UMAP_1, UMAP_2, genotype2)],
    m,
    by = "cell_id",
    all = FALSE
  )
  q
}

plot_lineage_root_and_query <- function(ref_dt, qry_dt, lineage_name, root_label, out_plots) {
  p_root <- ggplot(ref_dt, aes(UMAP_1, UMAP_2, color = pseudotime_ref_scaled)) +
    geom_point(size = 0.22, alpha = 0.9) +
    geom_point(
      data = ref_dt[root_flag == TRUE],
      aes(UMAP_1, UMAP_2),
      inherit.aes = FALSE,
      shape = 21, fill = "white", color = "black", stroke = 0.35, size = 0.9
    ) +
    scale_color_viridis_c(option = "plasma", limits = c(0, 1), name = "Ref PT") +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Reference pseudotime: ", lineage_name),
      subtitle = paste0("Root label: ", root_label, " (outlined)"),
      x = "UMAP_1", y = "UMAP_2"
    )

  q2 <- qry_dt[!is.na(genotype2)]
  if (nrow(q2) > 0 && uniqueN(q2$genotype2) >= 2) {
    p_q <- ggplot(q2, aes(UMAP_1, UMAP_2, color = pseudotime_transfer_scaled)) +
      geom_point(size = 0.22, alpha = 0.9) +
      scale_color_viridis_c(option = "plasma", limits = c(0, 1), name = "Mapped PT") +
      facet_wrap(~ genotype2) +
      theme_classic(base_size = 12) +
      labs(
        title = paste0("E18 query mapped pseudotime: ", lineage_name),
        subtitle = "Split by genotype",
        x = "UMAP_1", y = "UMAP_2"
      )
  } else {
    p_q <- ggplot(qry_dt, aes(UMAP_1, UMAP_2, color = pseudotime_transfer_scaled)) +
      geom_point(size = 0.22, alpha = 0.9) +
      scale_color_viridis_c(option = "plasma", limits = c(0, 1), name = "Mapped PT") +
      theme_classic(base_size = 12) +
      labs(
        title = paste0("E18 query mapped pseudotime: ", lineage_name),
        subtitle = "Genotype split unavailable",
        x = "UMAP_1", y = "UMAP_2"
      )
  }

  save_pdf(p_root, file.path(out_plots, paste0("09f_", lineage_name, "_reference_root_pseudotime")), w = 7.4, h = 5.8)
  save_pdf(p_q, file.path(out_plots, paste0("09f_", lineage_name, "_query_pseudotime_split_genotype")), w = 10.2, h = 5.3)
  save_side_by_side(
    p_left = p_root + theme(legend.position = "right"),
    p_right = p_q + theme(legend.position = "right"),
    path_no_ext = file.path(out_plots, paste0("09f_", lineage_name, "_root_vs_query_split")),
    w = 14.2,
    h = 5.8
  )
}

main <- function() {
  base <- CFG$paths$out_dir
  out_base <- file.path(base, "figures_09f")
  out_plots <- file.path(out_base, "plots")
  out_tbl <- file.path(out_base, "tables")
  dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)

  msg("Loading 09b mapped tables.")
  map <- read_mapping_tables(base)

  msg("Loading roots summary from 09a.")
  roots <- read_roots_summary(base)

  msg("Loading E18 query object.")
  query <- readRDS(CFG$paths$query_e18_rds)
  query_dt <- make_query_combined_df(query_obj = query, map_best = map$best)

  msg("Plotting query combined pseudotime UMAPs.")
  plot_query_combined(query_dt, out_plots)

  root_rows <- vector("list", nrow(roots))
  for (i in seq_len(nrow(roots))) {
    lin <- roots$lineage_name[i]
    root_label <- as.character(roots$root_label[i])
    filter_col <- as.character(roots$filter_label_col[i])

    msg("[09f:", lin, "] plotting reference roots + query split panels.")
    ref_dt <- make_reference_root_df(base, lineage_name = lin, root_label = root_label, filter_label_col = filter_col)
    qry_dt <- make_query_lineage_df(query_dt = query_dt, map_all = map$all, lineage_name = lin)

    if (is.null(qry_dt) || nrow(qry_dt) == 0) {
      msg("[09f:", lin, "] no mapped query cells (skip query panel).")
      next
    }

    plot_lineage_root_and_query(
      ref_dt = ref_dt,
      qry_dt = qry_dt,
      lineage_name = lin,
      root_label = root_label,
      out_plots = out_plots
    )

    root_rows[[i]] <- data.table(
      lineage_name = lin,
      root_label = root_label,
      start_cluster = as.character(roots$start_cluster[i]),
      filter_label_col = filter_col,
      n_ref_cells = nrow(ref_dt),
      n_root_cells = sum(ref_dt$root_flag, na.rm = TRUE),
      n_query_mapped = nrow(qry_dt),
      n_query_wt = sum(qry_dt$genotype2 == "WT", na.rm = TRUE),
      n_query_mut = sum(qry_dt$genotype2 == "MUT", na.rm = TRUE)
    )
  }

  root_tbl <- rbindlist(root_rows, fill = TRUE)
  fwrite(root_tbl, file.path(out_tbl, "09f_roots_used.csv"))

  comb_tbl <- query_dt[, .(
    n_query_cells_total = .N,
    n_query_cells_mapped = sum(!is.na(pseudotime_transfer)),
    n_query_cells_unmapped = sum(is.na(pseudotime_transfer)),
    n_wt = sum(genotype2 == "WT", na.rm = TRUE),
    n_mut = sum(genotype2 == "MUT", na.rm = TRUE)
  )]
  fwrite(comb_tbl, file.path(out_tbl, "09f_query_combined_summary.csv"))

  msg("09f complete. Outputs under: ", out_base)
}

main()
