# ============================================================
# 05_withinstate_DE_GO.R — within-state DE (mut vs wt) + GO
# ============================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

use_kept_only <- TRUE
use_cortical_striatal_scope <- TRUE
analysis_scope_label <- "cortical/striatal analysis scope"

in_obj_full <- file.path(OUT$objects, "03_seu_annotated.rds")
in_obj_kept <- file.path(OUT$objects, "03_seu_annotated_KEPT.rds")
in_obj_cortstr <- file.path(OUT$objects, "03_seu_annotated_CORTSTR.rds")
in_obj_tel  <- file.path(OUT$objects, "03_seu_annotated_TEL.rds")
in_obj <- if (use_cortical_striatal_scope && file.exists(in_obj_cortstr)) {
  in_obj_cortstr
} else if (use_cortical_striatal_scope && file.exists(in_obj_tel)) {
  in_obj_tel
} else if (use_kept_only) {
  in_obj_kept
} else {
  in_obj_full
}
stopifnot(file.exists(in_obj))

seu <- readRDS(in_obj)
message("Loaded DE object: ", in_obj)
stopifnot("celltype_fine" %in% colnames(seu@meta.data))
stopifnot("celltype_broad" %in% colnames(seu@meta.data))

# Keep only analysis-eligible cells.
if ("include_in_analysis" %in% colnames(seu@meta.data)) {
  seu <- subset(seu, subset = include_in_analysis)
}

DefaultAssay(seu) <- "SCT"
# Re-correct SCT models once globally to avoid per-subset model-size mismatch errors.
seu <- PrepSCTFindMarkers(seu)

# ---- Volcano helper ----
volcano_plot <- function(df, title){
  stopifnot(all(c("gene","p_val_adj","avg_log2FC") %in% colnames(df)))
  df %>%
    mutate(
      de_class = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0 ~ "Up in mut",
        p_val_adj < 0.05 & avg_log2FC < 0 ~ "Down in mut",
        TRUE ~ "Non-significant"
      )
    ) %>%
    ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = de_class), size=1.2, alpha=0.8) +
    theme_classic(13) +
    labs(title=title, x="avg_log2FC (mut vs wt)", y="-log10(FDR)") +
    scale_color_manual(
      values = c(
        "Non-significant" = "gray70",
        "Up in mut" = "#d62728",
        "Down in mut" = "#1f77b4"
      )
    )
}

# ---- GO helper (BP) ----
run_go_bp <- function(genes, universe, prefix){
  eg <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  uni <- bitr(universe, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  
  if (nrow(eg) < 10 || nrow(uni) < 100) return(NULL)
  
  enr <- enrichGO(
    gene = eg$ENTREZID,
    universe = uni$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  if (is.null(enr) || nrow(enr@result)==0) return(NULL)
  
  write_csv(as.data.frame(enr), file.path(OUT$tables, paste0(prefix, "_GO_BP.csv")))
  
  p <- dotplot(enr, showCategory = 15) + ggtitle(paste0(prefix, " | ", analysis_scope_label))
  folder <- if (startsWith(prefix, "celltype_broad_")) {
    "de_broad"
  } else if (startsWith(prefix, "global_")) {
    "de_global"
  } else {
    "de_fine"
  }
  save_plot(p, paste0(folder, "/GO_", prefix), w=9, h=7)
  enr
}

# ---- Within-state DE loops (fine + broad) ----
min_cells_per_group_fine <- 50
min_cells_per_group_broad <- 30
min_genes_for_go <- 10
min_universe_for_go <- 100

run_de_by_level <- function(level_col, level_tag, min_cells_per_group) {
  seu[[level_col]] <- droplevels(as.factor(seu[[level_col, drop = TRUE]]))
  Idents(seu) <- level_col
  celltypes <- levels(seu[[level_col, drop = TRUE]])
  # Exclude placeholder/non-analysis labels if present.
  celltypes <- celltypes[!grepl("^Excluded", celltypes)]
  summary_list <- list()
  run_log <- list()
  de_results <- list()
  
  for (ct in celltypes) {
    message("DE [", level_col, "]: ", ct)
    tryCatch({
      if (!(ct %in% levels(Idents(seu)))) {
        msg <- "skip (identity not found in object levels)"
        message("  ", msg)
        run_log[[ct]] <- data.frame(annotation_level = level_col, celltype = ct, status = "skipped", reason = msg)
        next
      }
      
      obj <- subset(seu, idents = ct)
      tab <- table(obj$condition)
      
      if (!all(c("wt","mut") %in% names(tab)) || any(tab[c("wt","mut")] < min_cells_per_group)) {
        msg <- paste("skip (cells wt/mut too low):", paste(names(tab), tab, collapse=" / "))
        message("  ", msg)
        run_log[[ct]] <- data.frame(annotation_level = level_col, celltype = ct, status = "skipped", reason = msg)
      } else {
        de <- FindMarkers(
          obj,
          assay = "SCT",
          group.by = "condition",
          ident.1 = "mut",
          ident.2 = "wt",
          test.use = "wilcox",
          logfc.threshold = 0,
          min.pct = 0.1,
          recorrect_umi = FALSE
        )
        
        de <- de %>%
          as.data.frame() %>%
          tibble::rownames_to_column("gene") %>%
          arrange(p_val_adj)
        
        safe_ct <- gsub("[^A-Za-z0-9]+","_", ct)
        write_csv(de, file.path(OUT$tables, paste0("DE_", level_tag, "_", safe_ct, "_mut_vs_wt.csv")))
        de_results[[ct]] <- de
        
        pvol <- volcano_plot(de, paste0(level_col, ": ", ct, " (mut vs wt) | ", analysis_scope_label))
        de_folder <- if (level_tag == "celltype_broad") "de_broad" else "de_fine"
        save_plot(pvol, paste0(de_folder, "/Volcano_", level_tag, "_", safe_ct), w=8, h=6)
        
        sig <- de %>% filter(p_val_adj < 0.05)
        up <- sig %>% filter(avg_log2FC > 0) %>% pull(gene)
        down <- sig %>% filter(avg_log2FC < 0) %>% pull(gene)
        universe <- de$gene
        
        if (length(up) >= min_genes_for_go && length(universe) >= min_universe_for_go) {
          run_go_bp(up, universe, paste0(level_tag, "_", safe_ct, "_UP_in_mut"))
        }
        if (length(down) >= min_genes_for_go && length(universe) >= min_universe_for_go) {
          run_go_bp(down, universe, paste0(level_tag, "_", safe_ct, "_DOWN_in_mut"))
        }
        
        summary_list[[ct]] <- data.frame(
          annotation_level = level_col,
          celltype = ct,
          n_wt = as.integer(tab["wt"]),
          n_mut = as.integer(tab["mut"]),
          n_sig = nrow(sig),
          n_up = length(up),
          n_down = length(down)
        )
        run_log[[ct]] <- data.frame(annotation_level = level_col, celltype = ct, status = "ok", reason = NA_character_)
      }
    }, error = function(e) {
      msg <- paste("error:", conditionMessage(e))
      message("  ", msg)
      run_log[[ct]] <- data.frame(annotation_level = level_col, celltype = ct, status = "error", reason = msg)
    })
  }
  
  summary_df <- bind_rows(summary_list)
  write_csv(summary_df, file.path(OUT$tables, paste0("DE_summary_by_", level_tag, ".csv")))
  write_csv(bind_rows(run_log), file.path(OUT$tables, paste0("DE_runlog_", level_tag, ".csv")))
  invisible(list(summary = summary_df, runlog = bind_rows(run_log), de_results = de_results))
}

res_broad <- run_de_by_level("celltype_broad", "celltype_broad", min_cells_per_group_broad)
res_fine <- run_de_by_level("celltype_fine", "celltype_fine", min_cells_per_group_fine)

# ---- Global DE across all included cells (mut vs wt) ----
Idents(seu) <- "condition"
tab_all <- table(seu$condition)
if (all(c("wt", "mut") %in% names(tab_all)) && all(tab_all[c("wt", "mut")] > 0)) {
  message("DE [global all cells]: mut vs wt")
  
  # Unadjusted all-cell DE (can capture strong composition+state shifts).
  de_global_wilcox <- FindMarkers(
    seu,
    assay = "SCT",
    ident.1 = "mut",
    ident.2 = "wt",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0.1,
    recorrect_umi = FALSE
  ) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj)
  write_csv(de_global_wilcox, file.path(OUT$tables, "DE_global_allcells_wilcox_mut_vs_wt.csv"))
  save_plot(
    volcano_plot(de_global_wilcox, paste0("All cells (global): mut vs wt [Wilcoxon] | ", analysis_scope_label)),
    "de_global/Volcano_global_allcells_wilcox_mut_vs_wt", w = 8, h = 6
  )
  
  sig_w <- de_global_wilcox %>% filter(p_val_adj < 0.05)
  up_w <- sig_w %>% filter(avg_log2FC > 0) %>% pull(gene)
  down_w <- sig_w %>% filter(avg_log2FC < 0) %>% pull(gene)
  if (length(up_w) >= min_genes_for_go && nrow(de_global_wilcox) >= min_universe_for_go) {
    run_go_bp(up_w, de_global_wilcox$gene, "global_allcells_wilcox_UP_in_mut")
  }
  if (length(down_w) >= min_genes_for_go && nrow(de_global_wilcox) >= min_universe_for_go) {
    run_go_bp(down_w, de_global_wilcox$gene, "global_allcells_wilcox_DOWN_in_mut")
  }
  
  # Celltype-adjusted global DE (tests condition effect while accounting for celltype label).
  latent_vars <- c()
  if ("celltype_broad" %in% colnames(seu@meta.data)) latent_vars <- c(latent_vars, "celltype_broad")
  if ("nCount_RNA" %in% colnames(seu@meta.data)) latent_vars <- c(latent_vars, "nCount_RNA")
  if ("percent.mt" %in% colnames(seu@meta.data)) latent_vars <- c(latent_vars, "percent.mt")
  
  if (length(latent_vars) > 0) {
    de_global_lr <- FindMarkers(
      seu,
      assay = "SCT",
      ident.1 = "mut",
      ident.2 = "wt",
      test.use = "LR",
      latent.vars = latent_vars,
      logfc.threshold = 0,
      min.pct = 0.1,
      recorrect_umi = FALSE
    ) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("gene") %>%
      arrange(p_val_adj)
    
    write_csv(de_global_lr, file.path(OUT$tables, "DE_global_allcells_LR_celltypeadjusted_mut_vs_wt.csv"))
    save_plot(
      volcano_plot(de_global_lr, paste0("All cells (global): mut vs wt [LR + celltype-adjusted] | ", analysis_scope_label)),
      "de_global/Volcano_global_allcells_LR_celltypeadjusted_mut_vs_wt", w = 8, h = 6
    )
    
    sig_lr <- de_global_lr %>% filter(p_val_adj < 0.05)
    up_lr <- sig_lr %>% filter(avg_log2FC > 0) %>% pull(gene)
    down_lr <- sig_lr %>% filter(avg_log2FC < 0) %>% pull(gene)
    if (length(up_lr) >= min_genes_for_go && nrow(de_global_lr) >= min_universe_for_go) {
      run_go_bp(up_lr, de_global_lr$gene, "global_allcells_LR_celltypeadjusted_UP_in_mut")
    }
    if (length(down_lr) >= min_genes_for_go && nrow(de_global_lr) >= min_universe_for_go) {
      run_go_bp(down_lr, de_global_lr$gene, "global_allcells_LR_celltypeadjusted_DOWN_in_mut")
    }
  }
}

# ---- Cross-celltype consistency summary (fine labels) ----
if (length(res_fine$de_results) > 0) {
  de_long <- bind_rows(lapply(names(res_fine$de_results), function(ct) {
    res_fine$de_results[[ct]] %>%
      transmute(
        celltype = ct,
        gene = gene,
        p_val_adj = p_val_adj,
        avg_log2FC = avg_log2FC,
        sig = !is.na(p_val_adj) & p_val_adj < 0.05,
        sig_up = !is.na(p_val_adj) & p_val_adj < 0.05 & avg_log2FC > 0,
        sig_down = !is.na(p_val_adj) & p_val_adj < 0.05 & avg_log2FC < 0
      )
  }))
  write_csv(de_long, file.path(OUT$tables, "DE_fine_longform_gene_by_celltype.csv"))
  
  de_consensus <- de_long %>%
    group_by(gene) %>%
    summarise(
      n_celltypes_tested = n(),
      n_sig = sum(sig, na.rm = TRUE),
      n_sig_up = sum(sig_up, na.rm = TRUE),
      n_sig_down = sum(sig_down, na.rm = TRUE),
      frac_sig = n_sig / n_celltypes_tested,
      frac_sig_up = n_sig_up / n_celltypes_tested,
      frac_sig_down = n_sig_down / n_celltypes_tested,
      median_log2FC = median(avg_log2FC, na.rm = TRUE),
      mean_log2FC = mean(avg_log2FC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      direction_consistency = case_when(
        n_sig_up > 0 & n_sig_down == 0 ~ "consistent_up",
        n_sig_down > 0 & n_sig_up == 0 ~ "consistent_down",
        n_sig_up > 0 & n_sig_down > 0 ~ "mixed",
        TRUE ~ "not_significant"
      )
    ) %>%
    arrange(desc(n_sig), desc(abs(median_log2FC)))
  
  write_csv(de_consensus, file.path(OUT$tables, "DE_fine_gene_consistency_summary.csv"))
  
  min_consistent_celltypes <- max(2, ceiling(0.30 * length(res_fine$de_results)))
  universe_cons <- unique(de_long$gene)
  up_cons <- de_consensus %>% filter(n_sig_up >= min_consistent_celltypes, n_sig_down == 0) %>% pull(gene)
  down_cons <- de_consensus %>% filter(n_sig_down >= min_consistent_celltypes, n_sig_up == 0) %>% pull(gene)
  
  if (length(up_cons) >= min_genes_for_go && length(universe_cons) >= min_universe_for_go) {
    run_go_bp(up_cons, universe_cons, paste0("global_consistent_fine_UP_nge", min_consistent_celltypes))
  }
  if (length(down_cons) >= min_genes_for_go && length(universe_cons) >= min_universe_for_go) {
    run_go_bp(down_cons, universe_cons, paste0("global_consistent_fine_DOWN_nge", min_consistent_celltypes))
  }
}
