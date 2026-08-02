# ============================================================
# 04_chipseq_deg_integration.R — integrate CTCF peak loss with E18.5 snRNA-seq DEGs
# ============================================================

source(file.path("scripts", "chipseq", "00_chipseq_setup.R"))
chipseq_require(c(
  "GenomicRanges", "GenomicFeatures", "AnnotationDbi", "org.Mm.eg.db",
  "TxDb.Mmusculus.UCSC.mm10.knownGene", "ggplot2", "dplyr", "readr", "tidyr", "forcats", "purrr"
))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
class_obj_path <- file.path(CHIP_DIRS$objects, "chipseq_02_classified_peaks.rds")
anno_obj_path <- file.path(CHIP_DIRS$objects, "chipseq_03_peak_annotations.rds")
if (!file.exists(class_obj_path) || !file.exists(anno_obj_path)) {
  stop("Run scripts 02 and 03 first: classified peak and annotation objects are required.", call. = FALSE)
}
class_obj <- readRDS(class_obj_path)
anno_obj <- readRDS(anno_obj_path)

key_DOWN_global <- c(
  "Grin1","Grin2a","Grin2b","Gria1","Gria3","Gria4","Shank1","Shank2","Shank3",
  "Dlgap1","Dlgap2","Camk2b","Sparcl1","Bcan","Sst","Sox6","Lhx6","Pvalb",
  "Nrxn1","Nrxn2","Nrxn3","Nlgn1","Nlgn3","Syn2","Syp","Stx1a","Snap25","Scn1a","Kcnq2"
)
key_UP_global <- c(
  "Robo1","Robo2","Efna5","Sema6d","Epha4","Epha7","Dcx","Dlx1","Dlx2","Dlx5",
  "Dlx6","Lin28b","Vim","Sox2","Nes","Ttr","Igfbp2"
)
pseudotime_interaction_genes <- c("Sst","Sox6","Robo1","Camk2b","Grin1","Shank3")
spn_identity_genes <- c("Foxp1","Isl1","Ebf1","Drd2","Penk","Gpr88","Adora2a")
cpcdh_positive_control_genes <- c("Pcdha1","Pcdhb1","Pcdhg")
key_genes <- c(
  "Grin1","Grin2a","Shank1","Shank2","Shank3","Camk2b","Sst","Sox6","Robo1","Dlx1","Dlx2",
  "Sparcl1","Bcan","Nrxn1","Pcdha1","Pcdhb1","Pcdhg","Foxp1","Isl1","Ebf1","Gpr88","Adora2a"
)

background_genes <- unique(readr::read_csv(DE_FILES$global, show_col_types = FALSE)$gene)
background_genes <- sort(unique(background_genes[!is.na(background_genes) & background_genes != ""]))

pseudotime_hits_union <- unique(c(
  read_pseudotime_hits(DE_FILES$pseudo_ventral),
  read_pseudotime_hits(DE_FILES$pseudo_dorsal),
  pseudotime_interaction_genes
))

celltype_gene_sets <- list(
  Global_DOWN = read_deg_directional(DE_FILES$global, "down"),
  Global_UP = read_deg_directional(DE_FILES$global, "up"),
  Synaptic_maturation_curated_DOWN = key_DOWN_global,
  Immature_UP_curated = key_UP_global,
  Key_DOWN_curated = key_DOWN_global,
  Key_UP_curated = key_UP_global,
  Pseudotime_key_curated = pseudotime_interaction_genes,
  Pseudotime_sig_union = pseudotime_hits_union,
  SPN_identity_curated = spn_identity_genes,
  cPcdh_positive_control = cpcdh_positive_control_genes,
  SPN_broad_DOWN = read_deg_directional(DE_FILES$broad_spn, "down"),
  SPN_broad_UP = read_deg_directional(DE_FILES$broad_spn, "up"),
  MGE_broad_DOWN = read_deg_directional(DE_FILES$broad_mge, "down"),
  MGE_broad_UP = read_deg_directional(DE_FILES$broad_mge, "up"),
  Deep_layer_EN_DOWN = read_deg_directional(DE_FILES$broad_deep, "down"),
  Upper_layer_EN_DOWN = read_deg_directional(DE_FILES$broad_upper, "down"),
  Immature_Astrocytes_DOWN = read_deg_directional(DE_FILES$broad_astro, "down"),
  LGE_IN_prec_DOWN = read_deg_directional(DE_FILES$broad_lge_in, "down"),
  Migrating_CGE_IN_DOWN = read_deg_directional(DE_FILES$broad_cge, "down"),
  ExtendedAmygdala_GABA_DOWN = read_deg_directional(DE_FILES$broad_ext, "down"),
  SPN_D1_DOWN = read_deg_directional(DE_FILES$fine_spn_d1, "down"),
  SPN_D2_DOWN = read_deg_directional(DE_FILES$fine_spn_d2, "down"),
  MGE_SSTPV_DOWN = read_deg_directional(DE_FILES$fine_mge, "down")
)
celltype_gene_sets <- lapply(celltype_gene_sets, function(gs) sort(unique(intersect(gs, background_genes))))

manifest_df <- data.frame(
  gene_set = names(celltype_gene_sets),
  n_genes = vapply(celltype_gene_sets, length, integer(1)),
  stringsAsFactors = FALSE
) %>% dplyr::arrange(dplyr::desc(n_genes))
chipseq_write_csv(manifest_df, "chipseq_04_gene_set_manifest.csv")

all_genes_gr <- GenomicFeatures::genes(txdb)
window_sizes <- c(10000L, 50000L, 100000L)
peak_class_list <- list(
  LOST = class_obj$lost_peaks,
  GAINED = class_obj$gained_peaks,
  MAINTAINED = class_obj$maintained_peaks
)

windows_by_kb <- lapply(window_sizes, function(w) build_gene_windows(txdb, upstream = w, downstream = w))
names(windows_by_kb) <- paste0(window_sizes / 1000L, "kb")

get_nearby_symbols <- function(gene_window_gr, peak_gr) {
  hits <- GenomicRanges::countOverlaps(gene_window_gr, peak_gr) > 0
  ids <- names(gene_window_gr)[hits]
  syms <- safe_map_ids(ids, column = "SYMBOL", keytype = "ENTREZID")
  sort(unique(as.character(stats::na.omit(syms))))
}

nearby_symbol_cache <- list()
for (peak_class in names(peak_class_list)) {
  nearby_symbol_cache[[peak_class]] <- list()
  for (win_name in names(windows_by_kb)) {
    nearby_symbol_cache[[peak_class]][[win_name]] <- get_nearby_symbols(windows_by_kb[[win_name]], peak_class_list[[peak_class]])
  }
}

run_fisher_enrichment <- function(target_symbols, background_symbols, nearby_symbols) {
  target_symbols <- sort(unique(intersect(target_symbols, background_symbols)))
  non_target_symbols <- sort(setdiff(background_symbols, target_symbols))
  if (length(target_symbols) == 0 || length(non_target_symbols) == 0) {
    return(data.frame(
      n_target = length(target_symbols),
      n_background = length(background_symbols),
      n_target_with_peak = NA_integer_,
      n_target_without_peak = NA_integer_,
      n_bg_with_peak = NA_integer_,
      n_bg_without_peak = NA_integer_,
      pct_target_with_peak = NA_real_,
      pct_bg_with_peak = NA_real_,
      odds_ratio = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = NA_real_
    ))
  }
  a <- sum(target_symbols %in% nearby_symbols)
  b <- length(target_symbols) - a
  c <- sum(non_target_symbols %in% nearby_symbols)
  d <- length(non_target_symbols) - c
  mat <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE)
  ft <- tryCatch(stats::fisher.test(mat, alternative = "greater"), error = function(e) NULL)
  data.frame(
    n_target = length(target_symbols),
    n_background = length(background_symbols),
    n_target_with_peak = a,
    n_target_without_peak = b,
    n_bg_with_peak = c,
    n_bg_without_peak = d,
    pct_target_with_peak = a / length(target_symbols),
    pct_bg_with_peak = c / length(non_target_symbols),
    odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
    conf_low = if (is.null(ft)) NA_real_ else ft$conf.int[1],
    conf_high = if (is.null(ft)) NA_real_ else ft$conf.int[2],
    p_value = if (is.null(ft)) NA_real_ else ft$p.value
  )
}

enrichment_rows <- list()
idx <- 0L
for (gene_set_name in names(celltype_gene_sets)) {
  genes <- celltype_gene_sets[[gene_set_name]]
  for (peak_class in names(peak_class_list)) {
    for (win_name in names(windows_by_kb)) {
      idx <- idx + 1L
      enrichment_rows[[idx]] <- cbind(
        data.frame(
          gene_set = gene_set_name,
          peak_class = peak_class,
          window_kb = as.integer(gsub("kb", "", win_name)),
          stringsAsFactors = FALSE
        ),
        run_fisher_enrichment(genes, background_genes, nearby_symbol_cache[[peak_class]][[win_name]])
      )
    }
  }
}
enrichment_df <- dplyr::bind_rows(enrichment_rows) %>%
  dplyr::group_by(window_kb) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sig_label = sig_label(p_adj))
chipseq_write_csv(enrichment_df, "chipseq_04_DEG_CTCF_peak_enrichment_Fisher.csv")

p_forest <- enrichment_df %>%
  dplyr::filter(window_kb == 50, gene_set %in% c("Global_DOWN", "Global_UP", "Synaptic_maturation_curated_DOWN", "Immature_UP_curated", "Pseudotime_key_curated", "SPN_identity_curated", "SPN_broad_DOWN", "MGE_broad_DOWN")) %>%
  dplyr::mutate(
    gene_set = forcats::fct_reorder(gene_set, odds_ratio, .na_rm = FALSE),
    peak_class = factor(peak_class, levels = c("LOST", "GAINED", "MAINTAINED"))
  ) %>%
  ggplot(aes(x = odds_ratio, y = gene_set, color = peak_class)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55") +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.18, position = position_dodge(width = 0.55), linewidth = 0.7) +
  geom_point(position = position_dodge(width = 0.55), size = 2.2) +
  scale_color_manual(values = c(
    "LOST" = chipseq_palette$lost,
    "GAINED" = chipseq_palette$gained,
    "MAINTAINED" = chipseq_palette$maintained
  )) +
  theme_classic(base_size = 13) +
  labs(
    title = "DE gene sets near CTCF peak classes",
    subtitle = "50 kb TSS window Fisher enrichment (odds ratio > 1 indicates over-representation near that peak class)",
    x = "Fisher odds ratio",
    y = NULL,
    color = "Peak class"
  )
chipseq_save_plot(p_forest, "chipseq_04_peak_enrichment_forest_50kb", subdir = "de_integration", w = 10.4, h = 6.6)

p_bubble <- enrichment_df %>%
  dplyr::filter(window_kb == 50, gene_set %in% c("Global_DOWN", "Global_UP", "Synaptic_maturation_curated_DOWN", "Immature_UP_curated", "Pseudotime_key_curated", "SPN_identity_curated", "SPN_broad_DOWN", "MGE_broad_DOWN")) %>%
  ggplot(aes(x = peak_class, y = gene_set)) +
  geom_point(aes(size = n_target, color = odds_ratio, alpha = pct_target_with_peak), stroke = 1.1) +
  scale_color_gradient2(low = chipseq_palette$wt, mid = chipseq_palette$light, high = chipseq_palette$coral, midpoint = 1, na.value = "grey80") +
  scale_size_continuous(range = c(2.8, 10)) +
  scale_alpha_continuous(range = c(0.45, 1)) +
  theme_classic(base_size = 13) +
  labs(
    title = "Fraction of DEG sets near LOST / GAINED / MAINTAINED peaks",
    subtitle = "50 kb window; bubble size = gene-set size, color = enrichment odds ratio, alpha = fraction of target genes near peaks",
    x = NULL,
    y = NULL,
    color = "Odds ratio",
    size = "n genes",
    alpha = "% with peak"
  )
chipseq_save_plot(p_bubble, "chipseq_04_peak_enrichment_bubble_50kb", subdir = "de_integration", w = 9.6, h = 6.2)

key_gene_symbols <- sort(unique(c(key_genes, key_DOWN_global, key_UP_global, pseudotime_interaction_genes)))
key_entrez <- safe_map_ids(key_gene_symbols, column = "ENTREZID", keytype = "SYMBOL")
key_map_df <- data.frame(
  gene = names(key_entrez),
  entrez = as.character(key_entrez),
  stringsAsFactors = FALSE
) %>% dplyr::filter(!is.na(entrez), entrez != "")

key_gene_gr <- all_genes_gr[key_map_df$entrez]
if (length(key_gene_gr) > 0) {
  names(key_gene_gr) <- key_map_df$gene[match(names(key_gene_gr), key_map_df$entrez)]
  key_tss_gr <- promoters(key_gene_gr, upstream = 0, downstream = 1)
  nearest_hits <- GenomicRanges::distanceToNearest(key_tss_gr, class_obj$lost_peaks, ignore.strand = TRUE)
  nearest_df <- data.frame(
    gene = names(key_tss_gr)[S4Vectors::queryHits(nearest_hits)],
    nearest_idx = S4Vectors::subjectHits(nearest_hits),
    distance_to_TSS_bp = S4Vectors::mcols(nearest_hits)$distance,
    stringsAsFactors = FALSE
  )

  key_info_df <- data.frame(
    gene = names(key_gene_gr),
    chr = as.character(seqnames(key_gene_gr)),
    gene_start = start(key_gene_gr),
    gene_end = end(key_gene_gr),
    strand = as.character(strand(key_gene_gr)),
    gene_width = width(key_gene_gr),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(tss = ifelse(strand == "-", gene_end, gene_start))

  flags_df <- cbind(
    gene = names(key_tss_gr),
    window_overlap_flags(key_tss_gr, class_obj$lost_peaks, distances = c(10000L, 50000L, 100000L, 500000L))
  )
  gene_body_overlap <- data.frame(
    gene = names(key_gene_gr),
    lost_peak_in_gene_body = GenomicRanges::countOverlaps(key_gene_gr, class_obj$lost_peaks, ignore.strand = TRUE) > 0,
    stringsAsFactors = FALSE
  )

  lost_df <- granges_to_df(class_obj$lost_peaks)
  nearest_lost_df <- nearest_df %>%
    dplyr::mutate(
      nearest_LOST_peak_chr = as.character(seqnames(class_obj$lost_peaks))[nearest_idx],
      nearest_LOST_peak_start = start(class_obj$lost_peaks)[nearest_idx],
      nearest_LOST_peak_end = end(class_obj$lost_peaks)[nearest_idx],
      nearest_LOST_peak_signalValue = class_obj$lost_peaks$signalValue[nearest_idx],
      nearest_LOST_peak_qValue = class_obj$lost_peaks$qValue[nearest_idx]
    ) %>%
    dplyr::select(-nearest_idx)

  key_lost_summary <- key_info_df %>%
    dplyr::left_join(nearest_lost_df, by = "gene") %>%
    dplyr::left_join(gene_body_overlap, by = "gene") %>%
    dplyr::left_join(flags_df, by = "gene")
} else {
  key_lost_summary <- data.frame()
}
chipseq_write_csv(key_lost_summary, "chipseq_04_key_gene_nearest_LOST_peak.csv")

if (nrow(key_lost_summary) > 0) {
  p_key <- key_lost_summary %>%
    dplyr::mutate(
      gene = forcats::fct_reorder(gene, distance_to_TSS_bp, .desc = TRUE),
      label_50kb = ifelse(within_50kb, "LOST peak within 50 kb", "No LOST peak within 50 kb")
    ) %>%
    ggplot(aes(x = distance_to_TSS_bp, y = gene, color = label_50kb)) +
    geom_segment(aes(x = 0, xend = distance_to_TSS_bp, y = gene, yend = gene), linewidth = 1) +
    geom_point(size = 2.4) +
    scale_color_manual(values = c("LOST peak within 50 kb" = chipseq_palette$coral, "No LOST peak within 50 kb" = chipseq_palette$navy)) +
    theme_classic(base_size = 13) +
    labs(
      title = "Nearest LOST CTCF peak to key snRNA-seq genes",
      subtitle = "Distance is measured from the gene TSS to the nearest WT-specific (LOST) CTCF peak",
      x = "Distance to nearest LOST peak (bp)",
      y = NULL,
      color = NULL
    )
  chipseq_save_plot(p_key, "chipseq_04_key_gene_nearest_LOST_peak", subdir = "de_integration", w = 9.4, h = 7.2)
}

chipseq_msg("Script 04 complete: DEG integration tables and enrichment plots written.")
