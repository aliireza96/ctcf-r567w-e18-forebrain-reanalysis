#!/usr/bin/env Rscript
# ============================================================
# 09_chipseq_peak_signal_gene_fc_correlation.R
#
# Continuous peak signal log2FC (mutant / WT) versus nearest-gene
# log2FC (snRNA-seq DE). Stratified by: (i) genomic feature class
# (promoter vs distal), (ii) peak-signal quartile, (iii) TAD
# boundary proximity.
#
# Inputs:
#   - chipseq_02_classified_peaks.rds   peak class membership
#   - chipseq_03_peak_annotations.rds   genomic annotation per peak
#   - GSM6614265_Brain_CTCF_wt.mm10.bw
#   - GSM6614266_Brain_CTCF_homo.mm10.bw
#   - tad_boundary_insulation_summary.csv (from script 16)
#   - global DE CSV (from 05_withinstate_DE_GO.R)
#
# Outputs:
#   tables/chipseq_09_peak_signal_fc_per_peak.csv
#   tables/chipseq_09_peak_gene_fc_correlation.csv
#   plots/chipseq_09_peak_gene_fc_scatter_by_feature.pdf
#   plots/chipseq_09_peak_gene_fc_scatter_by_quartile.pdf
#   plots/chipseq_09_peak_gene_fc_scatter_tad_boundary.pdf
#   plots/chipseq_09_peak_gene_rho_barplot.pdf
# ============================================================

source(file.path("scripts", "chipseq", "00_chipseq_setup.R"))

# Load all Bioconductor class packages before any RDS deserialisation.
# Order matters: GenomeInfoDb must precede GenomicRanges so Seqinfo is
# registered before any GRanges class dispatch occurs.
chipseq_require(c(
  "GenomeInfoDb", "BiocGenerics", "S4Vectors", "IRanges",
  "GenomicRanges", "GenomicFeatures", "rtracklayer",
  "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db",
  "AnnotationDbi", "ggplot2", "dplyr", "readr", "tidyr", "scales"
))

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
msg <- chipseq_msg

# ── Safe GRanges -> data.frame extractor ─────────────────────────────────────
# Converts a GRanges to a plain data.frame immediately after readRDS so no
# S4 class operations (c, reduce, findOverlaps) are ever called on the raw
# serialised objects. This avoids the Seqinfo class-registry error.
gr_as_df <- function(gr, class_name = NA_character_) {
  if (is.null(gr) || length(gr) == 0) return(NULL)
  data.frame(
    chrom      = as.character(GenomicRanges::seqnames(gr)),
    start      = as.integer(BiocGenerics::start(gr)),
    end        = as.integer(BiocGenerics::end(gr)),
    peak_class = class_name,
    stringsAsFactors = FALSE
  )
}

# ── 1. Load and immediately flatten classified peaks ─────────────────────────
msg("Loading peak objects...")
class_obj <- readRDS(file.path(CHIP_DIRS$objects, "chipseq_02_classified_peaks.rds"))
anno_obj  <- readRDS(file.path(CHIP_DIRS$objects, "chipseq_03_peak_annotations.rds"))

# Extract each class to a data frame right away -- no further S4 ops on these
lost_df  <- gr_as_df(class_obj$lost_peaks,       "LOST")
mnt_df   <- gr_as_df(class_obj$maintained_peaks, "MAINTAINED")
gain_df  <- gr_as_df(class_obj$gained_peaks,     "GAINED")

# Remove class_obj from environment to avoid accidental S4 dispatch later
rm(class_obj)

# Union data frame: LOST + MAINTAINED + GAINED
all_df <- rbind(lost_df, mnt_df, gain_df)
# Deduplicate by coordinates (keep first assignment = priority LOST > MAINTAINED > GAINED)
all_df <- all_df[!duplicated(all_df[c("chrom", "start", "end")]), ]
all_df <- all_df[order(all_df$chrom, all_df$start, all_df$end), ]
all_df$peak_id <- paste0("peak_", seq_len(nrow(all_df)))

msg("  Total peaks in union: ", nrow(all_df),
    " (LOST=", sum(all_df$peak_class == "LOST"),
    ", MAINTAINED=", sum(all_df$peak_class == "MAINTAINED"),
    ", GAINED=", sum(all_df$peak_class == "GAINED"), ")")

# Build a fresh GRanges from the data frame (no Seqinfo inherited from RDS)
all_peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(
  all_df,
  seqnames.field    = "chrom",
  start.field       = "start",
  end.field         = "end",
  keep.extra.columns = FALSE
)
names(all_peaks_gr) <- all_df$peak_id

# ── 2. BigWig signal per peak ─────────────────────────────────────────────────
BW_WT   <- CHIP_FILES$wt_bw
BW_HOMO <- CHIP_FILES$mut_bw
if (!file.exists(BW_WT))   stop("WT bigWig not found: ",   BW_WT)
if (!file.exists(BW_HOMO)) stop("HOMO bigWig not found: ", BW_HOMO)

msg("Scoring WT bigWig signal over ", length(all_peaks_gr), " peaks...")
sig_wt   <- bw_mean_on_ranges(BW_WT,   all_peaks_gr)
msg("Scoring HOMO bigWig signal...")
sig_homo <- bw_mean_on_ranges(BW_HOMO, all_peaks_gr)

# log2FC with pseudocount 0.01
peak_fc <- all_df
peak_fc$sig_wt      <- sig_wt
peak_fc$sig_homo    <- sig_homo
peak_fc$log2fc_peak <- log2((sig_homo + 0.01) / (sig_wt + 0.01))

# WT signal quartile
wt_q <- quantile(peak_fc$sig_wt, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
peak_fc$wt_signal_quartile <- as.character(cut(
  peak_fc$sig_wt,
  breaks = wt_q,
  labels = c("Q1_weak", "Q2", "Q3", "Q4_strong"),
  include.lowest = TRUE
))

msg("  Peak signal computed for ", nrow(peak_fc), " peaks")

# ── 3. Nearest-gene assignment (50 kb window) ─────────────────────────────────
msg("Assigning peaks to nearest gene within 50 kb...")
txdb     <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
genes_gr <- suppressMessages(GenomicFeatures::genes(txdb))

gene_sym <- tryCatch(
  AnnotationDbi::mapIds(
    org.Mm.eg.db::org.Mm.eg.db,
    keys    = as.character(genes_gr$gene_id),
    column  = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  ),
  error = function(e) { msg("  mapIds failed: ", conditionMessage(e)); NULL }
)
if (!is.null(gene_sym)) {
  genes_gr$symbol <- unname(gene_sym[as.character(genes_gr$gene_id)])
} else {
  genes_gr$symbol <- as.character(genes_gr$gene_id)
}

# distanceToNearest works on fresh all_peaks_gr (no Seqinfo from RDS)
hits      <- GenomicRanges::distanceToNearest(all_peaks_gr, genes_gr, ignore.strand = TRUE)
close_idx <- which(S4Vectors::mcols(hits)$distance <= 50000)

peak_fc$nearest_gene <- NA_character_
peak_fc$dist_to_gene <- NA_real_
if (length(close_idx) > 0) {
  qh <- S4Vectors::queryHits(hits)[close_idx]
  sh <- S4Vectors::subjectHits(hits)[close_idx]
  peak_fc$nearest_gene[qh] <- genes_gr$symbol[sh]
  peak_fc$dist_to_gene[qh] <- S4Vectors::mcols(hits)$distance[close_idx]
}
msg("  Peaks with gene within 50 kb: ",
    sum(!is.na(peak_fc$nearest_gene)), " / ", nrow(peak_fc))

# ── 4. Genomic feature class from annotation object ──────────────────────────
peak_fc$genomic_feature <- NA_character_

if (!is.null(anno_obj) && !is.null(anno_obj$anno_df)) {
  anno_df  <- anno_obj$anno_df
  ann_col  <- intersect(c("annotation", "Annotation", "feature", "Feature"),
                        colnames(anno_df))
  if (length(ann_col) > 0) {
    ann_col <- ann_col[1]
    feature_map <- function(ann) {
      ann <- tolower(as.character(ann))
      ifelse(grepl("promoter", ann), "Promoter",
             ifelse(grepl("5'|3'|exon|intron", ann), "Genic", "Distal"))
    }
    # Try to match by coordinates rather than assuming row order
    if (all(c("chrom","start","end") %in% colnames(anno_df))) {
      anno_key      <- paste0(anno_df$chrom, ":", anno_df$start, "-", anno_df$end)
      peak_key      <- paste0(peak_fc$chrom, ":", peak_fc$start, "-", peak_fc$end)
      feat_vec      <- feature_map(anno_df[[ann_col]])
      names(feat_vec) <- anno_key
      peak_fc$genomic_feature <- unname(feat_vec[peak_key])
    }
  }
}
# Fallback: call promoter if dist_to_gene < 2000
na_feat <- is.na(peak_fc$genomic_feature)
if (any(na_feat)) {
  peak_fc$genomic_feature[na_feat] <- ifelse(
    !is.na(peak_fc$dist_to_gene[na_feat]) & peak_fc$dist_to_gene[na_feat] < 2000,
    "Promoter", "Distal"
  )
}

# ── 5. TAD boundary proximity ─────────────────────────────────────────────────
TAD_CSV <- file.path(E18P5_ROOT, "results", "discovery_hic_14", "tables",
                     "tad_boundary_insulation_summary.csv")
peak_fc$near_tad_boundary      <- NA
peak_fc$near_lost_tad_boundary <- NA
peak_fc$dist_to_tad            <- NA_real_

if (file.exists(TAD_CSV)) {
  msg("  Loading TAD boundary positions from Hi-C insulation analysis...")
  tad_df <- readr::read_csv(TAD_CSV, show_col_types = FALSE)
  # Build fresh GRanges from CSV (no Seqinfo baggage)
  tad_gr <- GenomicRanges::makeGRangesFromDataFrame(
    tad_df,
    seqnames.field = "chrom",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE
  )
  hits_tad <- GenomicRanges::distanceToNearest(all_peaks_gr, tad_gr, ignore.strand = TRUE)
  dist_vec <- S4Vectors::mcols(hits_tad)$distance
  peak_fc$near_tad_boundary <- FALSE
  peak_fc$dist_to_tad[S4Vectors::queryHits(hits_tad)] <- dist_vec
  near_tad_qh <- S4Vectors::queryHits(hits_tad)[dist_vec <= 50000]
  peak_fc$near_tad_boundary[near_tad_qh] <- TRUE

  if ("boundary_lost_in_homo" %in% colnames(tad_df)) {
    lost_mask <- tad_df$boundary_lost_in_homo %in% c(TRUE, "TRUE", "true")
    lost_tad_gr <- tad_gr[lost_mask]
    if (length(lost_tad_gr) > 0) {
      hits_lost <- GenomicRanges::distanceToNearest(all_peaks_gr, lost_tad_gr,
                                                     ignore.strand = TRUE)
      dist_lost <- S4Vectors::mcols(hits_lost)$distance
      peak_fc$near_lost_tad_boundary <- FALSE
      nl_qh <- S4Vectors::queryHits(hits_lost)[dist_lost <= 50000]
      peak_fc$near_lost_tad_boundary[nl_qh] <- TRUE
    }
    msg("  TAD boundaries: ", nrow(tad_df), " total | ",
        sum(lost_mask, na.rm = TRUE), " lost in HOMO")
    msg("  Peaks near any TAD boundary (<=50kb): ",
        sum(peak_fc$near_tad_boundary, na.rm = TRUE))
    msg("  Peaks near LOST TAD boundary (<=50kb): ",
        sum(peak_fc$near_lost_tad_boundary, na.rm = TRUE))
  }
} else {
  msg("  TAD boundary CSV not found; run 16_hic_insulation_tad_boundaries.py first.")
  msg("  Proceeding without TAD stratification.")
}

# ── 6. Load global DE log2FC ──────────────────────────────────────────────────
msg("Loading snRNA-seq DE results...")
de_global <- tryCatch(
  readr::read_csv(DE_FILES$global, show_col_types = FALSE),
  error = function(e) {
    msg("  Could not load global DE file: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(de_global)) {
  gene_col  <- intersect(c("gene", "gene_name", "Gene"),  colnames(de_global))[1]
  fc_col    <- intersect(c("avg_log2FC", "log2FoldChange", "logFC"), colnames(de_global))[1]
  padj_col  <- intersect(c("p_val_adj", "padj", "FDR"),   colnames(de_global))[1]
  if (!is.na(gene_col) && !is.na(fc_col) && !is.na(padj_col)) {
    de_global <- de_global %>%
      dplyr::rename(gene        = !!gene_col,
                    log2fc_gene = !!fc_col,
                    padj_gene   = !!padj_col) %>%
      dplyr::group_by(gene) %>%
      dplyr::slice_min(padj_gene, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, log2fc_gene, padj_gene)
    msg("  Loaded DE results for ", nrow(de_global), " genes")
  } else {
    msg("  DE file column names unrecognised; skipping gene FC merge")
    de_global <- NULL
  }
}

# ── 7. Merge peak FC with gene FC ────────────────────────────────────────────
df <- peak_fc %>% dplyr::filter(!is.na(nearest_gene))
if (!is.null(de_global)) {
  df <- dplyr::left_join(df, de_global, by = c("nearest_gene" = "gene"))
}

# Save per-peak table
chipseq_write_csv(df, "chipseq_09_peak_signal_fc_per_peak.csv")
msg("  Saved per-peak FC table: ", nrow(df), " peaks with assigned gene")

# ── 8. Correlation analysis ───────────────────────────────────────────────────
msg("Running Spearman correlation analyses...")

run_cor <- function(sub_df, label) {
  sub_df <- sub_df %>% dplyr::filter(!is.na(log2fc_peak) & !is.na(log2fc_gene))
  if (nrow(sub_df) < 20) return(NULL)
  ct <- cor.test(sub_df$log2fc_peak, sub_df$log2fc_gene,
                  method = "spearman", exact = FALSE)
  data.frame(
    group   = label,
    n_peaks = nrow(sub_df),
    rho     = round(ct$estimate, 4),
    p_value = signif(ct$p.value, 4),
    stringsAsFactors = FALSE
  )
}

cor_results <- list()

# Overall
cor_results[["all"]] <- run_cor(df, "All peaks")

# By peak class
for (cls in c("LOST", "MAINTAINED", "GAINED")) {
  cor_results[[cls]] <- run_cor(
    df %>% dplyr::filter(peak_class == cls),
    paste0("Peak class: ", cls)
  )
}

# By genomic feature x peak class
for (feat in c("Promoter", "Distal")) {
  for (cls in c("LOST", "MAINTAINED", "GAINED")) {
    lbl <- paste0(feat, " | ", cls)
    cor_results[[lbl]] <- run_cor(
      df %>% dplyr::filter(genomic_feature == feat & peak_class == cls), lbl
    )
  }
}

# By WT signal quartile
for (q in c("Q1_weak", "Q2", "Q3", "Q4_strong")) {
  cor_results[[q]] <- run_cor(
    df %>% dplyr::filter(wt_signal_quartile == q),
    paste0("WT signal: ", q)
  )
}

# TAD boundary stratification
if (!all(is.na(df$near_tad_boundary))) {
  cor_results[["near_TAD"]] <- run_cor(
    df %>% dplyr::filter(near_tad_boundary == TRUE),
    "Near TAD boundary (<=50kb)"
  )
  cor_results[["distal_TAD"]] <- run_cor(
    df %>% dplyr::filter(near_tad_boundary == FALSE),
    "Not near TAD boundary"
  )
  if ("near_lost_tad_boundary" %in% colnames(df) &&
      !all(is.na(df$near_lost_tad_boundary))) {
    cor_results[["near_lost_TAD"]] <- run_cor(
      df %>% dplyr::filter(near_lost_tad_boundary == TRUE),
      "Near LOST TAD boundary (<=50kb)"
    )
  }
}

cor_tbl <- dplyr::bind_rows(Filter(Negate(is.null), cor_results))
chipseq_write_csv(cor_tbl, "chipseq_09_peak_gene_fc_correlation.csv")
msg("  Correlation table (", nrow(cor_tbl), " comparisons):")
print(cor_tbl, row.names = FALSE)

# ── 9. Plots ──────────────────────────────────────────────────────────────────
msg("Plotting...")

class_colors <- c(LOST = "#C44E52", MAINTAINED = "#4878CF", GAINED = "#6ACC65")

base_scatter <- function(sub_df, col_var = "peak_class", title = NULL) {
  sub_df <- sub_df %>% dplyr::filter(!is.na(log2fc_peak) & !is.na(log2fc_gene))
  if (nrow(sub_df) < 10) return(NULL)
  ggplot2::ggplot(sub_df,
    ggplot2::aes(.data[["log2fc_peak"]], .data[["log2fc_gene"]],
                 color = .data[[col_var]])) +
    ggplot2::geom_point(size = 0.4, alpha = 0.35) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.7,
                         ggplot2::aes(group = 1), color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    ggplot2::scale_color_manual(values = class_colors, na.value = "grey60") +
    ggplot2::labs(x = "Peak log2FC (HOMO/WT)",
                  y = "Gene log2FC (snRNA-seq HOMO/WT)",
                  title = title) +
    ggplot2::theme_classic(base_size = 11)
}

# Feature-stratified
p_feat <- df %>%
  dplyr::filter(!is.na(genomic_feature) & !is.na(log2fc_gene)) %>%
  base_scatter(title = "Peak log2FC vs gene log2FC")
if (!is.null(p_feat)) {
  p_feat <- p_feat + ggplot2::facet_wrap(~ genomic_feature, ncol = 2)
  ggplot2::ggsave(
    file.path(CHIP_DIRS$plots, "chipseq_09_peak_gene_fc_scatter_by_feature.pdf"),
    p_feat, width = 10, height = 5
  )
}

# Quartile-stratified
p_q <- df %>%
  dplyr::filter(!is.na(wt_signal_quartile) & !is.na(log2fc_gene)) %>%
  base_scatter(title = "Peak log2FC vs gene log2FC by WT signal quartile")
if (!is.null(p_q)) {
  p_q <- p_q + ggplot2::facet_wrap(~ wt_signal_quartile, ncol = 4)
  ggplot2::ggsave(
    file.path(CHIP_DIRS$plots, "chipseq_09_peak_gene_fc_scatter_by_quartile.pdf"),
    p_q, width = 14, height = 4
  )
}

# TAD boundary scatter
if (!all(is.na(df$near_tad_boundary))) {
  p_tad <- df %>%
    dplyr::mutate(
      tad_label = dplyr::case_when(
        near_tad_boundary == TRUE  ~ "Near TAD boundary",
        near_tad_boundary == FALSE ~ "Interior",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(log2fc_gene) & !is.na(tad_label)) %>%
    base_scatter(title = "Peak log2FC vs gene log2FC by TAD boundary proximity")
  if (!is.null(p_tad)) {
    p_tad <- p_tad + ggplot2::facet_wrap(~ tad_label, ncol = 2)
    ggplot2::ggsave(
      file.path(CHIP_DIRS$plots, "chipseq_09_peak_gene_fc_scatter_tad_boundary.pdf"),
      p_tad, width = 10, height = 5
    )
  }
}

# Summary bar plot: rho by group
if (nrow(cor_tbl) > 0) {
  p_rho <- ggplot2::ggplot(
    cor_tbl %>% dplyr::filter(!is.na(rho)),
    ggplot2::aes(y = reorder(group, rho), x = rho, fill = rho > 0)
  ) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#4878CF", "FALSE" = "#C44E52")) +
    ggplot2::labs(
      x     = "Spearman rho (peak log2FC vs gene log2FC)",
      y     = NULL,
      title = "Peak-gene FC correlation by subgroup"
    ) +
    ggplot2::theme_classic(base_size = 11)
  ggplot2::ggsave(
    file.path(CHIP_DIRS$plots, "chipseq_09_peak_gene_rho_barplot.pdf"),
    p_rho, width = 8, height = max(4, 0.4 * nrow(cor_tbl) + 2)
  )
}

msg("09 complete. Outputs:")
msg("  tables/chipseq_09_peak_signal_fc_per_peak.csv")
msg("  tables/chipseq_09_peak_gene_fc_correlation.csv")
msg("  plots/chipseq_09_peak_gene_fc_scatter_by_feature.pdf")
msg("  plots/chipseq_09_peak_gene_fc_scatter_by_quartile.pdf")
msg("  plots/chipseq_09_peak_gene_fc_scatter_tad_boundary.pdf")
msg("  plots/chipseq_09_peak_gene_rho_barplot.pdf")
