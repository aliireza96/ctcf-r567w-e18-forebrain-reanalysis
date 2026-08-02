# ============================================================
# 03_chipseq_annotation.R — genomic annotation of peak classes
# ============================================================

source(file.path("scripts", "chipseq", "00_chipseq_setup.R"))
chipseq_require(c(
  "GenomicRanges", "ChIPseeker", "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db", "ggplot2", "dplyr", "readr", "tidyr", "forcats"
))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
class_obj_path <- file.path(CHIP_DIRS$objects, "chipseq_02_classified_peaks.rds")
if (!file.exists(class_obj_path)) {
  stop("Run 02_chipseq_peak_classification.R first: classified peak object not found.", call. = FALSE)
}
obj <- readRDS(class_obj_path)

chipseq_msg("Annotating LOST / GAINED / MAINTAINED peaks with ChIPseeker.")
anno_lost <- ChIPseeker::annotatePeak(obj$lost_peaks, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
anno_gained <- ChIPseeker::annotatePeak(obj$gained_peaks, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")
anno_maintained <- ChIPseeker::annotatePeak(obj$maintained_peaks, tssRegion = c(-3000, 3000), TxDb = txdb, annoDb = "org.Mm.eg.db")

chipseq_save_rds(list(lost = anno_lost, gained = anno_gained, maintained = anno_maintained), "chipseq_03_peak_annotations.rds")

lost_df <- as.data.frame(anno_lost) %>% dplyr::mutate(peak_class = "LOST")
gained_df <- as.data.frame(anno_gained) %>% dplyr::mutate(peak_class = "GAINED")
maintained_df <- as.data.frame(anno_maintained) %>% dplyr::mutate(peak_class = "MAINTAINED")

chipseq_write_csv(lost_df, "chipseq_03_LOST_peaks_annotated.csv")
chipseq_write_csv(gained_df, "chipseq_03_GAINED_peaks_annotated.csv")
chipseq_write_csv(maintained_df, "chipseq_03_MAINTAINED_peaks_annotated.csv")

anno_df <- dplyr::bind_rows(lost_df, gained_df, maintained_df) %>%
  dplyr::mutate(genomic_feature = collapse_annotation_class(annotation))

feature_summary <- anno_df %>%
  dplyr::count(peak_class, genomic_feature, name = "n_peaks") %>%
  dplyr::group_by(peak_class) %>%
  dplyr::mutate(frac = n_peaks / sum(n_peaks)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    genomic_feature = forcats::fct_relevel(
      genomic_feature,
      "Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic", "Other"
    )
  )
chipseq_write_csv(feature_summary, "chipseq_03_annotation_feature_summary.csv")

p_anno <- ggplot(feature_summary, aes(x = peak_class, y = frac, fill = genomic_feature)) +
  geom_col(position = "fill", width = 0.72) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_brewer(palette = "Set3") +
  theme_classic(base_size = 13) +
  labs(
    title = "Genomic annotation of peak classes",
    subtitle = "LOST peaks can be contrasted against MAINTAINED and GAINED peaks for promoter-proximal versus distal enrichment",
    x = NULL,
    y = "Fraction of peaks",
    fill = "Genomic feature"
  )
chipseq_save_plot(p_anno, "chipseq_03_annotation_stacked_bar", subdir = "genomic_annotation", w = 9.2, h = 5.8)

p_feature_side <- ggplot(feature_summary, aes(x = genomic_feature, y = frac, fill = peak_class)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_manual(values = c(
    "LOST" = chipseq_palette$lost,
    "MAINTAINED" = chipseq_palette$maintained,
    "GAINED" = chipseq_palette$gained
  )) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(
    title = "Peak classes differ in genomic context",
    subtitle = "Side-by-side view of promoter, genic, and distal intergenic annotations",
    x = NULL,
    y = "Fraction of peaks",
    fill = NULL
  )
chipseq_save_plot(p_feature_side, "chipseq_03_annotation_grouped_bar", subdir = "genomic_annotation", w = 10.5, h = 5.6)

anno_tss <- anno_df %>%
  dplyr::transmute(
    peak_class = peak_class,
    abs_distance_to_tss = abs(distanceToTSS),
    signed_distance_to_tss = distanceToTSS
  )
chipseq_write_csv(anno_tss, "chipseq_03_distance_to_TSS.csv")

p_tss <- ggplot(anno_tss, aes(x = peak_class, y = abs_distance_to_tss, fill = peak_class)) +
  geom_violin(trim = FALSE, alpha = 0.75, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c(
    "LOST" = chipseq_palette$lost,
    "MAINTAINED" = chipseq_palette$maintained,
    "GAINED" = chipseq_palette$gained
  )) +
  scale_y_log10(labels = scales::comma_format()) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Distance to nearest TSS by peak class",
    subtitle = "Absolute TSS distance (log10 scale) highlights whether LOST peaks skew more distal than MAINTAINED peaks",
    x = NULL,
    y = "Absolute distance to TSS (bp)"
  )
chipseq_save_plot(p_tss, "chipseq_03_distance_to_TSS_violin", subdir = "genomic_annotation", w = 7.6, h = 5.4)

nearest_gene_summary <- anno_df %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::count(peak_class, SYMBOL, sort = TRUE, name = "n_peaks")
chipseq_write_csv(nearest_gene_summary, "chipseq_03_nearest_gene_peak_counts.csv")

chipseq_msg("Script 03 complete: annotated peak tables and genomic context plots written.")
