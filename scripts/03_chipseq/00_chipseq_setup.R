# ============================================================
# 00_chipseq_setup.R — shared paths, helpers, and plotting style
# ============================================================

set.seed(42)
options(stringsAsFactors = FALSE)

E18P5_ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(file.path(dirname(getwd()))))
CHIP_RAW   <- file.path(dirname(E18P5_ROOT), "ChiP_seq_brain_reanalysis")
CHIP_OUT   <- file.path(E18P5_ROOT, "results/chipseq")

CHIP_DIRS <- list(
  root = CHIP_OUT,
  objects = file.path(CHIP_OUT, "objects"),
  tables = file.path(CHIP_OUT, "tables"),
  plots = file.path(CHIP_OUT, "plots"),
  logs = file.path(CHIP_OUT, "logs"),
  peak_overview = file.path(CHIP_OUT, "plots/peak_overview"),
  motif = file.path(CHIP_OUT, "plots/motif"),
  genomic_annotation = file.path(CHIP_OUT, "plots/genomic_annotation"),
  de_integration = file.path(CHIP_OUT, "plots/de_integration"),
  locus_tracks = file.path(CHIP_OUT, "plots/locus_tracks"),
  heatmaps = file.path(CHIP_OUT, "plots/heatmaps")
)

invisible(lapply(CHIP_DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

CHIP_FILES <- list(
  wt_peaks  = file.path(CHIP_RAW, "GSM6614265_Brain_CTCF_wt_peaks.narrowPeak.gz"),
  mut_peaks = file.path(CHIP_RAW, "GSM6614266_Brain_CTCF_homo_peaks.narrowPeak.gz"),
  wt_bw     = file.path(CHIP_RAW, "GSM6614265_Brain_CTCF_wt.mm10.bw"),
  mut_bw    = file.path(CHIP_RAW, "GSM6614266_Brain_CTCF_homo.mm10.bw")
)
stopifnot(all(vapply(CHIP_FILES, file.exists, logical(1))))

DE_FILES <- list(
  global = file.path(E18P5_ROOT, "results/tables/DE_global_allcells_LR_celltypeadjusted_mut_vs_wt.csv"),
  broad_spn = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_SPNs_mut_vs_wt.csv"),
  broad_mge = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_MGE_IN_mut_vs_wt.csv"),
  broad_deep = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_Deep_layer_EN_mut_vs_wt.csv"),
  broad_upper = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_Upper_layer_EN_mut_vs_wt.csv"),
  broad_astro = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_Immature_Astrocytes_mut_vs_wt.csv"),
  broad_lge_in = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_LGE_IN_prec_mut_vs_wt.csv"),
  broad_cge = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_Migrating_CGE_derived_IN_mut_vs_wt.csv"),
  broad_ext = file.path(E18P5_ROOT, "results/tables/DE_celltype_broad_ExtendedAmygdala_GABA_mut_vs_wt.csv"),
  fine_spn_d1 = file.path(E18P5_ROOT, "results/tables/DE_celltype_fine_Striatal_SPN_D1_striosome_like__mut_vs_wt.csv"),
  fine_spn_d2 = file.path(E18P5_ROOT, "results/tables/DE_celltype_fine_Striatal_SPN_D2_indirect_pathway__mut_vs_wt.csv"),
  fine_mge = file.path(E18P5_ROOT, "results/tables/DE_celltype_fine_MGE_derived_interneurons_SST_PV_lineage__mut_vs_wt.csv"),
  pseudo_ventral = file.path(E18P5_ROOT, "results/projection_09/results/tables/09c_ventral_lge_spn_mgcv_interaction_ranked.csv"),
  pseudo_dorsal = file.path(E18P5_ROOT, "results/projection_09/results/tables/09c_dorsal_rg_ipc_exc_mgcv_interaction_ranked.csv")
)

chipseq_palette <- list(
  navy = "#1C1C2E",
  teal = "#028090",
  coral = "#E85D5D",
  amber = "#F5A623",
  light = "#F2F4F8",
  wt = "#1f77b4",
  mut = "#d95f02",
  lost = "#E85D5D",
  gained = "#F5A623",
  maintained = "#028090"
)

chipseq_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), ..., "\n", sep = "")
}

chipseq_require <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      ". Please install them before running this chip-seq script.",
      call. = FALSE
    )
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}

chipseq_save_plot <- function(p, fname, subdir = NULL, w = 8, h = 6, ext = "pdf", dpi = 300) {
  out_dir <- if (is.null(subdir)) CHIP_DIRS$plots else file.path(CHIP_DIRS$plots, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, paste0(fname, ".", ext))
  ggplot2::ggsave(filename = out_file, plot = p, width = w, height = h, dpi = dpi)
  invisible(out_file)
}

chipseq_write_csv <- function(df, fname) {
  readr::write_csv(df, file.path(CHIP_DIRS$tables, fname))
}

chipseq_write_tsv <- function(df, fname) {
  readr::write_tsv(df, file.path(CHIP_DIRS$tables, fname))
}

chipseq_save_rds <- function(obj, fname) {
  saveRDS(obj, file.path(CHIP_DIRS$objects, fname))
}

chipseq_append_log <- function(script_name) {
  log_file <- file.path(CHIP_DIRS$logs, paste0(script_name, ".log"))
  chipseq_msg("Logging to ", log_file)
  if (sink.number(type = "output") == 0) {
    sink(log_file, split = TRUE)
  }
  if (sink.number(type = "message") == 2) {
    sink(log_file, type = "message", append = TRUE)
  }
  invisible(log_file)
}

load_narrowpeak <- function(path) {
  if (requireNamespace("rtracklayer", quietly = TRUE)) {
    gr <- tryCatch(
      rtracklayer::import(
        path,
        format = "BED",
        extraCols = c(
          signalValue = "numeric",
          pValue = "numeric",
          qValue = "numeric",
          peak = "integer"
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(gr)) {
      return(gr)
    }
  }

  col_names <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  df <- utils::read.table(gzfile(path), header = FALSE, sep = "\t", stringsAsFactors = FALSE, col.names = col_names)
  GenomicRanges::makeGRangesFromDataFrame(
    df,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )
}

standardize_peak_gr <- function(gr, sample_label = NULL) {
  if (!is.null(sample_label)) {
    gr$sample <- sample_label
  }
  if (is.null(mcols(gr)$name)) {
    gr$name <- paste0("peak_", seq_along(gr))
  }
  if (is.null(mcols(gr)$score)) {
    gr$score <- 0
  }
  if (is.null(mcols(gr)$signalValue)) {
    gr$signalValue <- NA_real_
  }
  if (is.null(mcols(gr)$pValue)) {
    gr$pValue <- NA_real_
  }
  if (is.null(mcols(gr)$qValue)) {
    gr$qValue <- NA_real_
  }
  if (is.null(mcols(gr)$peak)) {
    gr$peak <- NA_integer_
  }
  gr
}

granges_to_df <- function(gr) {
  df <- as.data.frame(gr)
  if (!"seqnames" %in% colnames(df)) {
    df$seqnames <- as.character(GenomicRanges::seqnames(gr))
  }
  if (!"start" %in% colnames(df)) {
    df$start <- BiocGenerics::start(gr)
  }
  if (!"end" %in% colnames(df)) {
    df$end <- BiocGenerics::end(gr)
  }
  if (!"width" %in% colnames(df)) {
    df$width <- BiocGenerics::width(gr)
  }
  df
}

export_peak_bed <- function(gr, fname) {
  df <- granges_to_df(gr)
  bed_df <- data.frame(
    chr = as.character(df$seqnames),
    start = df$start - 1L,
    end = df$end,
    name = dplyr::coalesce(as.character(df$name), paste0("peak_", seq_len(nrow(df)))),
    score = dplyr::coalesce(as.numeric(df$score), 0),
    strand = dplyr::coalesce(as.character(df$strand), "."),
    signalValue = dplyr::coalesce(as.numeric(df$signalValue), NA_real_),
    pValue = dplyr::coalesce(as.numeric(df$pValue), NA_real_),
    qValue = dplyr::coalesce(as.numeric(df$qValue), NA_real_),
    peak = dplyr::coalesce(as.integer(df$peak), NA_integer_),
    stringsAsFactors = FALSE
  )
  extra_cols <- setdiff(colnames(df), c("seqnames", "start", "end", "width", "strand", "name", "score", "signalValue", "pValue", "qValue", "peak"))
  if (length(extra_cols) > 0) {
    bed_df <- cbind(bed_df, df[, extra_cols, drop = FALSE])
  }
  utils::write.table(
    bed_df,
    file = file.path(CHIP_DIRS$tables, fname),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )
  invisible(bed_df)
}

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

collapse_annotation_class <- function(annotation) {
  dplyr::case_when(
    grepl("Promoter", annotation, ignore.case = TRUE) ~ "Promoter",
    grepl("5' UTR", annotation, ignore.case = TRUE) ~ "5' UTR",
    grepl("3' UTR", annotation, ignore.case = TRUE) ~ "3' UTR",
    grepl("Exon", annotation, ignore.case = TRUE) ~ "Exon",
    grepl("Intron", annotation, ignore.case = TRUE) ~ "Intron",
    grepl("Downstream", annotation, ignore.case = TRUE) ~ "Downstream",
    grepl("Intergenic", annotation, ignore.case = TRUE) ~ "Distal Intergenic",
    TRUE ~ "Other"
  )
}

safe_map_ids <- function(keys, column, keytype, multiVals = "first") {
  if (length(keys) == 0) {
    return(character())
  }
  keys <- unique(as.character(keys))
  out <- tryCatch(
    AnnotationDbi::mapIds(
      org.Mm.eg.db::org.Mm.eg.db,
      keys = keys,
      column = column,
      keytype = keytype,
      multiVals = multiVals
    ),
    error = function(e) NULL
  )
  if (!is.null(out)) {
    return(out)
  }
  fallback <- vapply(keys, function(k) {
    res <- tryCatch(
      AnnotationDbi::mapIds(
        org.Mm.eg.db::org.Mm.eg.db,
        keys = k,
        column = column,
        keytype = keytype,
        multiVals = multiVals
      ),
      error = function(e) NA_character_
    )
    as.character(res[[1]])
  }, character(1))
  stats::setNames(fallback, keys)
}

read_deg_directional <- function(path, direction = c("down", "up"), padj_cutoff = 0.05) {
  direction <- match.arg(direction)
  df <- readr::read_csv(path, show_col_types = FALSE)
  sign_ok <- if (direction == "down") df$avg_log2FC < 0 else df$avg_log2FC > 0
  unique(df$gene[!is.na(df$gene) & !is.na(df$p_val_adj) & df$p_val_adj < padj_cutoff & sign_ok])
}

read_pseudotime_hits <- function(path, fdr_cutoff = 0.05) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  unique(df$gene[!is.na(df$gene) & !is.na(df$fdr_interaction) & df$fdr_interaction < fdr_cutoff])
}

build_gene_windows <- function(txdb, upstream = 50000L, downstream = 50000L) {
  genes_gr <- GenomicFeatures::genes(txdb)
  tss_gr <- promoters(genes_gr, upstream = upstream, downstream = downstream)
  names(tss_gr) <- names(genes_gr)
  tss_gr
}

find_best_overlap_pairs <- function(query_gr, subject_gr, hits) {
  if (length(hits) == 0) {
    return(data.frame())
  }
  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)
  ov_width <- IRanges::width(IRanges::pintersect(query_gr[qh], subject_gr[sh]))
  df <- data.frame(query = qh, subject = sh, overlap_width = ov_width)
  df <- df[order(df$query, -df$overlap_width), , drop = FALSE]
  df <- df[!duplicated(df$query), , drop = FALSE]
  df
}

bw_mean_on_ranges <- function(bw_path, gr) {
  if (length(gr) == 0) {
    return(numeric())
  }
  vals <- tryCatch(
    rtracklayer::import(bw_path, which = gr, as = "NumericList"),
    error = function(e) NULL
  )
  if (is.null(vals)) {
    return(rep(NA_real_, length(gr)))
  }
  out <- vapply(vals, function(x) {
    x_num <- as.numeric(x)
    if (length(x_num) == 0) {
      return(NA_real_)
    }
    mean(x_num, na.rm = TRUE)
  }, numeric(1))
  out
}

window_overlap_flags <- function(query_gr, peak_gr, distances = c(10000L, 50000L, 100000L, 500000L)) {
  out <- list()
  for (d in distances) {
    win_gr <- GenomicRanges::resize(query_gr, width = 2L * d + 1L, fix = "center")
    out[[paste0("within_", d / 1000L, "kb")]] <- GenomicRanges::countOverlaps(win_gr, peak_gr) > 0
  }
  as.data.frame(out)
}

find_candidate_tad_file <- function() {
  candidates <- list.files(E18P5_ROOT, recursive = TRUE, full.names = TRUE,
                           pattern = "(TAD|tad).*(bed|tsv|txt)$")
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  candidates[[1]]
}

chipseq_summary_markdown <- function(lines) {
  out <- file.path(CHIP_OUT, "CHIPSEQ_INTEGRATION_SUMMARY.md")
  writeLines(lines, con = out)
  invisible(out)
}
