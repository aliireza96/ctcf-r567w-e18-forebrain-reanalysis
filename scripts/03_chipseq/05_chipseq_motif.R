# ============================================================
# 05_chipseq_motif.R — signal proxy plus targeted HOMER motif analysis
# ============================================================

source(file.path("scripts", "chipseq", "00_chipseq_setup.R"))
chipseq_require(c("ggplot2", "dplyr", "readr", "scales", "XML"))

class_obj_path <- file.path(CHIP_DIRS$objects, "chipseq_02_classified_peaks.rds")
if (!file.exists(class_obj_path)) {
  stop("Run 02_chipseq_peak_classification.R first: classified peak object not found.", call. = FALSE)
}
obj <- readRDS(class_obj_path)

find_homer_binary <- function(bin_name) {
  candidates <- unique(c(
    Sys.which(bin_name),
    file.path(Sys.getenv("HOMER_BIN", unset = Sys.getenv("CONDA_PREFIX")), "bin", bin_name)
  ))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  normalizePath(candidates[[1]], mustWork = TRUE)
}

run_shell_logged <- function(cmd, args, log_file) {
  quoted <- c(shQuote(cmd), vapply(args, shQuote, character(1)))
  cmd_string <- paste(quoted, collapse = " ")
  status <- system(paste0("(", cmd_string, ") >> ", shQuote(log_file), " 2>&1"))
  if (!identical(status, 0L)) {
    stop("Command failed: ", cmd_string, call. = FALSE)
  }
  invisible(status)
}

run_system2_checked <- function(cmd, args, stdout_file = NULL, stderr_file = NULL) {
  status <- system2(cmd, args = args, stdout = stdout_file, stderr = stderr_file)
  if (!identical(status, 0L)) {
    stop("Command failed: ", paste(c(cmd, args), collapse = " "), call. = FALSE)
  }
  invisible(status)
}

copy_dir_recursive <- function(src, dst) {
  if (!dir.exists(src)) {
    return(invisible(FALSE))
  }
  unlink(dst, recursive = TRUE, force = TRUE)
  dir.create(dst, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(src, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(entries) == 0) {
    return(invisible(TRUE))
  }
  ok <- file.copy(entries, dst, recursive = TRUE, overwrite = TRUE)
  invisible(all(ok))
}

write_homer_motif_file_from_sequences <- function(seqs, motif_names, out_file, major_prob = 0.97, minor_prob = 0.01) {
  stopifnot(length(seqs) == length(motif_names))
  nts <- c("A", "C", "G", "T")
  lines <- character()
  for (i in seq_along(seqs)) {
    seq_i <- toupper(seqs[[i]])
    chars <- strsplit(seq_i, "")[[1]]
    block <- c(paste0(">", seq_i, "\t", motif_names[[i]], "\t", nchar(seq_i)))
    for (nt in chars) {
      probs <- rep(minor_prob, 4L)
      idx <- match(nt, nts)
      if (!is.na(idx)) {
        probs[idx] <- major_prob
      }
      block <- c(block, paste(format(probs, digits = 3, trim = TRUE), collapse = "\t"))
    }
    lines <- c(lines, block)
  }
  writeLines(lines, out_file)
  invisible(out_file)
}

extract_homer_motif_subset <- function(known_motif_db, out_file, pattern) {
  lines <- readLines(known_motif_db)
  starts <- which(startsWith(lines, ">"))
  ends <- c(starts[-1] - 1L, length(lines))
  keep <- logical(length(starts))
  for (i in seq_along(starts)) {
    header <- lines[starts[[i]]]
    keep[[i]] <- grepl(pattern, header, ignore.case = TRUE, perl = TRUE)
  }
  out_lines <- unlist(Map(function(s, e) lines[s:e], starts[keep], ends[keep]), use.names = FALSE)
  writeLines(out_lines, out_file)
  invisible(out_file)
}

read_known_results <- function(path, scan_label) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  df <- readr::read_tsv(path, show_col_types = FALSE)
  if (nrow(df) == 0) {
    return(data.frame())
  }
  colnames(df) <- make.names(colnames(df))
  q_col <- grep("^q.value", colnames(df), value = TRUE)[1]
  target_n_col <- grep("Target.Sequences.with.Motif.of", colnames(df), value = TRUE)[1]
  target_pct_col <- grep("Target.Sequences.with.Motif$", colnames(df), value = TRUE)[1]
  background_n_col <- grep("Background.Sequences.with.Motif.of", colnames(df), value = TRUE)[1]
  background_pct_col <- grep("Background.Sequences.with.Motif$", colnames(df), value = TRUE)[1]
  out <- data.frame(
    scan_label = scan_label,
    motif_name = df$Motif.Name,
    consensus = df$Consensus,
    p_value = as.numeric(df$P.value),
    log_p_value = as.numeric(df$Log.P.value),
    q_value = as.numeric(df[[q_col]]),
    target_n = readr::parse_number(as.character(df[[target_n_col]])),
    target_pct = readr::parse_number(as.character(df[[target_pct_col]])),
    background_n = readr::parse_number(as.character(df[[background_n_col]])),
    background_pct = readr::parse_number(as.character(df[[background_pct_col]])),
    stringsAsFactors = FALSE
  )
  out$delta_pct <- out$target_pct - out$background_pct
  out
}

read_denovo_results <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  tabs <- XML::readHTMLTable(path, stringsAsFactors = FALSE)
  if (length(tabs) == 0 || nrow(tabs[[1]]) == 0) {
    return(data.frame())
  }
  df <- tabs[[1]]
  colnames(df) <- make.names(colnames(df))
  target_pct_col <- grep("Targets$", colnames(df), value = TRUE)[1]
  background_pct_col <- grep("Background$", colnames(df), value = TRUE)[1]
  positional_col <- grep("^STD.Bg.STD", colnames(df), value = TRUE)[1]
  best_match_col <- grep("^Best.Match.Details", colnames(df), value = TRUE)[1]
  best_match_raw <- as.character(df[[best_match_col]])
  out <- data.frame(
    rank = readr::parse_number(as.character(df$Rank)),
    motif_logo_consensus = as.character(df$Motif),
    p_value = as.numeric(df$P.value),
    log_p_value = as.numeric(df$log.P.pvalue),
    target_pct = readr::parse_number(as.character(df[[target_pct_col]])),
    background_pct = readr::parse_number(as.character(df[[background_pct_col]])),
    positional_bias = as.character(df[[positional_col]]),
    best_match = trimws(sub("\\(.*$", "", best_match_raw)),
    best_match_similarity = suppressWarnings(as.numeric(sub(".*\\(([0-9.]+)\\).*", "\\1", best_match_raw))),
    stringsAsFactors = FALSE
  )
  out$delta_pct <- out$target_pct - out$background_pct
  out$motif_id <- paste0("motif", out$rank)
  out
}

write_fasta_from_homer_prep <- function(seq_file, group_file, target_fasta, background_fasta) {
  seq_df <- readr::read_tsv(seq_file, col_names = c("peak_id", "sequence"), show_col_types = FALSE)
  group_df <- readr::read_tsv(group_file, col_names = c("peak_id", "group", "weight"), show_col_types = FALSE)
  merged <- dplyr::left_join(group_df, seq_df, by = "peak_id") %>%
    dplyr::filter(!is.na(sequence), nzchar(sequence))

  write_fasta <- function(df, path) {
    con <- file(path, open = "wt")
    on.exit(close(con), add = TRUE)
    for (i in seq_len(nrow(df))) {
      writeLines(paste0(">", df$peak_id[[i]]), con)
      writeLines(df$sequence[[i]], con)
    }
  }

  write_fasta(merged %>% dplyr::filter(group == 1), target_fasta)
  write_fasta(merged %>% dplyr::filter(group == 0), background_fasta)
  invisible(list(target = target_fasta, background = background_fasta))
}

build_homer_prep_from_extracts <- function(target_extract, background_extract, seq_out, group_out, target_fasta, background_fasta) {
  target_df <- readr::read_tsv(target_extract, col_names = c("peak_id", "sequence"), show_col_types = FALSE)
  background_df <- readr::read_tsv(background_extract, col_names = c("peak_id", "sequence"), show_col_types = FALSE) %>%
    dplyr::mutate(peak_id = paste0("BG_", peak_id))

  combined_df <- dplyr::bind_rows(target_df, background_df) %>%
    dplyr::filter(!is.na(sequence), nzchar(sequence))

  utils::write.table(
    combined_df,
    file = seq_out,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  group_df <- dplyr::bind_rows(
    data.frame(peak_id = target_df$peak_id, group = 1, weight = 1, stringsAsFactors = FALSE),
    data.frame(peak_id = background_df$peak_id, group = 0, weight = 1, stringsAsFactors = FALSE)
  )
  utils::write.table(
    group_df,
    file = group_out,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
  )

  write_fasta <- function(df, path) {
    con <- file(path, open = "wt")
    on.exit(close(con), add = TRUE)
    for (i in seq_len(nrow(df))) {
      writeLines(paste0(">", df$peak_id[[i]]), con)
      writeLines(df$sequence[[i]], con)
    }
  }
  write_fasta(target_df, target_fasta)
  write_fasta(background_df, background_fasta)

  invisible(list(seq = seq_out, group = group_out, target_fasta = target_fasta, background_fasta = background_fasta))
}

homer_bins <- list(
  homerTools = find_homer_binary("homerTools"),
  findKnownMotifs = find_homer_binary("findKnownMotifs.pl"),
  findMotifs = find_homer_binary("findMotifs.pl")
)
homer_available <- all(vapply(homer_bins, function(x) nzchar(x) && !is.na(x), logical(1)))

homer_root_candidates <- unique(c(
  file.path(dirname(dirname(homer_bins$findKnownMotifs))),
  normalizePath(file.path(dirname(homer_bins$findKnownMotifs), ".."), mustWork = FALSE),
  Sys.getenv("HOMER_HOME", unset = file.path(Sys.getenv("CONDA_PREFIX"), "share/homer"))
))
known_motif_candidates <- unique(file.path(
  homer_root_candidates,
  "data", "knownTFs", "vertebrates", "known.motifs"
))
known_motif_candidates <- known_motif_candidates[file.exists(known_motif_candidates)]
known_motif_db <- if (length(known_motif_candidates) > 0) known_motif_candidates[[1]] else NA_character_
genome_dir_candidates <- unique(file.path(homer_root_candidates, "data", "genomes", "mm10"))
genome_dir_candidates <- genome_dir_candidates[file.exists(genome_dir_candidates)]
homer_genome_dir <- if (length(genome_dir_candidates) > 0) genome_dir_candidates[[1]] else NA_character_
known_db_available <- homer_available && file.exists(known_motif_db) && file.exists(homer_genome_dir)

analysis_mode <- if (homer_available && known_db_available) {
  "signal_proxy_plus_targeted_homer"
} else {
  "signal_proxy_only"
}

lost_df <- granges_to_df(obj$lost_peaks) %>% dplyr::mutate(peak_class = "LOST")
maintained_df <- granges_to_df(obj$maintained_peaks) %>% dplyr::mutate(peak_class = "MAINTAINED")

wt_proxy_df <- dplyr::bind_rows(lost_df, maintained_df) %>%
  dplyr::mutate(
    wt_neg_log10_q = qValue / 10,
    peak_width = end - start + 1
  )

wilcox_signal <- suppressWarnings(stats::wilcox.test(signalValue ~ peak_class, data = wt_proxy_df))
wilcox_q <- suppressWarnings(stats::wilcox.test(wt_neg_log10_q ~ peak_class, data = wt_proxy_df))
threshold_q <- 1.3

proxy_summary <- wt_proxy_df %>%
  dplyr::group_by(peak_class) %>%
  dplyr::summarise(
    n_peaks = dplyr::n(),
    median_signalValue = stats::median(signalValue, na.rm = TRUE),
    mean_signalValue = mean(signalValue, na.rm = TRUE),
    median_neg_log10_q = stats::median(wt_neg_log10_q, na.rm = TRUE),
    frac_q_ge_1p3 = mean(wt_neg_log10_q >= threshold_q, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(analysis_mode = analysis_mode)

stats_summary <- data.frame(
  analysis_mode = analysis_mode,
  comparison = c("WT signalValue: LOST vs MAINTAINED", "WT -log10(q): LOST vs MAINTAINED"),
  p_value = c(wilcox_signal$p.value, wilcox_q$p.value),
  median_lost = c(
    stats::median(lost_df$signalValue, na.rm = TRUE),
    stats::median(lost_df$qValue / 10, na.rm = TRUE)
  ),
  median_maintained = c(
    stats::median(maintained_df$signalValue, na.rm = TRUE),
    stats::median(maintained_df$qValue / 10, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

chipseq_write_csv(proxy_summary, "chipseq_05_peak_strength_proxy_summary.csv")
chipseq_write_csv(stats_summary, "chipseq_05_signal_proxy_statistics.csv")

p_signal <- ggplot(wt_proxy_df, aes(x = peak_class, y = signalValue, fill = peak_class)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("LOST" = chipseq_palette$lost, "MAINTAINED" = chipseq_palette$maintained)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "WT signal at LOST versus MAINTAINED CTCF peaks",
    subtitle = paste0(
      "Signal-strength proxy for motif dependence. Wilcoxon p = ",
      signif(wilcox_signal$p.value, 3)
    ),
    x = NULL,
    y = "WT signalValue"
  )
chipseq_save_plot(p_signal, "chipseq_05_signal_proxy_violin", subdir = "motif", w = 7, h = 5.4)

p_q <- ggplot(wt_proxy_df, aes(x = peak_class, y = wt_neg_log10_q, fill = peak_class)) +
  geom_violin(trim = FALSE, alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white") +
  geom_hline(yintercept = threshold_q, linetype = "dashed", color = "grey45") +
  scale_fill_manual(values = c("LOST" = chipseq_palette$lost, "MAINTAINED" = chipseq_palette$maintained)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "WT confidence of LOST versus MAINTAINED peaks",
    subtitle = paste0("Dashed line marks -log10(q) = 1.3. Wilcoxon p = ", signif(wilcox_q$p.value, 3)),
    x = NULL,
    y = "WT -log10(q)"
  )
chipseq_save_plot(p_q, "chipseq_05_qvalue_proxy_violin", subdir = "motif", w = 7, h = 5.4)

p_summary <- ggplot(proxy_summary, aes(x = peak_class, y = median_signalValue, fill = peak_class)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = sprintf("%.2f", median_signalValue)), vjust = -0.25, fontface = "bold") +
  scale_fill_manual(values = c("LOST" = chipseq_palette$lost, "MAINTAINED" = chipseq_palette$maintained)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none") +
  labs(
    title = "Median WT signal at LOST versus MAINTAINED peaks",
    subtitle = paste0("Analysis mode: ", analysis_mode),
    x = NULL,
    y = "Median WT signalValue"
  )
chipseq_save_plot(p_summary, "chipseq_05_signal_proxy_bar", subdir = "motif", w = 6.5, h = 5)

known_ctcf_results <- data.frame()
known_u_results <- data.frame()
denovo_results <- data.frame()

if (homer_available && known_db_available) {
  chipseq_msg("HOMER detected — running targeted CTCF and U-motif analyses.")

  log_file <- file.path(CHIP_DIRS$logs, "05_chipseq_motif_homer.log")
  unlink(log_file, force = TRUE)
  writeLines(
    c(
      paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] 05_chipseq_motif.R HOMER log"),
      paste0("homerTools = ", homer_bins$homerTools),
      paste0("findKnownMotifs.pl = ", homer_bins$findKnownMotifs),
      paste0("findMotifs.pl = ", homer_bins$findMotifs),
      paste0("Genome dir = ", homer_genome_dir)
    ),
    con = log_file
  )

  homer_bin_dir <- dirname(homer_bins$homerTools)
  Sys.setenv(PATH = paste(unique(c(homer_bin_dir, Sys.getenv("PATH"))), collapse = .Platform$path.sep))

  lost_bed <- file.path(CHIP_DIRS$tables, "chipseq_CTCF_WT_only_LOST_peaks.bed")
  maintained_bed <- file.path(CHIP_DIRS$tables, "chipseq_CTCF_MAINTAINED_peaks.bed")
  if (!file.exists(lost_bed)) {
    export_peak_bed(obj$lost_peaks, "chipseq_CTCF_WT_only_LOST_peaks.bed")
  }
  if (!file.exists(maintained_bed)) {
    export_peak_bed(obj$maintained_peaks, "chipseq_CTCF_MAINTAINED_peaks.bed")
  }

  ctcf_motif_file <- file.path(CHIP_DIRS$objects, "chipseq_05_ctcf_like_known.motifs")
  u_motif_file <- file.path(CHIP_DIRS$objects, "chipseq_05_u_motif_like_custom.motifs")
  extract_homer_motif_subset(
    known_motif_db,
    ctcf_motif_file,
    pattern = "CTCF|BORIS|SatelliteElement"
  )
  write_homer_motif_file_from_sequences(
    seqs = c("CTGCAGTTTC", "CAGCTGTTCC", "CTGCAGTTCC"),
    motif_names = c("U_motif_IL3_like", "U_motif_Sox2_like", "U_motif_consensus_like"),
    out_file = u_motif_file
  )

  object_prep_dir <- file.path(CHIP_DIRS$objects, "chipseq_05_homer_prep")
  object_ctcf_dir <- file.path(CHIP_DIRS$objects, "chipseq_05_homer_ctcf_only")
  object_u_dir <- file.path(CHIP_DIRS$objects, "chipseq_05_homer_u_motif_only")
  object_denovo_dir <- file.path(CHIP_DIRS$objects, "chipseq_05_homer_denovo_small")

  temp_root <- tempfile(pattern = "chipseq_homer_05_", tmpdir = "/tmp")
  dir.create(temp_root, recursive = TRUE, showWarnings = FALSE)
  prep_dir <- file.path(temp_root, "prep")
  ctcf_dir <- file.path(temp_root, "ctcf_only")
  u_dir <- file.path(temp_root, "u_motif_only")
  denovo_dir <- file.path(temp_root, "denovo_small")
  temp_lost_bed <- file.path(prep_dir, "lost_peaks.bed")
  temp_maintained_bed <- file.path(prep_dir, "maintained_peaks.bed")
  temp_ctcf_motif_file <- file.path(temp_root, "ctcf_like_known.motifs")
  temp_u_motif_file <- file.path(temp_root, "u_motif_like_custom.motifs")
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(normalizePath(lost_bed, mustWork = TRUE), temp_lost_bed, overwrite = TRUE)
  file.copy(normalizePath(maintained_bed, mustWork = TRUE), temp_maintained_bed, overwrite = TRUE)
  file.copy(ctcf_motif_file, temp_ctcf_motif_file, overwrite = TRUE)
  file.copy(u_motif_file, temp_u_motif_file, overwrite = TRUE)

  prep_ready <- dir.exists(object_prep_dir) &&
    file.exists(file.path(object_prep_dir, "seq.tsv")) &&
    file.exists(file.path(object_prep_dir, "group.adj")) &&
    file.exists(file.path(object_prep_dir, "target.fa")) &&
    file.exists(file.path(object_prep_dir, "background.fa"))
  if (prep_ready) {
    chipseq_msg("Reusing cached HOMER sequence preparation.")
    copy_dir_recursive(object_prep_dir, prep_dir)
  } else {
    chipseq_msg("Preparing HOMER target/background sequence sets with homerTools extract.")
    target_extract <- file.path(prep_dir, "target_extract.tsv")
    background_extract <- file.path(prep_dir, "background_extract.tsv")
    run_system2_checked(
      homer_bins$homerTools,
      c("extract", temp_lost_bed, homer_genome_dir, "-mask"),
      stdout_file = target_extract,
      stderr_file = log_file
    )
    run_system2_checked(
      homer_bins$homerTools,
      c("extract", temp_maintained_bed, homer_genome_dir, "-mask"),
      stdout_file = background_extract,
      stderr_file = log_file
    )
    build_homer_prep_from_extracts(
      target_extract = target_extract,
      background_extract = background_extract,
      seq_out = file.path(prep_dir, "seq.tsv"),
      group_out = file.path(prep_dir, "group.adj"),
      target_fasta = file.path(prep_dir, "target.fa"),
      background_fasta = file.path(prep_dir, "background.fa")
    )
    copy_dir_recursive(prep_dir, object_prep_dir)
  }

  ctcf_ready <- dir.exists(object_ctcf_dir) && file.exists(file.path(object_ctcf_dir, "knownResults.txt"))
  if (ctcf_ready) {
    chipseq_msg("Reusing cached HOMER CTCF-like known motif scan.")
    copy_dir_recursive(object_ctcf_dir, ctcf_dir)
  } else {
    chipseq_msg("Running targeted CTCF-like known motif scan.")
    run_shell_logged(
      homer_bins$findKnownMotifs,
      c(
        "-s", file.path(prep_dir, "seq.tsv"),
        "-g", file.path(prep_dir, "group.adj"),
        "-o", ctcf_dir,
        "-m", temp_ctcf_motif_file,
        "-homer2",
        "-p", "8"
      ),
      log_file = log_file
    )
    copy_dir_recursive(ctcf_dir, object_ctcf_dir)
  }

  u_ready <- dir.exists(object_u_dir) && file.exists(file.path(object_u_dir, "knownResults.txt"))
  if (u_ready) {
    chipseq_msg("Reusing cached HOMER U-motif known motif scan.")
    copy_dir_recursive(object_u_dir, u_dir)
  } else {
    chipseq_msg("Running targeted U-motif-like known motif scan.")
    run_shell_logged(
      homer_bins$findKnownMotifs,
      c(
        "-s", file.path(prep_dir, "seq.tsv"),
        "-g", file.path(prep_dir, "group.adj"),
        "-o", u_dir,
        "-m", temp_u_motif_file,
        "-homer2",
        "-p", "8"
      ),
      log_file = log_file
    )
    copy_dir_recursive(u_dir, object_u_dir)
  }

  target_fasta <- file.path(prep_dir, "target.fa")
  background_fasta <- file.path(prep_dir, "background.fa")

  denovo_ready <- dir.exists(object_denovo_dir) && file.exists(file.path(object_denovo_dir, "homerResults.html"))
  if (denovo_ready) {
    chipseq_msg("Reusing cached HOMER de novo motif scan.")
    copy_dir_recursive(object_denovo_dir, denovo_dir)
  } else {
    chipseq_msg("Running compact HOMER de novo motif discovery (8/10/12-mers).")
    run_shell_logged(
      homer_bins$findMotifs,
      c(
        target_fasta,
        "fasta",
        denovo_dir,
        "-fastaBg", background_fasta,
        "-len", "8,10,12",
        "-S", "10",
        "-noknown",
        "-mcheck", temp_ctcf_motif_file,
        "-nogo",
        "-bits"
      ),
      log_file = log_file
    )
    copy_dir_recursive(denovo_dir, object_denovo_dir)
  }

  known_ctcf_results <- read_known_results(file.path(ctcf_dir, "knownResults.txt"), "known_ctcf_like")
  known_u_results <- read_known_results(file.path(u_dir, "knownResults.txt"), "known_u_like")
  denovo_results <- read_denovo_results(file.path(denovo_dir, "homerResults.html"))

  if (nrow(known_ctcf_results) > 0) {
    chipseq_write_csv(known_ctcf_results, "chipseq_05_motif_known_ctcf_like.csv")
  }
  if (nrow(known_u_results) > 0) {
    chipseq_write_csv(known_u_results, "chipseq_05_motif_known_u_like.csv")
  }
  if (nrow(denovo_results) > 0) {
    chipseq_write_csv(denovo_results, "chipseq_05_motif_denovo_summary.csv")
  }

  top_logo_ranks <- unique(stats::na.omit(head(denovo_results$rank, 3)))
  for (rank_i in top_logo_ranks) {
    src_logo <- file.path(denovo_dir, "homerResults", paste0("motif", rank_i, ".logo.svg"))
    dst_logo <- file.path(CHIP_DIRS$motif, paste0("chipseq_05_denovo_motif", rank_i, ".svg"))
    if (file.exists(src_logo)) {
      file.copy(src_logo, dst_logo, overwrite = TRUE)
    }
  }
}

combined_motif_parts <- list(
  proxy_summary %>%
    dplyr::transmute(
      analysis_type = "peak_strength_proxy",
      feature_name = peak_class,
      p_value = NA_real_,
      q_value = NA_real_,
      target_pct = NA_real_,
      background_pct = NA_real_,
      delta_pct = NA_real_,
      effect_value = median_signalValue,
      supporting_note = paste0("Median WT signal=", round(median_signalValue, 2))
    )
)
if (nrow(known_ctcf_results) > 0) {
  combined_motif_parts[[length(combined_motif_parts) + 1L]] <- known_ctcf_results %>%
    dplyr::transmute(
      analysis_type = scan_label,
      feature_name = motif_name,
      p_value = p_value,
      q_value = q_value,
      target_pct = target_pct,
      background_pct = background_pct,
      delta_pct = delta_pct,
      effect_value = delta_pct,
      supporting_note = consensus
    )
}
if (nrow(known_u_results) > 0) {
  combined_motif_parts[[length(combined_motif_parts) + 1L]] <- known_u_results %>%
    dplyr::transmute(
      analysis_type = scan_label,
      feature_name = motif_name,
      p_value = p_value,
      q_value = q_value,
      target_pct = target_pct,
      background_pct = background_pct,
      delta_pct = delta_pct,
      effect_value = delta_pct,
      supporting_note = consensus
    )
}
if (nrow(denovo_results) > 0) {
  combined_motif_parts[[length(combined_motif_parts) + 1L]] <- denovo_results %>%
    dplyr::transmute(
      analysis_type = "de_novo",
      feature_name = paste0("motif", rank),
      p_value = p_value,
      q_value = NA_real_,
      target_pct = target_pct,
      background_pct = background_pct,
      delta_pct = delta_pct,
      effect_value = delta_pct,
      supporting_note = best_match
    )
}
combined_motif_summary <- dplyr::bind_rows(combined_motif_parts)
chipseq_write_csv(combined_motif_summary, "chipseq_05_motif_enrichment_at_peak_classes.csv")

if (nrow(known_ctcf_results) > 0) {
  ctcf_plot_df <- rbind(
    data.frame(motif_name = known_ctcf_results$motif_name, source = "LOST peaks", pct = known_ctcf_results$target_pct),
    data.frame(motif_name = known_ctcf_results$motif_name, source = "MAINTAINED peaks", pct = known_ctcf_results$background_pct)
  )
  p_ctcf <- ggplot(ctcf_plot_df, aes(x = reorder(motif_name, pct), y = pct, fill = source)) +
    geom_col(position = position_dodge(width = 0.72), width = 0.65) +
    coord_flip() +
    scale_fill_manual(values = c("LOST peaks" = chipseq_palette$lost, "MAINTAINED peaks" = chipseq_palette$maintained)) +
    theme_classic(base_size = 12) +
    labs(
      title = "Canonical CTCF-like motifs are not enriched in LOST peaks",
      subtitle = "HOMER targeted scan using CTCF/BORIS/Satellite-like vertebrate motifs",
      x = NULL,
      y = "Percent of sequences with motif"
    )
  chipseq_save_plot(p_ctcf, "chipseq_05_ctcf_like_known_motif_scan", subdir = "motif", w = 8.6, h = 4.8)
}

if (nrow(known_u_results) > 0) {
  u_plot_df <- known_u_results %>%
    dplyr::mutate(
      label = paste0("q=", ifelse(is.na(q_value), "NA", format(signif(q_value, 2), scientific = TRUE)))
    )
  p_u <- ggplot(u_plot_df, aes(x = reorder(motif_name, delta_pct), y = delta_pct, fill = delta_pct > 0)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = label), hjust = -0.05, size = 3.4) +
    coord_flip(clip = "off") +
    scale_fill_manual(values = c("TRUE" = chipseq_palette$coral, "FALSE" = chipseq_palette$teal), guide = "none") +
    theme_classic(base_size = 12) +
    theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +
    labs(
      title = "U-motif-like sequences show a suggestive enrichment trend in LOST peaks",
      subtitle = "Primary mechanistic test for ZF11-sensitive upstream CTCF binding sequences",
      x = NULL,
      y = "LOST minus MAINTAINED motif prevalence (percentage points)"
    )
  chipseq_save_plot(p_u, "chipseq_05_u_motif_like_known_scan", subdir = "motif", w = 8.8, h = 4.6)
}

if (nrow(denovo_results) > 0) {
  denovo_plot_df <- denovo_results %>%
    dplyr::slice_min(order_by = p_value, n = min(8, nrow(denovo_results))) %>%
    dplyr::mutate(best_match_short = sub("/.*$", "", best_match))
  p_denovo <- ggplot(denovo_plot_df, aes(x = reorder(motif_id, delta_pct), y = delta_pct, fill = best_match_short)) +
    geom_col(width = 0.68) +
    geom_text(aes(label = best_match_short), hjust = -0.05, size = 3.2) +
    coord_flip(clip = "off") +
    theme_classic(base_size = 12) +
    theme(plot.margin = margin(5.5, 120, 5.5, 5.5), legend.position = "none") +
    labs(
      title = "Top de novo motifs at LOST peaks are CTCF-family-like",
      subtitle = "Compact HOMER discovery against MAINTAINED-peak background",
      x = NULL,
      y = "LOST minus MAINTAINED motif prevalence (percentage points)"
    )
  chipseq_save_plot(p_denovo, "chipseq_05_denovo_motif_summary", subdir = "motif", w = 10.5, h = 5.6)
}

chipseq_msg("Script 05 complete: signal proxy and targeted motif outputs written.")
