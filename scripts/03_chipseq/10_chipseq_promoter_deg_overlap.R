# ============================================================
# 10_chipseq_promoter_deg_overlap.R
# Promoter-proximal CTCF peak turnover versus global snRNA-seq DE
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tidyr)
})

ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(file.path(dirname(getwd()))))
CHIP_TABLE_DIR <- file.path(ROOT, "results/chipseq/tables")
CHIP_PLOT_DIR <- file.path(ROOT, "results/chipseq/plots/de_integration")
dir.create(CHIP_PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

DE_FILE <- file.path(ROOT, "results/tables/05b_cortstr_plus_hypoglut_DE_global_allcells_LR_celltypeadjusted_mut_vs_wt.csv")
if (!file.exists(DE_FILE)) {
  DE_FILE <- file.path(ROOT, "results/tables/DE_global_allcells_LR_celltypeadjusted_mut_vs_wt.csv")
}
stopifnot(file.exists(DE_FILE))

COL_DARK <- "#1F2937"
COL_GRID <- "#E5E7EB"
COL_LIGHT <- "#F3F4F6"
COL_DOWN <- "#177E89"
COL_UP <- "#BA4A3C"
COL_LOST <- "#E85D5D"
COL_GAINED <- "#F5A623"
COL_MAINTAINED <- "#028090"

save_plot_dual <- function(p, name, w = 7, h = 4.5, dpi = 450) {
  ggsave(file.path(CHIP_PLOT_DIR, paste0(name, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(CHIP_PLOT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = dpi, bg = "white")
}

theme_chip_slide <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = COL_DARK, size = base_size + 3),
      plot.subtitle = element_text(color = "#4B5563", size = base_size - 1),
      axis.title = element_text(color = COL_DARK),
      axis.text = element_text(color = COL_DARK),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_text(color = COL_DARK, face = "bold"),
      legend.text = element_text(color = COL_DARK),
      strip.background = element_rect(fill = "white", color = COL_DARK, linewidth = 0.45),
      strip.text = element_text(face = "bold", color = COL_DARK),
      panel.grid = element_blank()
    )
}

fmt_p <- function(x) {
  case_when(
    is.na(x) ~ "NA",
    x < 0.001 ~ "<0.001",
    x < 0.01 ~ sprintf("%.3f", x),
    TRUE ~ sprintf("%.2f", x)
  )
}

read_annotated_peaks <- function(cls) {
  read_csv(
    file.path(CHIP_TABLE_DIR, paste0("chipseq_03_", cls, "_peaks_annotated.csv")),
    show_col_types = FALSE
  ) %>%
    mutate(event_class = cls)
}

peak_promoters <- bind_rows(lapply(c("LOST", "GAINED", "MAINTAINED"), read_annotated_peaks)) %>%
  mutate(
    gene = SYMBOL,
    is_promoter = str_detect(annotation, regex("Promoter", ignore_case = TRUE))
  ) %>%
  filter(is_promoter, !is.na(gene), gene != "")

promoter_gene_counts <- peak_promoters %>%
  count(gene, event_class, name = "n_promoter_peaks") %>%
  pivot_wider(
    names_from = event_class,
    values_from = n_promoter_peaks,
    values_fill = 0,
    names_prefix = "n_promoter_"
  )

promoter_gene_flags <- peak_promoters %>%
  distinct(gene, event_class) %>%
  mutate(has_event = TRUE) %>%
  pivot_wider(
    names_from = event_class,
    values_from = has_event,
    values_fill = FALSE,
    names_prefix = "has_promoter_"
  )

de <- read_csv(DE_FILE, show_col_types = FALSE) %>%
  filter(!is.na(gene), gene != "") %>%
  mutate(
    de_direction = case_when(
      p_val_adj < 0.05 & avg_log2FC < 0 ~ "RNA down in MUT",
      p_val_adj < 0.05 & avg_log2FC > 0 ~ "RNA up in MUT",
      TRUE ~ "Not FDR-significant"
    )
  )

gene_df <- de %>%
  left_join(promoter_gene_counts, by = "gene") %>%
  left_join(promoter_gene_flags, by = "gene") %>%
  mutate(
    across(starts_with("n_promoter_"), ~ replace_na(.x, 0L)),
    across(starts_with("has_promoter_"), ~ replace_na(.x, FALSE)),
    has_any_promoter_ctcf = has_promoter_LOST | has_promoter_GAINED | has_promoter_MAINTAINED,
    promoter_event_call = case_when(
      has_promoter_LOST & !has_promoter_GAINED ~ "LOST promoter",
      !has_promoter_LOST & has_promoter_GAINED ~ "GAINED promoter",
      has_promoter_LOST & has_promoter_GAINED ~ "Mixed LOST+GAINED promoter",
      has_promoter_MAINTAINED ~ "Maintained promoter only",
      TRUE ~ "No promoter CTCF peak"
    )
  )

write_csv(gene_df, file.path(CHIP_TABLE_DIR, "chipseq_10_promoter_event_gene_flags.csv"))

run_one_fisher <- function(df, event_col, direction_label, background_label) {
  event_present <- df[[event_col]]
  target <- df$de_direction == direction_label
  a <- sum(event_present & target, na.rm = TRUE)
  b <- sum(event_present & !target, na.rm = TRUE)
  c <- sum(!event_present & target, na.rm = TRUE)
  d <- sum(!event_present & !target, na.rm = TRUE)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "greater")
  tibble(
    background = background_label,
    event = recode(event_col, has_promoter_LOST = "LOST promoter", has_promoter_GAINED = "GAINED promoter"),
    rna_direction = direction_label,
    n_background = nrow(df),
    n_event_genes = sum(event_present, na.rm = TRUE),
    n_direction_genes = sum(target, na.rm = TRUE),
    n_overlap = a,
    expected_overlap = sum(event_present, na.rm = TRUE) * sum(target, na.rm = TRUE) / nrow(df),
    pct_event_genes_in_direction = a / max(sum(event_present, na.rm = TRUE), 1),
    odds_ratio = unname(ft$estimate),
    conf_low = ft$conf.int[[1]],
    conf_high = ft$conf.int[[2]],
    p_value = ft$p.value
  )
}

backgrounds <- list(
  `All tested genes` = gene_df,
  `Genes with any promoter CTCF peak` = filter(gene_df, has_any_promoter_ctcf)
)

overlap_summary <- bind_rows(lapply(names(backgrounds), function(background_name) {
  bg_df <- backgrounds[[background_name]]
  bind_rows(lapply(c("has_promoter_LOST", "has_promoter_GAINED"), function(event_col) {
    bind_rows(lapply(c("RNA down in MUT", "RNA up in MUT"), function(direction_label) {
      run_one_fisher(bg_df, event_col, direction_label, background_name)
    }))
  }))
})) %>%
  group_by(background) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    log2_or = log2(pmax(odds_ratio, 1e-6)),
    event = factor(event, levels = c("LOST promoter", "GAINED promoter")),
    rna_direction = factor(rna_direction, levels = c("RNA down in MUT", "RNA up in MUT")),
    label = paste0(
      "n=", n_overlap,
      "\nOR=", sprintf("%.2f", odds_ratio),
      "\nq=", fmt_p(p_adj)
    )
  )

write_csv(overlap_summary, file.path(CHIP_TABLE_DIR, "chipseq_10_promoter_event_DE_overlap_summary.csv"))

heatmap_df <- overlap_summary %>%
  filter(background == "Genes with any promoter CTCF peak") %>%
  mutate(log2_or_plot = pmax(pmin(log2_or, 1.2), -1.2))

p_overlap <- ggplot(heatmap_df, aes(x = rna_direction, y = event, fill = log2_or_plot)) +
  geom_tile(color = "white", linewidth = 1.4, width = 0.96, height = 0.92) +
  geom_text(aes(label = label), color = COL_DARK, size = 4.1, lineheight = 0.92, fontface = "bold") +
  scale_fill_gradient2(
    low = "#D9DEE7",
    mid = "white",
    high = "#263A78",
    midpoint = 0,
    limits = c(-1.2, 1.2),
    name = "log2(OR)"
  ) +
  coord_equal() +
  theme_chip_slide(base_size = 13) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(8, 12, 8, 8)
  ) +
  labs(
    title = "Promoter overlap: directional, not global",
    subtitle = "Background: genes with promoter CTCF; DE: cell-type-adjusted global RNA",
    x = NULL,
    y = NULL
  )

save_plot_dual(p_overlap, "chipseq_10_promoter_DE_overlap_heatmap", w = 6.7, h = 3.9)

p_overlap_all_bg <- overlap_summary %>%
  mutate(
    background = factor(background, levels = c("All tested genes", "Genes with any promoter CTCF peak")),
    log2_or_plot = pmax(pmin(log2_or, 1.2), -1.2)
  ) %>%
  ggplot(aes(x = rna_direction, y = event, fill = log2_or_plot)) +
  geom_tile(color = "white", linewidth = 1.1, width = 0.96, height = 0.92) +
  geom_text(aes(label = label), color = COL_DARK, size = 3.35, lineheight = 0.92, fontface = "bold") +
  facet_wrap(~background, nrow = 1) +
  scale_fill_gradient2(
    low = "#D9DEE7",
    mid = "white",
    high = "#263A78",
    midpoint = 0,
    limits = c(-1.2, 1.2),
    name = "log2(OR)"
  ) +
  coord_equal() +
  theme_chip_slide(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", angle = 18, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Promoter CTCF turnover versus RNA direction",
    subtitle = "FDR-significant global DE, mutant versus wild type",
    x = NULL,
    y = NULL
  )

save_plot_dual(p_overlap_all_bg, "chipseq_10_promoter_DE_overlap_heatmap_with_background_sensitivity", w = 9.8, h = 4.4)

candidate_pool <- gene_df %>%
  filter(
    (n_promoter_LOST > 0 & de_direction == "RNA down in MUT") |
      (n_promoter_GAINED > 0 & de_direction == "RNA up in MUT")
  ) %>%
  mutate(
    candidate_class = case_when(
      n_promoter_LOST > 0 & de_direction == "RNA down in MUT" ~ "Promoter LOST + RNA down",
      n_promoter_GAINED > 0 & de_direction == "RNA up in MUT" ~ "Promoter GAINED + RNA up"
    ),
    abs_log2FC = abs(avg_log2FC)
  )

meaningful_gene_notes <- tibble::tribble(
  ~gene, ~selection_reason,
  "Atp1a3", "neuronal ion pump / activity; promoter LOST + RNA down",
  "Ptprn", "secretory-vesicle / neuronal maturation; promoter LOST + RNA down",
  "Dpp10", "neuronal excitability / potassium-channel modulation; promoter LOST + RNA down",
  "Rimbp2", "presynaptic active-zone scaffold; promoter LOST + RNA down",
  "Caskin1", "synaptic scaffold; promoter LOST + RNA down",
  "Sox11", "immature-neuron transcription factor; promoter GAINED + RNA up",
  "Meis2", "LGE / ventral identity transcription factor; promoter GAINED + RNA up",
  "Bcl11b", "SPN / neuronal identity factor; promoter GAINED + RNA up",
  "Dlx6os1", "DLX-locus ventral programme; promoter GAINED + RNA up",
  "Bcl11a", "neurodevelopmental transcription factor; multiple promoter GAINED peaks"
)

candidate_slide <- candidate_pool %>%
  inner_join(meaningful_gene_notes, by = "gene") %>%
  mutate(
    candidate_class = factor(candidate_class, levels = c("Promoter LOST + RNA down", "Promoter GAINED + RNA up")),
    gene_order = match(gene, meaningful_gene_notes$gene)
  ) %>%
  arrange(gene_order)

write_csv(candidate_pool %>% arrange(candidate_class, p_val_adj), file.path(CHIP_TABLE_DIR, "chipseq_10_promoter_candidate_genes_all.csv"))
write_csv(candidate_slide, file.path(CHIP_TABLE_DIR, "chipseq_10_promoter_candidate_genes_slide.csv"))

candidate_long <- candidate_slide %>%
  mutate(
    gene = factor(gene, levels = rev(candidate_slide$gene)),
    rna_label = sprintf("%+.2f", avg_log2FC),
    lost_label = if_else(n_promoter_LOST > 0, as.character(n_promoter_LOST), ""),
    gained_label = if_else(n_promoter_GAINED > 0, as.character(n_promoter_GAINED), "")
  ) %>%
  select(gene, candidate_class, avg_log2FC, n_promoter_LOST, n_promoter_GAINED, rna_label, lost_label, gained_label) %>%
  pivot_longer(
    cols = c(avg_log2FC, n_promoter_LOST, n_promoter_GAINED),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("avg_log2FC", "n_promoter_LOST", "n_promoter_GAINED"),
      labels = c("RNA\nlog2FC", "LOST\npromoter", "GAINED\npromoter")
    ),
    label = case_when(
      metric == "RNA\nlog2FC" ~ rna_label,
      metric == "LOST\npromoter" ~ lost_label,
      metric == "GAINED\npromoter" ~ gained_label,
      TRUE ~ ""
    ),
    tile_fill = case_when(
      metric == "RNA\nlog2FC" & value < 0 ~ COL_DOWN,
      metric == "RNA\nlog2FC" & value > 0 ~ COL_UP,
      metric == "LOST\npromoter" & value > 0 ~ COL_LOST,
      metric == "GAINED\npromoter" & value > 0 ~ COL_GAINED,
      TRUE ~ COL_LIGHT
    ),
    tile_alpha = case_when(
      metric == "RNA\nlog2FC" ~ pmin(abs(value) / 1.1, 1),
      metric != "RNA\nlog2FC" & value > 0 ~ pmin(value / 4, 1),
      TRUE ~ 0.35
    )
  )

p_candidates <- ggplot(candidate_long, aes(x = metric, y = gene)) +
  geom_tile(aes(fill = tile_fill, alpha = tile_alpha), color = "white", linewidth = 0.75, width = 0.96, height = 0.88) +
  geom_text(aes(label = label), color = COL_DARK, size = 3.55, fontface = "bold") +
  facet_grid(candidate_class ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_identity() +
  scale_alpha(range = c(0.18, 0.95), guide = "none") +
  theme_chip_slide(base_size = 12) +
  theme(
    axis.text.x = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(face = "bold", size = 10.5),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 9.7, face = "bold"),
    strip.background = element_rect(fill = "white", color = "white"),
    plot.margin = margin(8, 10, 8, 8)
  ) +
  labs(
    title = "Promoter-proximal examples",
    subtitle = "Numbers show promoter peak counts; RNA log2FC is MUT - WT",
    x = NULL,
    y = NULL
  )

save_plot_dual(p_candidates, "chipseq_10_promoter_candidate_gene_heatmap", w = 6.6, h = 4.9)

candidate_inset_genes <- meaningful_gene_notes$gene
candidate_inset_long <- candidate_long %>%
  filter(as.character(gene) %in% candidate_inset_genes) %>%
  mutate(
    gene = factor(as.character(gene), levels = rev(candidate_inset_genes)),
    candidate_class = factor(
      candidate_class,
      levels = c("Promoter LOST + RNA down", "Promoter GAINED + RNA up"),
      labels = c("LOST + down", "GAINED + up")
    )
  )

p_candidates_inset <- ggplot(candidate_inset_long, aes(x = metric, y = gene)) +
  geom_tile(aes(fill = tile_fill, alpha = tile_alpha), color = "white", linewidth = 0.65, width = 0.96, height = 0.88) +
  geom_text(aes(label = label), color = COL_DARK, size = 3.0, fontface = "bold") +
  facet_grid(candidate_class ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_identity() +
  scale_alpha(range = c(0.18, 0.95), guide = "none") +
  theme_chip_slide(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "bold", size = 8.8),
    axis.text.y = element_text(face = "bold", size = 9.2),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 0.5, size = 8.5, face = "bold"),
    strip.background = element_rect(fill = "white", color = "white"),
    plot.title = element_text(face = "bold", color = COL_DARK, size = 11.5),
    plot.subtitle = element_blank(),
    plot.margin = margin(2, 4, 2, 2)
  ) +
  labs(
    title = "Promoter examples",
    x = NULL,
    y = NULL
  )

save_plot_dual(p_candidates_inset, "chipseq_10_promoter_candidate_gene_heatmap_inset", w = 4.95, h = 4.1)

p_overlap_inset <- p_overlap +
  labs(title = NULL, subtitle = NULL) +
  theme(
    axis.text.x = element_text(face = "bold", size = 9.5),
    axis.text.y = element_text(face = "bold", size = 9.5),
    legend.position = "none",
    plot.margin = margin(2, 4, 2, 2)
  )

save_plot_dual(p_overlap_inset, "chipseq_10_promoter_DE_overlap_heatmap_inset", w = 3.9, h = 2.65)

count_inset_df <- overlap_summary %>%
  filter(
    background == "Genes with any promoter CTCF peak",
    (event == "LOST promoter" & rna_direction == "RNA down in MUT") |
      (event == "GAINED promoter" & rna_direction == "RNA up in MUT")
  ) %>%
  mutate(
    group = factor(
      case_when(
        event == "LOST promoter" ~ "LOST promoter\n+ RNA down",
        TRUE ~ "GAINED promoter\n+ RNA up"
      ),
      levels = c("LOST promoter\n+ RNA down", "GAINED promoter\n+ RNA up")
    ),
    fill_col = if_else(event == "LOST promoter", COL_LOST, COL_GAINED),
    count_label = paste0(
      n_overlap, "\n",
      "genes"
    )
  )

p_counts_inset <- ggplot(count_inset_df, aes(x = group, y = 1)) +
  geom_tile(aes(fill = fill_col), color = "white", linewidth = 0.8, width = 0.94, height = 0.88, alpha = 0.48) +
  geom_text(aes(label = count_label), color = COL_DARK, size = 3.2, lineheight = 0.92, fontface = "bold") +
  scale_fill_identity() +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 10) +
  theme(
    axis.text.x = element_text(color = COL_DARK, face = "bold", size = 9.3),
    plot.title = element_text(face = "bold", color = COL_DARK, size = 11.5, hjust = 0),
    plot.margin = margin(2, 4, 2, 2)
  ) +
  scale_x_discrete(position = "top") +
  labs(title = "Directional promoter-overlap candidates")

p_counts_plus_candidates <- p_counts_inset / p_candidates_inset +
  plot_layout(heights = c(0.48, 1.52))

save_plot_dual(
  p_counts_plus_candidates,
  "chipseq_10_promoter_counts_plus_meaningful_genes_inset",
  w = 5.05,
  h = 5.7
)

p_slidebody <- p_overlap / p_candidates +
  plot_layout(heights = c(0.72, 1.28)) +
  plot_annotation(
    title = "Promoter CTCF turnover provides candidate RNA links",
    theme = theme(
      plot.title = element_text(face = "bold", color = COL_DARK, size = 18),
      plot.margin = margin(6, 8, 6, 8)
    )
  )

save_plot_dual(p_slidebody, "chipseq_10_promoter_overlap_plus_candidates_slidebody", w = 7.1, h = 9.2)

cat("Wrote promoter-overlap ChIP/RNA outputs to:\n")
cat("  ", file.path(CHIP_TABLE_DIR, "chipseq_10_promoter_event_DE_overlap_summary.csv"), "\n", sep = "")
cat("  ", file.path(CHIP_PLOT_DIR, "chipseq_10_promoter_DE_overlap_heatmap.pdf"), "\n", sep = "")
cat("  ", file.path(CHIP_PLOT_DIR, "chipseq_10_promoter_candidate_gene_heatmap.pdf"), "\n", sep = "")
cat("  ", file.path(CHIP_PLOT_DIR, "chipseq_10_promoter_overlap_plus_candidates_slidebody.pdf"), "\n", sep = "")
