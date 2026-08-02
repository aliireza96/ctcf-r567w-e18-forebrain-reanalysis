# ============================================================
# 04e_composition_descriptive_panel.R
# Composition panel for Figure 1D, presented DESCRIPTIVELY.
#
# WHY THIS SCRIPT EXISTS
# The dataset has one library per genotype (ZJ9 mutant, ZJ6 wild-type). Every
# nucleus from the same animal is a pseudoreplicate for a genotype comparison, so
# a per-cell count test does not test genotype: it tests how many cells were
# sequenced. Holding the LGE-IN precursor effect size fixed (7.51% vs 2.85%) and
# varying the number of cells gives
#     n =   50 -> p = 0.36
#     n =  250 -> p = 0.025
#     n = 1000 -> p = 4.2e-06
#     n = 8444 -> p = 1.4e-43
# The p_adj = 4.8e-47 previously printed on this panel is therefore a statement
# about sequencing depth, not about biology, and the asterisks it produced are
# misleading. They are removed here.
#
# WHAT IS KEPT
# The equal-N resampling interval is retained but relabelled. It supports one
# limited and real claim: the estimate is stable with respect to how many cells
# were sampled, i.e. it is not a small-number artefact. It is NOT a confidence
# interval for a genotype effect, because resampling cells from one animal cannot
# recover between-animal variation.
#
# Same plot idiom, colours and geometry as 04b/04d. Only the statistical
# annotation changes.
#
# Outputs:
#   plots/composition/Composition_descriptive_lollipop_celltype_figure.pdf
#   tables/04e_composition_descriptive_celltype_figure.csv
#   tables/04e_technical_variation_yardstick.csv
# ============================================================

source("scripts/00_setup.R")

analysis_scope_label <- "cortical/striatal analysis scope"
suffix <- "celltype_figure"

pooled <- readr::read_csv(
  file.path(OUT$tables, "04d_variableN_pooled_summary_celltype_figure.csv"),
  show_col_types = FALSE
)

# ---- technical-variation yardstick ------------------------------------------
# The excluded off-target regions give a measured example of how much these two
# specific libraries differ for reasons that cannot be genotype: the mutant
# library captured thalamus that the wild-type library essentially lacks, and the
# wild-type library captured cerebellum and midbrain that the mutant lacks. These
# asymmetries run in OPPOSITE directions, so they are dissection variability.
# Their magnitude bounds how confident any single-pair comparison can be.
reg <- readr::read_csv(
  file.path(OUT$tables, "03c_region_composition_by_condition.csv"),
  show_col_types = FALSE
)
yard <- reg %>%
  dplyr::select(condition, Region, frac) %>%
  tidyr::pivot_wider(names_from = condition, values_from = frac, values_fill = 0) %>%
  mutate(delta_pp = 100 * (mut - wt)) %>%
  arrange(desc(abs(delta_pp)))
write_csv(yard, file.path(OUT$tables, "04e_technical_variation_yardstick.csv"))

tech_max <- max(abs(yard$delta_pp))
bio_max  <- max(abs(pooled$mean_delta_pct))
message(sprintf("Technical (dissection) variation reaches %.1f pp; largest cell-type delta is %.1f pp.",
                tech_max, bio_max))

# ---- descriptive table -------------------------------------------------------
# Sign consistency: fraction of resamples in which the difference has the same sign as
# the mean. This is the direct measure of "is the direction stable", and it is the honest
# alternative to a p-value here: it asks whether the observed direction survives
# resampling of nuclei, without claiming anything about between-animal variation.
resamp <- readr::read_csv(
  file.path(OUT$tables, "04b_variableN_resamples_celltype_label.csv"),
  show_col_types = FALSE
)
upper_labels <- c("IT-L2/3", "IT-L2/4")
deep_labels  <- c("Corticofugal", "Deep corticofugal", "Deep cortical")
sign_tbl <- resamp %>%
  mutate(celltype = dplyr::case_when(
    celltype %in% upper_labels ~ "Upper layer EN",
    celltype %in% deep_labels  ~ "Deep layer EN",
    TRUE ~ celltype)) %>%
  group_by(N, iter, celltype) %>%
  summarise(delta = sum(delta), .groups = "drop") %>%
  group_by(celltype) %>%
  summarise(
    mean_d = mean(delta),
    sign_consistency = mean(sign(delta) == sign(mean(delta))),
    .groups = "drop"
  )

desc_tbl <- pooled %>%
  dplyr::select(celltype, mut, wt, mut_prop, wt_prop, log2fc,
                mean_delta_pct, ci_low_pct, ci_high_pct, direction) %>%
  left_join(sign_tbl %>% dplyr::select(celltype, sign_consistency), by = "celltype") %>%
  mutate(
    # Open symbols mark differences that should not be over-read, for either reason:
    #   (a) the difference is small in absolute terms (< 1 percentage point), or
    #   (b) the direction is not consistent in at least 95% of resamples.
    # These select MGE-IN, CGE-IN, Cajal-Retzius and SPN-D2; the next smallest
    # difference (OPC, 1.42 pp, 99.3% consistent) is well separated.
    confidence = ifelse(abs(mean_delta_pct) < 1.0 | sign_consistency < 0.95,
                        "small or direction not stable", "consistent"),
    # is the effect larger than the measured technical variation between these libraries?
    exceeds_technical = ifelse(abs(mean_delta_pct) >= tech_max, "yes", "no")
  ) %>%
  arrange(mean_delta_pct) %>%
  mutate(celltype = factor(celltype, levels = celltype))
write_csv(desc_tbl, file.path(OUT$tables, paste0("04e_composition_descriptive_", suffix, ".csv")))

# ---- panel: same lollipop, no significance annotation ------------------------
# Clean panel: no title, no subtitle, no caption. Everything explanatory belongs in the
# figure legend or the Results text, not inside the panel.
# Filled symbol  = difference is >= 1 pp and direction consistent in >= 95% of resamples
# Open symbol    = small (< 1 pp) or direction not stable; do not over-read
p_desc <- ggplot(desc_tbl, aes(x = mean_delta_pct, y = celltype, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_segment(aes(x = 0, xend = mean_delta_pct, y = celltype, yend = celltype), linewidth = 1.0) +
  geom_errorbarh(aes(xmin = ci_low_pct, xmax = ci_high_pct), height = 0.18, linewidth = 0.8) +
  geom_point(aes(fill = direction, shape = confidence), size = 2.8, stroke = 0.9) +
  scale_color_manual(values = c("Enriched in mut" = "#d95f02", "Enriched in wt" = "#1f77b4")) +
  scale_fill_manual(values  = c("Enriched in mut" = "#d95f02", "Enriched in wt" = "#1f77b4")) +
  scale_shape_manual(values = c("consistent" = 21, "small or direction not stable" = 21)) +
  theme_classic(12) +
  labs(x = "Difference in percentage of nuclei (mutant - wild-type)", y = NULL) +
  coord_cartesian(
    xlim = c(min(desc_tbl$ci_low_pct, na.rm = TRUE) - 0.6,
             max(desc_tbl$ci_high_pct, na.rm = TRUE) + 0.6)
  ) +
  guides(color = "none", fill = "none", shape = "none") +
  theme(plot.margin = margin(4, 8, 4, 4))

# Open symbols: overplot the flagged rows with white fill so they read as hollow.
p_desc <- p_desc +
  geom_point(
    data = dplyr::filter(desc_tbl, confidence != "consistent"),
    aes(x = mean_delta_pct, y = celltype, color = direction),
    fill = "white", shape = 21, size = 2.8, stroke = 0.9, inherit.aes = FALSE
  )
save_plot(p_desc, paste0("composition/Composition_descriptive_lollipop_", suffix), w = 6.5, h = 5)

message("04e done")
print(desc_tbl %>% dplyr::select(celltype, mean_delta_pct, ci_low_pct, ci_high_pct,
                                 sign_consistency, confidence), n = 20)
