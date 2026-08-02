# 07e_supp_S2_reference_root_dotplots.R
# ------------------------------------------------------------------------------
# Supplementary Figure S2 panels A-D: one dotplot per reference, all four in an
# IDENTICAL format so the panels read as a set.
#
# Each panel shows the annotated subclusters of one reference (rows, ordered
# progenitor -> differentiated) against its marker set (columns, grouped into
# blocks). Block captions name the ROOT ROLE of the cluster each block
# identifies, so the panel is simultaneously a marker check and the
# justification for the trajectory root.
#
# Marker sets and cluster orders are taken verbatim from the original
# per-reference scripts (07b_II, 07c_II, 07d_II, 08b) so nothing is re-chosen
# here; only the presentation is unified.
#
# Shared format decisions (the point of this script):
#   - scale = FALSE, one legend style, one colour ramp for all four panels
#   - identical dot.scale, font sizes, panel height per row, block divider style
#   - rows ordered progenitor -> differentiated in every panel
#   - excluded clusters (LGE microglia) kept and labelled, not dropped
# ------------------------------------------------------------------------------

source("scripts/00_setup.R")
suppressMessages({library(Seurat); library(ggplot2); library(dplyr)})

IN  <- "results/"
OUTP <- file.path(OUT$plots, "supp_s2_dotplots")
dir.create(OUTP, recursive = TRUE, showWarnings = FALSE)

# ── shared style ──────────────────────────────────────────────────────────────
DOT_SCALE   <- 6.5
GENE_SIZE   <- 7.5
ROW_SIZE    <- 8.0
BLOCK_SIZE  <- 6.8
RAMP        <- c("grey92", "#3B4CC0")   # one ramp for all four panels
# Upper colour limit shared by all four panels. Set from the observed maximum
# mean expression across the four references (CGE is the highest at ~4.2), so a
# given colour means the same expression level in every panel. Values above the
# limit are squished rather than dropped.
COL_MAX     <- 4.2

find_gene <- function(genes, universe) {
  vapply(genes, function(g) {
    if (g %in% universe) return(g)
    hit <- universe[tolower(universe) == tolower(g)]
    if (length(hit)) hit[1] else NA_character_
  }, character(1))
}

# ── per-reference specification, lifted from the original scripts ─────────────
# blocks: named list of marker vectors; the NAME is the block caption and states
# the root role of the cluster that block identifies.
specs <- list(
  MGE = list(
    obj      = paste0(IN, "mayer_ventral_reference_mge_inspection/objects/07b_mge_only_reclustered.rds"),
    label_col = NULL,                        # cluster ids only; mapped below
    id_col   = "MGE_cluster",
    id_map   = c("0" = "Cycling prog. (root)", "2" = "Transitional prog.", "1" = "Diff. interneuron"),
    row_order = c("Cycling prog. (root)", "Transitional prog.", "Diff. interneuron"),
    blocks   = list(
      "pan-MGE"                 = c("Gad2", "Nkx2-1"),
      "root: cycling prog."     = c("Mki67", "Top2a", "Ascl1"),
      "non-root: transitional"  = c("Hmgb2", "Nfix", "Hes5", "Ccnd2", "Olig2"),
      "non-root: interneuron"   = c("Gad1", "Lhx6", "Lhx8")
    ),
    title = "MGE reference (Mayer et al. 2018)"
  ),
  CGE = list(
    obj       = paste0(IN, "mayer_ventral_reference_cge_inspection/objects/07c_cge_only_annotated.rds"),
    label_col = "CGE_cluster_umap_label",
    row_order = c("CGE RG / apical prog", "CGE prog cycling", "vLGE-like iSPN"),
    row_rename = c("CGE RG / apical prog" = "RG / apical prog",
                   "CGE prog cycling"     = "Cycling prog. (root)",
                   "vLGE-like iSPN"       = "Immature neuron"),
    blocks    = list(
      "shared / CGE identity"     = c("Gad1", "Gad2", "Nr2f2", "Prox1", "Nr2f1"),
      "sensitivity root: RG"      = c("Fabp7", "Hes1", "Hes5", "Slc1a3", "Sox9"),
      "root: cycling prog."       = c("Dlx1", "Dlx2", "Ascl1", "Insm1", "Top2a", "Ccnd2"),
      "non-root: immature neuron" = c("Isl1", "Foxp1", "Bcl11b", "Stmn2", "Dcx")
    ),
    title = "CGE reference (Mayer et al. 2018)"
  ),
  LGE = list(
    obj       = paste0(IN, "mayer_ventral_reference_lge_inspection/objects/07d_lge_only_annotated.rds"),
    label_col = "LGE_cluster_umap_label",
    row_order = c("LGE RG / apical prog", "LGE prog cycling", "vLGE iSPN", "Microglia"),
    row_rename = c("LGE RG / apical prog" = "RG / apical prog",
                   "LGE prog cycling"     = "Cycling prog. (root)",
                   "vLGE iSPN"            = "Immature SPN",
                   "Microglia"            = "Microglia (excluded)"),
    blocks    = list(
      "shared ventral identity"   = c("Gad1", "Gad2", "Dlx1", "Dlx2"),
      "sensitivity root: RG"      = c("Fabp7", "Hes1", "Hes5", "Ptprz1", "Nes"),
      "root: cycling prog."       = c("Ascl1", "Insm1", "Top2a", "Ccnd2"),
      "non-root: immature SPN"    = c("Ebf1", "Isl1", "Zfp503", "Foxp1", "Bcl11b"),
      "excluded: microglia"       = c("C1qa", "P2ry12", "Tyrobp", "Csf1r")
    ),
    title = "LGE reference (Mayer et al. 2018)"
  ),
  Dorsal = list(
    obj       = "results/dibella_dorsal_reference/DiBella_dorsal_reference_SCTintegrated.rds",
    label_col = "celltype_broad_dorsal",
    row_order = c("RG", "IPC", "Excitatory", "Glia"),
    row_rename = c("RG" = "Radial glia (root)", "IPC" = "IPC",
                   "Excitatory" = "Excitatory neuron", "Glia" = "Glia (excluded)"),
    blocks    = list(
      "root: radial glia"        = c("Sox2", "Pax6", "Fabp7", "Hes1", "Hes5"),
      "non-root: IPC"            = c("Eomes", "Neurod1", "Insm1", "Top2a"),
      "non-root: excitatory"     = c("Tbr1", "Bcl11b", "Fezf2", "Satb2", "Cux1", "Cux2", "Rorb"),
      "excluded: glia"           = c("Slc1a3", "Aqp4")
    ),
    title = "Dorsal reference (Di Bella et al. 2021)"
  )
)

build_panel <- function(nm, sp) {
  message("[07e] ", nm)
  seu <- readRDS(sp$obj)

  # resolve row labels
  if (is.null(sp$label_col)) {
    lab <- unname(sp$id_map[as.character(seu@meta.data[[sp$id_col]])])
  } else {
    lab <- as.character(seu@meta.data[[sp$label_col]])
    if (!is.null(sp$row_rename)) {
      hit <- lab %in% names(sp$row_rename)
      lab[hit] <- unname(sp$row_rename[lab[hit]])
    }
  }
  ro <- sp$row_order
  if (!is.null(sp$row_rename)) ro <- unname(sp$row_rename[ro])
  stopifnot(all(ro %in% unique(lab)))
  # rows top-to-bottom = progenitor -> differentiated, so reverse for ggplot y
  seu$panel_row <- factor(lab, levels = rev(ro))

  genes <- unlist(sp$blocks, use.names = FALSE)
  genes <- unique(find_gene(genes, rownames(seu)))
  genes <- genes[!is.na(genes)]

  # block boundaries for dividers and captions
  lens <- vapply(sp$blocks, function(g) sum(find_gene(g, rownames(seu)) %in% genes), integer(1))
  ends <- cumsum(lens)
  divs <- head(ends, -1) + 0.5
  mids <- ends - (lens - 1) / 2

  p <- DotPlot(seu, features = genes, group.by = "panel_row",
               scale = FALSE, dot.scale = DOT_SCALE) +
    # Shared scales across all four panels: same colour limits and the same
    # percent-expressing breaks, so dot size and colour mean the same thing in
    # every panel. Expression is left UNSCALED (scale = FALSE), so the colour
    # limit is capped at the highest value seen across the four references.
    scale_colour_gradient(low = RAMP[1], high = RAMP[2], name = "Mean expression",
                          limits = c(0, COL_MAX), oob = scales::squish,
                          breaks = c(0, 1, 2, 3, 4)) +
    scale_size(range = c(0.2, DOT_SCALE), name = "% expressing",
               limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) +
    guides(colour = guide_colourbar(order = 1, barwidth = unit(0.22, "cm"),
                                    barheight = unit(1.5, "cm")),
           size = guide_legend(order = 2)) +
    geom_vline(xintercept = divs, linetype = "dashed", colour = "grey72", linewidth = 0.32) +
    labs(title = sp$title, x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      plot.title       = element_text(size = 9.5, face = "bold", hjust = 0),
      axis.text.x      = element_text(size = GENE_SIZE, angle = 45, hjust = 1, face = "italic"),
      axis.text.y      = element_text(size = ROW_SIZE),
      axis.line        = element_line(linewidth = 0.35),
      axis.ticks       = element_line(linewidth = 0.3),
      legend.key.size  = unit(0.30, "cm"),
      legend.text      = element_text(size = 6.4),
      legend.title     = element_text(size = 6.8),
      legend.position  = "right",
      plot.margin      = margin(t = 16, r = 4, b = 2, l = 2)
    ) +
    coord_cartesian(clip = "off")

  # block captions above the plotting area
  ytop <- length(ro) + 0.62
  for (i in seq_along(sp$blocks)) {
    p <- p + annotate("text", x = mids[i], y = ytop, label = names(sp$blocks)[i],
                      size = BLOCK_SIZE / .pt, colour = "grey25", vjust = 0)
  }

  # height scales with row count so dot spacing is equal across panels
  h <- 1.05 + 0.34 * length(ro)
  ggsave(file.path(OUTP, paste0("07e_", tolower(nm), "_reference_root_dotplot.pdf")),
         p, width = 1.6 + 0.30 * length(genes), height = h, useDingbats = FALSE)
  cat(sprintf("  %s: %d rows, %d genes, %d blocks -> %.2f x %.2f in\n",
              nm, length(ro), length(genes), length(sp$blocks), 1.6 + 0.30 * length(genes), h))
  rm(seu); gc(verbose = FALSE)
  invisible(NULL)
}

for (nm in names(specs)) build_panel(nm, specs[[nm]])
message("[07e] done")
