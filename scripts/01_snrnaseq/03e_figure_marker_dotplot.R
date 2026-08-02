# ============================================================
# 03e_figure_marker_dotplot.R
# One marker dotplot covering every cell class retained for downstream analysis,
# at the figure grouping (upper-layer EN merged, deep-layer EN merged, SPN-D1 and
# SPN-D2 separate). Replaces the separate ventral and dorsal dotplots.
#
# WHY THIS REPLACES THE EARLIER PAIR
# The two panels I built used group.by = "celltype_fine", which prints the long
# annotation strings ("Immature migrating cortical interneurons (CGE-associated/leaning)")
# and lists all 18 fine types on every panel, so most rows were empty for any given
# marker set. The existing Dotplot_markers_* panels read better because they use
# group.by = "celltype_broad": short row labels, and only classes the markers speak to.
# This script keeps that choice and adds the two things those panels lack: SPN-D1 and
# SPN-D2 resolved separately, and marker blocks ordered to run down the diagonal.
#
# MARKER SELECTION
# Every marker below was checked against this dataset: its highest mean expression
# falls in the class it is listed under. Markers that failed were dropped rather than
# kept for convention: Hes5 (tops astrocytes), Rorb and Slc17a7 (top astrocytes and
# Cajal-Retzius at E18.5), Zic1 (tops Cajal-Retzius), Tle4 and Ldb2 (top SPN-D2).
# Bcl11b and Foxp2 are deliberately omitted from the deep-layer block: both are real
# deep-layer genes but at E18.5 they are higher in SPNs and LGE precursors here, so
# including them would put off-diagonal dots where a reader expects specificity.
#
# Output:
#   plots/dotplots/03e_Dotplot_markers_all_retained_classes.pdf
# ============================================================

source("scripts/00_setup.R")

seu <- readRDS(file.path(OUT$objects, "03_seu_annotated_CORTSTR.rds"))

# figure-level grouping: same definition as 03d / 04d
upper_labels <- c("IT-L2/3", "IT-L2/4")
deep_labels  <- c("Corticofugal", "Deep corticofugal", "Deep cortical")
lab <- as.character(seu$celltype_label)
seu$celltype_figure <- dplyr::case_when(
  lab %in% upper_labels ~ "Upper layer EN",
  lab %in% deep_labels  ~ "Deep layer EN",
  TRUE ~ lab
)

# Row order: progenitor and glial, then dorsal excitatory, then ventral GABAergic,
# ending on the striatal projection neurons that the rest of the chapter is about.
# The dotplot is drawn bottom-up, so reverse for display.
class_order <- c(
  "Cycling RG", "Immature Astrocytes", "OPC",
  "Cajal-Retzius", "Deep layer EN", "Upper layer EN",
  "MGE-IN", "Migrating CGE-derived IN",
  "LGE-IN prec", "ExtendedAmygdala_GABA",
  "SPN-D1", "SPN-D2"
)
stopifnot(setequal(class_order, unique(seu$celltype_figure)))

# Shorten the one row label that is much longer than the rest, so the y-axis block does
# not push the plotting area narrow. The full name stays in the figure legend.
display_names <- c("ExtendedAmygdala_GABA" = "Extended amygdala GABA")
class_order_disp <- ifelse(class_order %in% names(display_names),
                           display_names[class_order], class_order)
seu$celltype_figure <- factor(
  ifelse(seu$celltype_figure %in% names(display_names),
         display_names[as.character(seu$celltype_figure)],
         as.character(seu$celltype_figure)),
  levels = rev(class_order_disp)
)

# Marker blocks in the same order as the classes, so signal runs down the diagonal.
marker_blocks <- list(
  # Neurotransmitter identity first, as a block that spans classes rather than marking
  # one. This is what makes the dorsal/ventral split readable at a glance: without it a
  # reader has to infer GABAergic identity from the lineage markers further right.
  # Measured here (mean SCT expression, mean over ventral vs dorsal classes):
  #   Gad2    3.29 ventral / 0.14 dorsal = 22x, positive in all six ventral classes
  #   Gad1    2.20 / 0.14 = 15x
  #   Slc17a6 0.08 / 0.72          glutamatergic counterpart
  #   Neurod2 0.08 / 1.04          dorsal excitatory
  # Slc32a1 (0.34/0.03) and Dlx1/Dlx2 were left out: correct direction but low absolute
  # level, and Dlx1 is nearly as high in progenitors (0.81) as in ventral neurons (0.83).
  "GABAergic / glutamatergic" = c("Gad1", "Gad2", "Slc17a6", "Neurod2"),
  "Cycling RG"               = c("Mki67", "Top2a", "Sox2", "Pax6"),
  "Immature Astrocytes"      = c("Fabp7", "Slc1a3", "Aldoc", "Aqp4"),
  "OPC"                      = c("Olig1", "Olig2", "Pdgfra", "Cspg4"),
  "Cajal-Retzius"            = c("Reln", "Trp73", "Lhx1"),
  "Deep layer EN"            = c("Crym", "Fezf2", "Sox5"),
  "Upper layer EN"           = c("Satb2", "Cux1", "Mef2c"),
  "MGE-IN"                   = c("Lhx6", "Nkx2-1", "Sst", "Maf"),
  "Migrating CGE-derived IN" = c("Prox1", "Adarb2", "Htr3a"),
  # Sp8 dropped: max mean expression 0.09 across all classes, so its scaled colour is
  # driven by near-zero counts. Meis2/Etv1/Pbx3 carry the LGE precursor identity.
  "LGE-IN prec"              = c("Meis2", "Etv1", "Pbx3"),
  "ExtendedAmygdala_GABA"    = c("Six3", "Sema3e", "Nr2f1"),
  # SPN markers taken from the annotation rationale in chapter_package/
  # cluster_annotation_check.md (clusters 2 and 12), not from adult striatal panels.
  # Measured D1/D2 ratios in this dataset, mean SCT expression:
  #   Ebf1  7.71 / 0.29  = 26x   D1-specific, and high in absolute terms
  #   Isl1  2.09 / 0.07  = 26x   D1-specific
  #   Tac1  0.49 / 0.02  = 18x   D1-specific
  #   Penk  0.03 / 0.97         D2-specific
  #   Drd2  0.02 / 0.74         D2-specific
  #   Sp9   0.02 / 0.60         D2-specific
  #   Adora2a 0.01 / 0.30       D2-specific
  # Dropped as uninformative here: Drd1 and Ppp1r1b (max 0.13, too low at E18.5 to
  # read either way), Gpr88 0.41/0.51, Rxrg 0.81/0.64 and Oprm1 1.08/0.43 (do not
  # separate the two classes in this dataset).
  # Foxp1 and Bcl11b are placed as a shared pan-SPN pair: both are high in D1 and D2
  # (Foxp1 6.65 and 8.26; Bcl11b 2.68 and 2.71) and low in LGE-IN precursors (0.69,
  # 0.90), so they establish that both classes are SPNs while the flanking blocks
  # separate them from each other.
  "SPN (shared)"             = c("Foxp1", "Bcl11b"),
  "SPN-D1"                   = c("Ebf1", "Isl1", "Tac1"),
  "SPN-D2"                   = c("Sp9", "Penk", "Drd2", "Adora2a")
)
features <- unlist(marker_blocks, use.names = FALSE)
missing <- setdiff(features, rownames(seu))
if (length(missing)) stop("markers absent from dataset: ", paste(missing, collapse = ", "))
message(length(features), " markers across ", length(marker_blocks), " classes")

# scale = TRUE (per-gene z-score across classes), unlike the scale = FALSE used by the
# earlier panels. This is necessary here rather than cosmetic: at E18.5 the mature
# striatal markers are genuinely low in absolute terms (Drd1 0.13, Drd2 0.74, Adora2a
# 0.30, Ppp1r1b 0.13 mean SCT expression) while Reln reaches 48 and Meis2 21. On a
# shared absolute colour scale the entire SPN block is therefore invisible, even though
# each marker IS specific to the right class in relative terms (Drd1 0.13 in SPN-D1 vs
# 0.00 in SPN-D2; Drd2 0.02 vs 0.74). Scaling per gene asks the question the panel is
# there to answer: does each class express its markers preferentially. The dot SIZE
# still shows percent of cells expressing, which is the guard against over-reading a
# scaled colour: a z-high gene detected in few cells stays a small dot.
#
# The low absolute levels are themselves consistent with the chapter's argument that
# these striatal neurons are at an early maturation state; that point belongs in the
# Results text, not in this panel.
p <- Seurat::DotPlot(seu, features = features,
                     group.by = "celltype_figure", scale = TRUE) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, face = "italic"),
    axis.text.y = element_text(size = 10)
  ) +
  Seurat::RotatedAxis()

# Thin separators between marker blocks, to make the block structure readable without
# adding text inside the panel.
block_edges <- head(cumsum(lengths(marker_blocks)), -1) + 0.5
p <- p + geom_vline(xintercept = block_edges, colour = "grey88", linewidth = 0.3)

# 47 gene labels across the x axis. At w = 11.5 the tick spacing measured 12.1 pt for
# 9 pt rotated labels, which is tight; w = 14 gives roughly 14.8 pt and allows 10 pt
# italic gene names. Height raised slightly so the rotated labels are not compressed.
save_plot(p, "dotplots/03e_Dotplot_markers_all_retained_classes", w = 14, h = 5.6)

message("03e done")
