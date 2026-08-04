#!/usr/bin/env Rscript
# Sex verification of the two E18.5 libraries.
#
# Establishes the sex of each profiled embryo from the X-inactivation transcript and
# Y-linked gene detection, so that a sex difference can be excluded as a source of
# the between-library composition and maturation differences.
#
# Writes S1_sex_markers.csv, read by the supplementary panel script.
#
# Usage:  Rscript 13_sex_verification.R [out_dir]

suppressMessages({
  library(Seurat)
  library(Matrix)
})

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "."
root    <- Sys.getenv("PROJECT_ROOT", unset = NA)
if (is.na(root)) stop("set PROJECT_ROOT to the analysis root directory")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

X_MARKER <- "Xist"
Y_GENES  <- c("Ddx3y", "Uty", "Eif2s3y", "Kdm5d")

seu <- readRDS(file.path(root, "E18p5_clean", "results", "objects", "03_seu_annotated.rds"))

# Detection rates are read from RNA log-normalised values, not SCT residuals, so that
# a zero means the transcript was not captured in that nucleus. The saved object
# defaults to SCT, so set this explicitly.
DefaultAssay(seu) <- "RNA"

# In Seurat v5 the RNA assay is stored as one layer per sample. Join them so that a
# single data matrix spanning both libraries can be read.
if (length(Layers(seu[["RNA"]], search = "data")) > 1)
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

md   <- seu@meta.data
expr <- GetAssayData(seu, layer = "data")

geno <- ifelse(md$condition %in% c("wt", "WT"), "WT", "MUT")
y_present <- intersect(Y_GENES, rownames(expr))
if (!X_MARKER %in% rownames(expr)) stop(sprintf("%s not in the object", X_MARKER))
message(sprintf("Y genes found: %s", paste(y_present, collapse = ", ")))

xist <- as.numeric(expr[X_MARKER, ])
Ymat <- if (length(y_present))
  as.matrix(expr[y_present, , drop = FALSE]) else matrix(0, 1, ncol(expr))
# Per-cell detection uses the maximum across Y genes, so a nucleus counts as Y-positive
# if any of them is detected. mean_Y is the grand mean over genes and cells, which is
# the quantity plotted in the panel.
ymax <- apply(Ymat, 2, max)

per_cell <- data.frame(cell = colnames(seu), geno = geno,
                       Xist = xist, Ymax = ymax)
write.csv(per_cell, file.path(out_dir, "d1_sex_markers_percell.csv"), row.names = FALSE)

LABEL <- c(WT = "wild type", MUT = "R567W homozygous")
out <- NULL
for (g in c("WT", "MUT")) {
  k <- geno == g
  out <- rbind(out, data.frame(
    library           = LABEL[[g]],
    n_nuclei          = sum(k),
    mean_Xist         = round(mean(xist[k]), 3),
    pct_Xist_positive = round(100 * mean(xist[k] > 0), 1),
    mean_Y            = formatC(mean(Ymat[, k, drop = FALSE]), format = "f", digits = 5),
    pct_Y_positive    = round(100 * mean(ymax[k] > 0), 2),
    sex_call          = ifelse(mean(xist[k] > 0) > 0.5 && mean(ymax[k] > 0) < 0.05,
                               "female", "check manually")))
}

hdr <- c(
  "# Sex verification of the two E18.5 libraries.",
  sprintf("# X-inactivation transcript: %s. Y-linked genes: %s.",
          X_MARKER, paste(y_present, collapse = ", ")),
  "# pct_Y_positive is the percentage of nuclei detecting ANY of the Y genes.",
  "# mean_Y is the grand mean of Y-gene expression over genes and nuclei.",
  "# A female call requires Xist in most nuclei and Y detection under 5%.")
f <- file.path(out_dir, "S1_sex_markers.csv")
writeLines(hdr, f)
suppressWarnings(write.table(out, f, sep = ",", row.names = FALSE, append = TRUE))

print(out, row.names = FALSE)
message(sprintf("\nwrote %s", f))
