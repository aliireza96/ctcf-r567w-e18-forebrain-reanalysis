# ============================================================
# 00_setup.R — global settings, paths, helpers
# ============================================================

set.seed(42)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(readr)
  library(scales)
})

# Project root. Override with the PROJECT_ROOT environment variable, or edit this line.
ROOT <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(file.path(dirname(getwd()))))
setwd(ROOT)

# The converted 10X folders live one level above ROOT
DATA_ROOT <- dirname(ROOT)

OUT <- list(
  objects = file.path(ROOT, "results/objects"),
  tables  = file.path(ROOT, "results/tables"),
  plots   = file.path(ROOT, "results/plots"),
  logs    = file.path(ROOT, "results/logs")
)
dir.create(file.path(ROOT, "results"), showWarnings = FALSE)
lapply(OUT, dir.create, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, fname, w=7, h=5, ext="pdf", dpi=300) {
  out_file <- file.path(OUT$plots, paste0(fname, ".", ext))
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = out_file,
         plot = p, width = w, height = h, dpi = dpi)
}

sink(file.path(OUT$logs, "sessionInfo.txt"))
print(sessionInfo())
sink()
