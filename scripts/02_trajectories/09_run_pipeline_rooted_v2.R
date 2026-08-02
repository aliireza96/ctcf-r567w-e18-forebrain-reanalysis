#!/usr/bin/env Rscript

# ============================================================
# 09_run_pipeline_rooted_v2.R
# Rooted projection_09 v2 pipeline.
#
# Uses PROJECTION_CONFIG=rooted_v2 and writes outputs to:
#   results/projection_09_rooted_v2
# ============================================================

Sys.setenv(PROJECTION_CONFIG = "rooted_v2")

fail_fast <- TRUE

run_step <- function(path, label = NULL) {
  message("\n============================================================")
  if (!is.null(label)) message(label)
  message("Running: ", path)
  message("PROJECTION_CONFIG=", Sys.getenv("PROJECTION_CONFIG"))
  message("============================================================")
  source(path, local = FALSE)
}

steps <- data.frame(
  path = c(
    "scripts/09a_build_reference_pseudotime.R",
    "scripts/09b_map_query_to_reference.R",
    "scripts/09c_analyze_genotype_pseudotime.R",
    "scripts/09e_visualize_mapped_pseudotime.R",
    "scripts/09f_combined_lineage_umaps.R",
    "scripts/09g_pseudotime_qc_occupancy.R",
    "scripts/09i_dorsal_subtype_pseudotime.R",
    "scripts/09j_pseudotime_variableN_sensitivity.R"
  ),
  label = c(
    "Step 1/8: Build rooted reference pseudotime",
    "Step 2/8: Map E18.5 query to rooted references",
    "Step 3/8: Analyze genotype shifts",
    "Step 4/8: Plot mapped pseudotime panels",
    "Step 5/8: Combined lineage/root UMAP panels",
    "Step 6/8: Pseudotime occupancy + mapping QC",
    "Step 7/8: Dorsal subtype pseudotime summaries",
    "Step 8/8: Variable-N pseudotime sensitivity"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(steps))) {
  s <- steps$path[i]
  lbl <- steps$label[i]
  if (!file.exists(s)) stop("Missing script: ", s)
  if (isTRUE(fail_fast)) {
    run_step(s, lbl)
  } else {
    tryCatch(
      run_step(s, lbl),
      error = function(e) {
        message("Step failed: ", s)
        message("Error: ", conditionMessage(e))
      }
    )
  }
}

message("\nRooted 09 v2 pipeline completed.")
