#!/usr/bin/env Rscript

# ============================================================
# 09_config_projection_rooted_v2.R
#
# Projection 09 v2: hand-rooted ventral references.
# Uses the inspected Mayer ventral lineage subsets and initializes
# Slingshot from the selected progenitor-like root clusters.
# Set PROJECTION_CONFIG=rooted_v2 before running 09a/09b/09c/etc.
# ============================================================

ROOT_LOCAL <- Sys.getenv("PROJECT_ROOT", unset = normalizePath(file.path(dirname(getwd()))))

CFG <- list(
  paths = list(
    ref_ventral_rds = file.path(ROOT_LOCAL, "results", "mayer_ventral_reference", "objects", "mayer_GSE103983_ventral_reference_clean.rds"),
    ref_dorsal_rds = file.path(ROOT_LOCAL, "results", "dibella_dorsal_reference", "DiBella_dorsal_reference_SCTintegrated.rds"),
    query_e18_rds = file.path(ROOT_LOCAL, "results", "objects", "03_seu_annotated_TEL.rds"),
    out_dir = file.path(ROOT_LOCAL, "results", "projection_09_rooted_v2")
  ),

  cols = list(
    genotype_col = "genotype",
    broad_col = "celltype_broad",
    query_sample_col = "sample_id",
    dorsal_label_col = "celltype_broad_dorsal",
    ventral_label_col = "region",
    stage_col = "stage"
  ),

  genotype = list(
    wt_values = c("WT", "wt", "control"),
    mut_values = c("R567W_R567W", "mut", "MUT", "r567w")
  ),

  lineages = list(
    ventral_lge_spn = list(
      reference_name = "ventral",
      ref_rds = file.path(ROOT_LOCAL, "results", "mayer_ventral_reference_lge_inspection", "objects", "07d_lge_only_annotated.rds"),
      ref_label_col = "LGE_cluster_umap_label",
      trajectory_cluster_col = "LGE_cluster",
      ref_keep_labels = c("LGE prog cycling", "LGE RG / apical prog", "vLGE iSPN"),
      root_label = "LGE prog cycling",
      root_cluster = "0",
      normalization_method = "LogNormalize",
      dims = 1:20,
      k_weight = 30,
      anchor_features = 1500L,
      query_keep_broad = c("SPNs", "LGE-IN prec", "LGE-IN precursors"),
      marker_genes = c("Gad1", "Gad2", "Isl1", "Ebf1", "Oprm1", "Foxp1", "Gpr88", "Drd2", "Adora2a", "Penk")
    ),

    ventral_mge_in = list(
      reference_name = "ventral",
      ref_rds = file.path(ROOT_LOCAL, "results", "mayer_ventral_reference_mge_inspection", "objects", "07b_mge_only_reclustered.rds"),
      ref_label_col = "MGE_cluster",
      ref_label_map = c(
        "0" = "Cycling Prog.",
        "1" = "Diff. Interneuron",
        "2" = "Trans. Prog."
      ),
      trajectory_cluster_col = "MGE_cluster",
      ref_keep_labels = c("Cycling Prog.", "Trans. Prog.", "Diff. Interneuron"),
      root_label = "Cycling Prog.",
      root_cluster = "0",
      normalization_method = "LogNormalize",
      dims = 1:20,
      k_weight = 30,
      anchor_features = 1500L,
      query_keep_broad = c("MGE-IN", "MGE-derived interneurons"),
      marker_genes = c("Gad1", "Gad2", "Lhx6", "Sox6", "Sst", "Arx", "Pvalb", "Nkx2-1", "Npy", "Nos1")
    ),

    ventral_cge_in = list(
      reference_name = "ventral",
      ref_rds = file.path(ROOT_LOCAL, "results", "mayer_ventral_reference_cge_inspection", "objects", "07c_cge_only_annotated.rds"),
      ref_label_col = "CGE_cluster_umap_label",
      trajectory_cluster_col = "CGE_cluster",
      # Keep vLGE-like iSPN cells in the CGE reference context, but root
      # the trajectory at the selected CGE progenitor cluster.
      ref_keep_labels = c("CGE prog cycling", "CGE RG / apical prog", "vLGE-like iSPN"),
      root_label = "CGE prog cycling",
      root_cluster = "0",
      normalization_method = "LogNormalize",
      dims = 1:20,
      k_weight = 30,
      anchor_features = 1500L,
      query_keep_broad = c("Maturing CGE-derived IN", "Migrating CGE-derived IN", "CGE-derived interneurons"),
      marker_genes = c("Gad1", "Gad2", "Nr2f2", "Htr3a", "Reln", "Nr3c2", "Adarb2", "Vip", "Cck", "Calb2")
    ),

    dorsal_rg_ipc_exc = list(
      reference_name = "dorsal",
      ref_label_col = "celltype_broad_dorsal",
      ref_keep_labels = c("RG", "IPC", "Excitatory"),
      root_label = "RG",
      normalization_method = "LogNormalize",
      dims = 1:10,
      k_weight = 30,
      anchor_features = 800L,
      query_keep_broad = c("Deep layer EN", "Deep layer En", "Upper layer EN", "Upper layer En", "Cycling RG", "Cycling progenitors", "Immature astrocytes", "Immature Astrocytes"),
      marker_trend_max_cells = 4000L,
      marker_genes = c("Sox2", "Pax6", "Eomes", "Tbr1", "Bcl11b", "Fezf2", "Tle4", "Satb2", "Cux1", "Cux2", "Rorb", "Neurod1", "Hes1", "Hes5", "Top2a")
    )
  ),

  mapping = list(
    normalization_method = "LogNormalize",
    dims = 1:30,
    k_weight = 50,
    anchor_features = 1500L,
    force_sct_reference_model = FALSE,
    allow_lognorm_fallback = TRUE,
    allow_manual_transfer = FALSE,
    run_mapquery_projection = FALSE
  ),

  integration_qc = list(
    do_integration_qc = FALSE,
    run_tags = c("all3", "dorsal_focus", "ventral_focus"),
    method = "harmony",
    use_all_cells = TRUE,
    harmony_use_sct = FALSE,
    harmony_theta = 2,
    harmony_nfeatures = 2000L,
    harmony_premerge_hvg_intersect = TRUE,
    harmony_hvg_per_dataset = 1500L,
    sct_ncells = 5000L,
    anchor_features = 2500L,
    anchor_dims = 1:30,
    umap_dims = 1:30,
    max_cells_per_dataset = 4000L,
    regress_percent_mt = TRUE,
    regress_cell_cycle = TRUE,
    seed = 1L
  ),

  analysis = list(
    downsample_iterations = 300L,
    late_quantile = 0.80,
    min_cells_per_genotype = 50L,
    run_mgcv = TRUE
  ),

  extra_genes_csv = NA_character_
)

cat("[config] Loaded: projection_09_rooted_v2\n")
cat("[config] Output dir:", CFG$paths$out_dir, "\n")
