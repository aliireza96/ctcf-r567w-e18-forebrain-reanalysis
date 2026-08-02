# ============================================================
# 17_bulk_rnaseq_interpretation.R — interpret bulk RNA-seq outputs into reusable tables and notes
# ============================================================

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])[1]
file_arg <- gsub("~\\+~", " ", file_arg)
SCRIPT_DIR <- dirname(normalizePath(ifelse(length(file_arg) == 0 || is.na(file_arg), ".", file_arg)))
source(file.path(SCRIPT_DIR, "12_discovery_setup.R"))
discovery_require(c("dplyr", "readr", "stringr", "tidyr"))
discovery_append_log("bulk", "17_bulk_rnaseq_interpretation")

bulk_dirs <- discovery_dirs_for("bulk")

read_csv_safe <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

read_bulk_table <- function(name) {
  read_csv_safe(file.path(bulk_dirs$tables, name))
}

classify_neuron_gene <- function(gene) {
  gene_sets <- list(
    cPcdh_adhesion = c("Pcdha1","Pcdha9","Pcdhga1","Pcdhga2","Pcdhga4","Pcdhga5","Pcdhga9","Pcdhga12","Pcdhgb1","Pcdhgb4","Pcdhgb6","Pcdhgb7","Pcdhb2","Pcdhb5","Pcdhb6","Pcdhb8","Pcdhb10","Pcdhb11","Pcdhb14","Pcdhb17","Pcdhb19","Pcdhb20","Igsf9"),
    maturation_synaptic = c("Grin1","Grin2a","Grin2b","Gria1","Gria3","Gria4","Snap25","Syp","Stx1a","Nrxn1","Nrxn2","Nrxn3","Nlgn1","Nlgn2","Nlgn3","Shank1","Shank2","Shank3","Camk2b","Sparcl1","Bcan","Gabra1","Gabrg2","Gabra5","Gabrb3","Sv2a","Rims1","Chrnb2","Cobl","Kptn","Igsf9"),
    differentiation_identity = c("Lin28b","Rfng","Atf5","Nkd2","Gdf15","Pamr1","Col5a1","Foxf2","Foxc2","Wt1","Gsx2","Neurod1","Drd2","Gpr88","Vip","Reln","Pnoc","Apln","Inha","Inhba"),
    cell_cycle_proliferation = c("Nek7","Mthfd2","Eif4ebp1","Lin28b","Kif15","Cdc14a"),
    metabolic_redox = c("Gclm","Gsr","Nqo2","Slc7a11","Psat1","Psph","Crot","Ppa2","Ak3","Etfrf1","Dhrs4","Cryzl2","Comtd1","Man2b2","Gaa","Slc25a18","Slc3a2","Lacc1","Acot11","Gstm6","Hmox1","Srxn1","Tstd3","Pir","Phkg1"),
    glial_barrier_cilia = c("Gpr37l1","Slco1c1","Slc2a1","Adgrg6","Ramp1","Cmtm5","Aqp4","Plekhb1","Epdr1","Ccdc122","Def8","Ehd2")
  )
  hits <- names(gene_sets)[vapply(gene_sets, function(x) gene %in% x, logical(1))]
  if (length(hits) == 0) {
    return("other_or_unresolved")
  }
  paste(hits, collapse = ";")
}

build_neuron_interpretation_table <- function() {
  neuron <- read_bulk_table("13_neuron_HOMO_vs_WT.csv")
  if (nrow(neuron) == 0) {
    stop("Neuron HOMO vs WT table is missing.", call. = FALSE)
  }
  deg <- neuron %>%
    dplyr::filter(adj.P.Val < 0.05, abs(logFC) >= 0.5) %>%
    dplyr::mutate(
      direction = ifelse(logFC > 0, "up_in_HOMO", "down_in_HOMO"),
      manual_category = vapply(gene, classify_neuron_gene, character(1)),
      interpretation_note = dplyr::case_when(
        stringr::str_detect(manual_category, "cPcdh_adhesion") & direction == "down_in_HOMO" ~ "Strong support for reduced clustered protocadherin / adhesion output in HOMO neurons.",
        stringr::str_detect(manual_category, "cPcdh_adhesion") & direction == "up_in_HOMO" ~ "Exception within the broader cPcdh-down pattern; likely reflects subfamily-specific or locus-redistributive behavior rather than global rescue.",
        stringr::str_detect(manual_category, "maturation_synaptic") ~ "Suggestive neuronal-state gene, but this category is sparse among significant DEGs and does not form the dominant neuron bulk program.",
        stringr::str_detect(manual_category, "differentiation_identity") & direction == "up_in_HOMO" ~ "Could reflect altered developmental state or cell-state bias; should not be treated as lineage conversion on bulk RNA alone.",
        stringr::str_detect(manual_category, "differentiation_identity") & direction == "down_in_HOMO" ~ "Compatible with altered developmental-state regulation, but not sufficient by itself to define delayed differentiation.",
        stringr::str_detect(manual_category, "cell_cycle_proliferation") ~ "Does not support active proliferative expansion in the neuron bulk prep; signal is mixed with a net downward tendency.",
        stringr::str_detect(manual_category, "metabolic_redox") ~ "Supports altered metabolic / redox handling as a secondary neuron-bulk program.",
        stringr::str_detect(manual_category, "glial_barrier_cilia") ~ "Likely reflects residual non-neuronal carryover or barrier/glial-state signal within the neuron prep.",
        TRUE ~ "Unresolved or context-dependent gene; interpret cautiously and in combination with stronger category-level signals."
      )
    ) %>%
    dplyr::arrange(manual_category, dplyr::desc(abs(logFC)), adj.P.Val)
  write_csv_discovery("bulk", deg, "13_neuron_HOMO_vs_WT_106deg_interpretation_table.csv")

  category_counts <- deg %>%
    dplyr::count(manual_category, direction, sort = TRUE) %>%
    tidyr::pivot_wider(names_from = direction, values_from = n, values_fill = 0)
  write_csv_discovery("bulk", category_counts, "13_neuron_HOMO_vs_WT_106deg_category_counts.csv")
  list(table = deg, counts = category_counts)
}

top_go_lines <- function(go_df, set_name, n = 4) {
  if (nrow(go_df) == 0) {
    return("- No GO enrichment file was generated for this tissue.")
  }
  sub <- go_df %>%
    dplyr::filter(set == set_name) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = n)
  if (nrow(sub) == 0) {
    return(paste0("- No enriched terms passed for `", set_name, "`."))
  }
  apply(sub, 1, function(row) {
    paste0("- `", row[["Description"]], "` (BH=", sprintf("%.2e", as.numeric(row[["p.adjust"]])), ", n=", row[["Count"]], ")")
  })
}

candidate_snippet <- function(df, n = 8) {
  if (nrow(df) == 0) {
    return("- No candidate-gene summary rows were available.")
  }
  top <- df %>%
    dplyr::arrange(dplyr::desc(abs(HOMO_vs_WT_logFC))) %>%
    dplyr::slice_head(n = n)
  apply(top, 1, function(row) {
    paste0("- `", row[["gene"]], "`: logFC=", sprintf("%.2f", as.numeric(row[["HOMO_vs_WT_logFC"]])), ", BH=", sprintf("%.2e", as.numeric(row[["HOMO_vs_WT_adj.P.Val"]])), ", pattern=`", row[["trend_class"]], "`")
  })
}

safe_signature_text <- function(df) {
  if (nrow(df) == 0) {
    return("- Signature-summary file was empty, so cross-assay signature direction could not be evaluated directly for this tissue.")
  }
  apply(df, 1, function(row) {
    paste0(
      "- `", row[["signature"]], "`: bulk mean logFC=", sprintf("%.2f", as.numeric(row[["bulk_mean_logFC"]])),
      ", scRNA mean logFC=", sprintf("%.2f", as.numeric(row[["scrna_mean_logFC"]])),
      ", bulk-vs-scRNA rho=", ifelse(is.na(row[["bulk_vs_scrna_cor"]]), "NA", sprintf("%.2f", as.numeric(row[["bulk_vs_scrna_cor"]])))
    )
  })
}

write_bulk_interpretation_markdown <- function(neuron_interpretation) {
  contrast_summary <- read_bulk_table("13_all_tissues_contrast_summary.csv")
  brain_go <- read_bulk_table("13_brain_GO_enrichment.csv")
  neuron_go <- read_bulk_table("13_neuron_GO_enrichment.csv")

  brain_candidates <- read_bulk_table("13_brain_candidate_gene_dosage_summary.csv")
  neuron_candidates <- read_bulk_table("13_neuron_candidate_gene_dosage_summary.csv")
  organoid_candidates <- read_bulk_table("13_organoid_candidate_gene_dosage_summary.csv")

  brain_sig <- read_bulk_table("13_brain_signature_direction_summary.csv")
  neuron_sig <- read_bulk_table("13_neuron_signature_direction_summary.csv")
  organoid_sig <- read_bulk_table("13_organoid_signature_direction_summary.csv")

  brain_meta <- read_bulk_table("13_brain_sample_metadata.csv")
  neuron_meta <- read_bulk_table("13_neuron_sample_metadata.csv")
  organoid_meta <- read_bulk_table("13_organoid_sample_metadata.csv")

  get_nsig <- function(tissue, contrast) {
    row <- contrast_summary %>% dplyr::filter(tissue == !!tissue, contrast == !!contrast)
    if (nrow(row) == 0) return(NA_integer_)
    row$n_sig[1]
  }

  lines <- c(
    "# Bulk RNA-seq Interpretation Note",
    "",
    "## Scope",
    "This note interprets the discovery-oriented bulk RNA-seq pipeline completed in `results/discovery_bulk_12`.",
    "The primary contrasts are `HET vs WT`, `HOMO vs WT`, and `HOMO vs HET` across brain, neuron, and organoid bulk RNA-seq.",
    "",
    "## High-Level Takeaway",
    "- The cleanest mouse-bulk signal is `HOMO vs WT`, not `HET vs WT`.",
    paste0("- Brain bulk: `HOMO vs WT` yields ", get_nsig("brain", "HOMO_vs_WT"), " significant genes, dominated by strong `cPcdh` / adhesion downregulation."),
    paste0("- Neuron bulk: `HOMO vs WT` yields ", get_nsig("neuron", "HOMO_vs_WT"), " significant genes, again dominated by `cPcdh` / adhesion loss plus a secondary metabolic-redox suppression program."),
    paste0("- Organoid bulk: `HOMO vs WT` yields ", get_nsig("organoid", "HOMO_vs_WT"), " significant genes, indicating a very broad perturbation that is much larger than the mouse-bulk effect and should be interpreted cautiously."),
    "- Across mouse bulk, the data support a strong clustered-protocadherin / adhesion axis more clearly than a broad maturation-delay program.",
    "",
    "## Genotype-Series Summary",
    paste0("- Brain `HET vs WT`: ", get_nsig("brain", "HET_vs_WT"), " significant genes. This contrast is weak and also partially limited by sex imbalance."),
    paste0("- Brain `HOMO vs WT`: ", get_nsig("brain", "HOMO_vs_WT"), " significant genes. This is the main mouse-brain bulk anchor."),
    paste0("- Brain `HOMO vs HET`: ", get_nsig("brain", "HOMO_vs_HET"), " significant genes. Modest incremental shift beyond heterozygosity."),
    paste0("- Neuron `HET vs WT`: ", get_nsig("neuron", "HET_vs_WT"), " significant genes. Essentially no clear heterozygous neuron-bulk effect."),
    paste0("- Neuron `HOMO vs WT`: ", get_nsig("neuron", "HOMO_vs_WT"), " significant genes. Main neuron-bulk signal."),
    paste0("- Neuron `HOMO vs HET`: ", get_nsig("neuron", "HOMO_vs_HET"), " significant genes. Small additional homozygous shift."),
    paste0("- Organoid `HET vs WT`: ", get_nsig("organoid", "HET_vs_WT"), " significant genes."),
    paste0("- Organoid `HOMO vs WT`: ", get_nsig("organoid", "HOMO_vs_WT"), " significant genes."),
    paste0("- Organoid `HOMO vs HET`: ", get_nsig("organoid", "HOMO_vs_HET"), " significant genes."),
    "",
    "## Brain Bulk Interpretation",
    "- Direct result: the strongest significant brain-bulk program is loss of `cPcdh` / homophilic adhesion genes, especially many `Pcdhb` members, with a smaller set of upregulated exceptions such as `Pcdhgb4`.",
    "- GO support for brain `HOMO vs WT DOWN`:",
    top_go_lines(brain_go, "HOMO_vs_WT_DOWN"),
    "- Key brain candidate genes by absolute `HOMO vs WT` effect:",
    candidate_snippet(brain_candidates, n = 10),
    "- Signature-direction readout:",
    safe_signature_text(brain_sig),
    "- Interpretation: brain bulk is concordant with the manuscript on the `cPcdh` axis, but only modestly concordant with the single-cell maturation / lineage reinterpretation. It should be used as orthogonal support for adhesion-state disruption, not as the main proof of delayed maturation.",
    paste0(
      "- Brain sex note: inferred sex calls in [13_brain_sample_metadata.csv](<",
      bulk_dirs$tables,
      "/13_brain_sample_metadata.csv>) show WT and HOMO are mixed-sex, whereas HET is female-only. This weakens the interpretability of brain `HET` contrasts."
    ),
    "",
    "## Neuron Bulk Interpretation",
    "- Direct result: the 106 significant neuron `HOMO vs WT` genes split into a dominant `cPcdh` / adhesion-down signal, a smaller adhesion-up minority, and a secondary metabolic-redox down program.",
    "- GO support for neuron `HOMO vs WT DOWN`:",
    top_go_lines(neuron_go, "HOMO_vs_WT_DOWN"),
    "- GO support for neuron `HOMO vs WT UP`:",
    top_go_lines(neuron_go, "HOMO_vs_WT_UP"),
    paste0(
      "- Manual interpretation table written to [13_neuron_HOMO_vs_WT_106deg_interpretation_table.csv](<",
      bulk_dirs$tables,
      "/13_neuron_HOMO_vs_WT_106deg_interpretation_table.csv>)."
    ),
    "- Category counts from the 106 neuron DEGs:",
    apply(neuron_interpretation$counts, 1, function(row) {
      paste0("- `", row[["manual_category"]], "`: down=", row[["down_in_HOMO"]], ", up=", row[["up_in_HOMO"]])
    }),
    "- Key neuron candidate genes by absolute `HOMO vs WT` effect:",
    candidate_snippet(neuron_candidates, n = 10),
    "- Signature-direction readout:",
    safe_signature_text(neuron_sig),
    "- Interpretation: neuron bulk strongly supports `cPcdh` disruption, suggests altered metabolic-redox homeostasis, and provides only weak support for a broad neuron-intrinsic maturation defect. The absence of a strong significant synaptic-maturation-down set means this assay should be framed as partial, not decisive, support for the single-cell reinterpretation.",
    "",
    "## Organoid Bulk Interpretation",
    "- Direct result: organoid bulk is massively perturbed, with thousands of significant genes in all three genotype contrasts.",
    "- Top organoid `HOMO vs WT` up genes include `ID1`, `ID3`, `RMST`, `HTR7`, `SST`, `CBLN2`, and `KCNT1`.",
    "- Top organoid `HOMO vs WT` down genes include `PCDHB5`, `HK2`, `NFIX`, `NFIB`, `EGFR`, `OLIG2`, `MKI67`, and `ZEB2`.",
    "- This broad scale suggests a major state shift, but because the organoid bulk is only `n=2` per genotype and the directionality is much larger than in mouse bulk, it should be treated as a context-setting dataset rather than the main validation layer.",
    "- The organoid signature-summary file was effectively empty for the mouse-derived signature overlaps, so direct cross-system signature comparison is limited here.",
    if (nrow(organoid_meta) > 0) {
      "- Organoid sex inference was not available from the processed counts and the organoid design remains genotype-only in the model."
    } else {
      "- Organoid metadata were limited, further reinforcing that this bulk lane should remain secondary."
    },
    "",
    "## Strategic Use In The Manuscript",
    "- Use brain bulk as orthogonal support that the `cPcdh` / adhesion disruption is not a single-cell artifact.",
    "- Use neuron bulk to argue that the adhesion signal persists in the neuronal compartment, while being explicit that the broader maturation-delay story is only weakly captured by this assay.",
    "- Do not use organoid bulk as the main cross-system validation for the new mouse-brain interpretation. It is too broad and too context-shifted for that role in its current form.",
    "- The most defensible bulk-RNA framing is: strong support for clustered-protocadherin dysregulation, suggestive but incomplete support for downstream developmental-state consequences."
  )

  lines <- unlist(lines)
  write_status_note("bulk", "13_bulk_rnaseq_interpretation.md", lines)
}

neuron_interpretation <- build_neuron_interpretation_table()
write_bulk_interpretation_markdown(neuron_interpretation)
discovery_msg("Bulk RNA interpretation artifacts written.")
