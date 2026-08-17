# CTCF R567W E18.5 forebrain reanalysis

Analysis code for a reanalysis of the mouse and human datasets published by
Zhang et al., *Nature Communications* **15**, 5524 (2024), covering E18.5
forebrain single-nucleus RNA-seq, brain CTCF ChIP-seq, bulk RNA-seq, 4C, Hi-C,
and human cortical organoid single-cell RNA-seq.

No new data were generated. All primary data are the deposited files from that
publication; external reference datasets are GSE103983 (Mayer et al., 2018) and
GSE153164 (Di Bella et al., 2021).

## Layout

```
scripts/
  01_snrnaseq/      quality control, clustering, annotation, composition, differential expression
  02_trajectories/  external reference construction, Slingshot trajectories, pseudotime transfer
  03_chipseq/       peak classification, annotation, motif and occupancy analysis
  04_orthogonal/    bulk RNA-seq, 4C and Hi-C at the clustered protocadherin locus
  05_figures/       figure scripts, each regenerating its panels from tables/
tables/             supplementary tables S1-S24, matching the manuscript's supplementary set
docs/               methods text and a table of which script produces which output
```

The repository preserves the analysis scripts and source tables used for the installed
manuscript package. Numerical prefixes record the analysis order; not every historical stage is
currently wired into a single portable runner.

## Regenerating package-level panels

The following scripts read packaged source tables and regenerate the indicated computational
panels without rerunning upstream analyses:

```bash
python scripts/05_figures/plot_figure1D.py  . figures/
python scripts/05_figures/plot_figure2D.py  . figures/
python scripts/05_figures/plot_figure3.py   . figures/
```

`plot_figure1D.py` and `plot_figure2D.py` regenerate the statistical summary panels D of their
respective figures, and `plot_figure3.py` regenerates the Figure 3 panels.

Figure 4 and Supplementary Figure S3 were rebuilt in August 2026 after an audit of the peak-class
properties, and their current scripts are the numbered ones in `scripts/05_figures/`: `140` and
`141` are the audit and the peak-call-free promoter measurement, and `145` to `154` build the
installed panels. They read from an analysis tree rather than from `tables/` alone, so they are
provided for provenance rather than as a self-contained runner. The previous `make_figure4.py`
built panels that the audit retracted and is kept only in `_superseded_2026-08-17/`. Other installed panels retain exact source-file provenance in
the corresponding figure-folder README. This repository has not yet passed a clean-room,
end-to-end execution test.

## Upstream analysis configuration

The upstream scripts expect the deposited matrices, tracks, peak calls and external references
to be available under a configured project root. Historical absolute paths and a small number of
stage-specific setup files still require consolidation before the complete pipeline can be run
portably. The intended inputs are:

```bash
export PROJECT_ROOT=/path/to/analysis
export REF_MAYER=/path/to/GSE103983_dropseq.csv.gz
export REF_DIBELLA=/path/to/DiBella_GSE153164
```

The trajectory scripts use the Mayer and Di Bella reference atlases; ChIP-seq scripts require
the deposited normalised bigWigs and peak calls plus mm10 sequence for the motif proxy. Consult
`docs/SCRIPT_INDEX.csv` for script-to-output provenance. A future archival release should replace
the remaining absolute paths with one configuration file and record a machine-readable session
environment before claiming full end-to-end reproducibility.

## Requirements

R (>= 4.2) with Seurat (v5), Signac, slingshot, ChIPseeker,
TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db, clusterProfiler, glmGamPoi,
dplyr, ggplot2, readr, patchwork.

Python (>= 3.9) with numpy, scipy, matplotlib, pyBigWig, pyfaidx, PyMuPDF.

HOMER is required only by `03_chipseq/05_chipseq_motif.R`; set `HOMER_HOME` and
`HOMER_BIN` if it is not on the default conda path. The motif analysis used for
the published figure is `03_chipseq/11_chipseq_umotif_by_occupancy.py`, which
does not require HOMER.

## Statistical conventions

The replication structure differs by assay. The single-nucleus and CTCF ChIP-seq comparisons
each use one library per genotype, so nucleus- or site-level P values do not estimate
between-animal variation. Bulk RNA-seq and 4C have biological replication. The deposited Hi-C
maps pool two brain libraries per genotype and are used descriptively. Effect sizes accompany
tests, resampling checks accompany the composition and developmental-position summaries, and
`tables/Table_S1_composition_FDR.csv` reports the variance-inflation factor at which each
composition result would stop passing correction.

Multiple-testing correction is applied within the prespecified composition, developmental-axis
and Gene Ontology families. The 90 displayed gene-programme combinations were selected after an
observed-data audit; their permutation P and q values are retained as descriptive diagnostics,
not as an independent confirmatory family.

## Citation

Please cite the original data publication:

Zhang, J. et al. CTCF mutation at R567 causes developmental disorders via 3D genome
rearrangement and abnormal neurodevelopment. *Nat. Commun.* **15**, 5524 (2024).

## Changelog

**2026-08-04.** Added the D2 spiny-projection-neuron subclustering analysis
(`12_d2_subclustering.R`, `plot_figure2FGH.py`, `Table_S12_D2_subclustering.csv`) and the
library sex verification (`13_sex_verification.R`, `S1_sex_markers.csv`). Supplementary
table filenames were realigned to the current figure numbering: the four `S4_*` and
`S_guidance_synaptic_tests.csv` files moved to `tables/_superseded/` and their current
equivalents are the `S3_*` names. 47 scripts, 29 tables.
