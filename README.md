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
tables/             supplementary tables S1-S10 and supporting result tables
docs/               methods text and a table of which script produces which output
```

42 scripts. The pipeline stages run in numerical order within each directory.

## Reproducing the figures

The figure scripts read only from `tables/` and regenerate the published panels
without rerunning the analysis pipeline:

```bash
python scripts/05_figures/plot_figure1D.py  . figures/
python scripts/05_figures/plot_figure2D.py  . figures/
python scripts/05_figures/plot_figure3.py   . figures/
python scripts/05_figures/make_figure4.py   . figures/
```

Each was verified to reproduce its installed panel exactly.

## Reproducing the analysis

Set the project root and reference locations, then run each stage in order:

```bash
export PROJECT_ROOT=/path/to/analysis
export REF_MAYER=/path/to/GSE103983_dropseq.csv.gz
export REF_DIBELLA=/path/to/DiBella_GSE153164
```

```r
source("scripts/01_snrnaseq/00_setup.R")
source("scripts/01_snrnaseq/01_load_qc_filter_doublets.R")
# ... in numerical order
```

The trajectory stage is driven by `02_trajectories/09_run_pipeline_rooted_v2.R`,
which reads its configuration from `09_config_projection_rooted_v2.R`. The
ChIP-seq stage is driven by the numbered scripts in `03_chipseq/` in order.

Running the full pipeline from raw counts takes several hours and requires the
deposited FASTQ-derived matrices, bigWig tracks and peak calls, plus an mm10
genome FASTA for the motif analysis.

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

Every genotype comparison rests on one library per genotype, so cells and nuclei
are pseudoreplicates and the P values describe how reliably the two libraries
differ. Effect sizes accompany every test, resampling checks accompany the
composition and maturation results, and `tables/Table_S1_composition_FDR.csv`
reports the variance-inflation factor at which each composition result would
stop passing correction.

Multiple-testing correction is applied within each analysis family: across the
twelve cell classes for composition, the nine classes for maturation, all 90
combinations for the gene-programme battery, and within each enrichment run for
gene ontology.

## Citation

Please cite the original data publication:

Zhang, Y. et al. Disruption of CTCF boundary-dependent gene regulation by a
neurodevelopmental disorder-associated variant. *Nat. Commun.* **15**, 5524 (2024).
