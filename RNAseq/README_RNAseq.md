# RNA-seq Analysis Scripts

Downstream analysis of hippocampal RNA-seq from non-transgenic (nTg) and
hTau/Mapt-KO (hTau) mice treated with the KAT activator RA013915. Twenty
libraries from sixteen animals, four per group, sequenced in two runs.

Scripts are standalone and numbered in the order they must be run: `01`
builds the model and writes the results table that every later script
reads, `02` computes the gene-set enrichment that two of the figure
scripts depend on, and `03`–`11` each produce one panel or one check.

## Contents

| Script | Figure | Description |
|-|-|-|
| `01_differential_expression_model.R` | source data | Builds the count matrix from the two provider runs, fits limma-voom with sequencing run in the design and the animal as a blocking factor, extracts four contrasts, writes the results table and the files for GEO deposition |
| `02_fgsea_both_contrasts.R` | source data | Gene-set enrichment on the disease and treatment contrasts, genes ranked by the moderated *t* statistic |
| `03_volcano_disease_Fig3A.R` | 3a | Genes altered by the tauopathy |
| `04_gsea_disease_processes_Fig3B.R` | 3b | Biological processes altered by the tauopathy |
| `05_scatter_disease_vs_treatment_Fig3C.R` | 3c | Transcriptome-wide relationship between the disease effect and the treatment effect |
| `06_heatmap_rescued_genes_Fig3D.R` | 3d | Expression of the rescued genes across the four groups |
| `07_scatter_treatment_specificity_Fig3E.R` | 3e | The treatment effect in nTg mice against the same effect in hTau mice |
| `08_gsea_mirror_Fig3F.R` | 3f | Gene sets that the pathology and the treatment move in opposite directions |
| `09_MAplot_nTg_FigED5A.R` | Extended Data 5a | RA013915 versus vehicle in non-transgenic mice |
| `10_qc_batch_and_technical_replicates_FigS1.R` | Supplementary 1 | Sequencing-run effect and its removal |
| `11_model_comparison.R` | source data | Differentially expressed gene counts under four alternative model specifications |

Run `01` first. `04` and `08` additionally require `02`. Every other
script reads only the outputs of `01`, so `03`, `05`, `06`, `07`, `09`,
`10` and `11` can be run in any order once `01` has completed.

All paths are relative to the working directory. Set it to the folder
holding the two provider matrices (`C3SVT8_results/` and
`VJ385X_results/`) before running `01`; every later script reads and
writes in that same folder.

## Notes on data provenance

* **Two sequencing runs, partly confounded with group.** nTg Vehicle and
  nTg RA013915 were sequenced only in run C3SVT8 and hTau RA013915 only in
  run VJ385X. The four hTau Vehicle animals were sequenced in *both* runs,
  and those repeated libraries are what make the run effect estimable.
  `01` verifies that the design matrix (`~ sequencing run + group`) is of
  full rank before fitting anything, and stops if it is not.

* **Repeated libraries are not collapsed.** Twenty libraries come from
  sixteen animals. Rather than summing the repeated hTau Vehicle
  libraries, which would destroy the information that makes the run
  effect estimable and would double that group's depth alone, the animal
  is declared as a blocking factor through `duplicateCorrelation`
  (consensus within-animal correlation 0.62). The hTau Vehicle group
  therefore contributes four animals, not eight libraries, to the
  residual degrees of freedom.

* **Rescue gene definition.** A gene must be significant (adjusted
  *P* < 0.05) in *both* the disease contrast (hTau Vehicle versus nTg
  Vehicle) and the treatment contrast (hTau RA013915 versus hTau Vehicle),
  with effects of opposite sign. This is a direction-aware cross-filter,
  not the intersection of two independently thresholded lists.

* **Genotype selectivity is tested formally.** Comparing the number of
  significant genes in two separate contrasts is not a test, because the
  difference between "significant" and "not significant" is not itself
  significant. `01` therefore fits a genotype × treatment interaction
  genome-wide and reports it as a fourth contrast.

* **Model choice.** `11` compares the model used here against four
  alternative specifications of the same data (pseudo-replicate and
  within-run designs, under both DESeq2 and limma). The specification
  used in the manuscript is among the most conservative of the five in
  both the disease and the treatment contrast.

## Software versions

* fastp v0.24.0 (Chen et al., 2018; RRID:SCR_016962), including UMI
  deduplication with UMICollapse v1.1.0
* STAR aligner v2.7.x (Dobin et al., 2013)
* featureCounts, Subread package v2.1.1 (Liao et al., 2014)
* R v4.5.1: limma, edgeR, fgsea, DESeq2 (variance-stabilizing transform
  for the heatmap only), pheatmap, ggplot2, ggrepel, dplyr, readr,
  openxlsx, writexl, patchwork
* Gene sets: MSigDB collection C5 (GO Biological Process)
