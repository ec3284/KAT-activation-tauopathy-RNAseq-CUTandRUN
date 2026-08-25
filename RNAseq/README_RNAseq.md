# RNA-seq Analysis Scripts

Downstream analysis scripts for hippocampal RNA-seq data (Plasmidsaurus
UMI-based single-end sequencing, two independent batches). Cleaned up from
an interactive lab-notebook session transcript into standalone, ordered
scripts, one per figure/panel.

## Contents

|Script|Figure(s)|Description|
|-|-|-|
|`01\_DESeq2\_collapsed\_rescue\_analysis.R`|3A (data), source data|Core pipeline: load both sequencing batches, collapse technical replicates (`collapseReplicates`), run DESeq2, define cross-filtered "rescue" genes|
|`02\_MAplots.R`|S6B, S6C|MA plots (nTg Veh vs 915; hTau Veh vs 915)|
|`03\_GO\_Fisher\_enrichment.R`|3B|GO enrichment of the disease signature; Fisher's exact test enrichment of disease-associated pathways among treatment-responsive genes|
|`04\_volcano\_plot.R`|3A|Categorized volcano plot (hTau: RA013915 vs Vehicle), points colored by pathway category, rescue genes labeled|
|`05\_heatmap\_rescue\_genes.R`|3C|Heatmap of the 22 rescue genes across all four experimental groups|
|`06\_GO\_2026\_barplot.R`|3D|GO Biological Process 2026 enrichment barplot of the 22 rescue genes (Enrichr)|
|`07\_TF\_enrichment\_barplot.R`|3E|Transcription factor enrichment barplot (ChEA, ARCHS4, TRRUST) of the 22 rescue genes|

Run `01` first; scripts `02`–`06` and `07` assume its output objects
(`df\_merged`, `rescued\_up`, `rescued\_down`, `dds`, `vsd`, `res\_nTg`,
`res\_treatment`) are in memory, or can be reloaded from
`RNAseq\_FINAL\_DESeq2\_results.xlsx`. Script `03` must be run before `04`
(supplies `mito\_all`, `neuronal\_all`, `synaptic\_all` for volcano coloring).
Notes on data provenance

* **Two sequencing batches**: nTg (Vehicle and RA013915) and one set of
hTau Vehicle replicates were sequenced together (batch C3SVT8); a second
set of hTau Vehicle replicates (the same biological RNA, re-sequenced)
and all hTau RA013915 replicates were sequenced in a separate batch
(batch VJ385X). Because a large, unexplained batch effect was detected
when modeling batch as a covariate (\~3,900 genes, likely reflecting
technical differences between sequencing runs), the two hTau Vehicle
technical replicate sets were instead collapsed (summed raw counts) using
`DESeq2::collapseReplicates()`, avoiding the need to statistically
estimate the batch effect.
* **Rescue gene definition**: genes must be significant (padj<0.05) in
*both* the disease comparison (hTau Vehicle vs nTg Vehicle) and the
treatment comparison (hTau RA013915 vs hTau Vehicle), with effects in
opposite directions. This is a direction-aware cross-filter, distinct
from simply intersecting two independently thresholded gene lists.
* **GO/TF enrichment values (scripts 06, 07)**: these scripts do not
re-run the enrichment; they plot values recorded from Enrichr
(maayanlab.cloud/Enrichr) web submissions of the 22 rescue gene symbols.
To regenerate, resubmit the gene list to Enrichr and export the relevant
library results (see comments at the top of each script).

## Software versions

* FastP v0.24.0 (Chen et al., 2018; RRID:SCR\_016962), including UMI
deduplication (UMICollapse v1.1.0)
* STAR aligner v2.7.x
* featureCounts (Subread package v2.1.1)
* R: DESeq2, clusterProfiler, org.Mm.eg.db, ggplot2, ggrepel, pheatmap,
dplyr, tidyr, openxlsx
* Enrichr (maayanlab.cloud/Enrichr): GO Biological Process 2026, ChEA
ChIP-seq, ARCHS4 TFs Coexp, TRRUST Transcription Factors 2019

