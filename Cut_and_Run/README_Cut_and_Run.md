# CUT&RUN Analysis Scripts

Downstream analysis scripts for H3K4ac CUT&RUN data (nf-core/cutandrun v3.2.2
output as starting point). Cleaned up from an interactive lab-notebook
session transcript into standalone, ordered scripts, one per figure/panel.

## Contents

| Script | Figure(s) | Description |
|---|---|---|
| `bash/01_bigwig_average_and_TSS_metaplot.sh` | 4A, 4B, S7A | Average bigWig across replicates; TSS metaplot; genome-wide promoter H3K4ac AUC (all 12 replicates) + paired Wilcoxon tests (inline Python) |
| `bash/02_consensus_peaks_union_heatmap.sh` | 4C, S7B | 4-way consensus peak union; heatmap |
| `R/03_feature_distribution_peaks_Fig4D.R` | 4D | Genomic feature distribution (promoter/exon/intron/intergenic) |
| `bash/04_RNAseq_integration_Fig4E.sh` + `python/04b_calc_auc_rnaseq_relaxed.py` | 4E | H3K4ac signal (AUC) at RNA-seq DEG-up loci; Wilcoxon test |
| `bash/05_overlap_and_rescued_peaks_Fig4F_G.sh` | 4F, 4G | Pairwise peak overlaps; identification of "rescued" peaks |
| `R/06_GO_rescued_peaks_dotplot_Fig4G.R` | 4G | GO Biological Process enrichment + dotplot of rescued peaks |
| `bash/07_tracks_hTau_Fos_Dusp1_Fig4H_I.sh` | 4H, 4I | Genome browser tracks (Fos, Dusp1), hTau samples |
| `python/08_venn_nTg_Fig_S7C.py` | S7C | nTg Veh vs nTg RA013915 peak Venn diagram |
| `bash/09_tracks_nTg_Fos_Dusp1_Fig_S7D_E.sh` | S7D, S7E | Genome browser tracks (Fos, Dusp1), nTg samples |
| `ini_configs/` | 4H, 4I, S7D, S7E | pyGenomeTracks configuration files for genome browser tracks |

Run in numerical order; each script assumes the outputs of the previous ones
are present in the working directory (`$WORKDIR`, default `./cutrun_work`).

See `../RNAseq/` for the RNA-seq analysis scripts (differential expression,
GO/TF enrichment of rescue genes).

## Notes on data provenance

- **Fig. S7C Venn diagram counts** (`python/08_venn_nTg_Fig_S7C.py`): the
  peak-overlap counts (122,922 / 10,272 / 13,636) are generated with
  `intersectBed` against the consensus peak BED files; see the commands
  in the comment at the top of the script.
- **pyGenomeTracks `.ini` configs**: all four files referenced by scripts 07
  and 09 (`tracks_hTau_Fos.ini`, `tracks_hTau_Dusp1.ini`, `tracks_nTg_Fos.ini`,
  `tracks_nTg_Dusp1.ini`) are included in `ini_configs/`, with consistent
  track structure, colors, and signal scale (matched to the corresponding
  hTau version for each locus) across nTg and hTau variants.
- **Figure 4B/S7A (genome-wide promoter H3K4ac AUC)**: computed inline
  (Python) at the end of `bash/01_bigwig_average_and_TSS_metaplot.sh`. Step
  3 runs a 12-replicate (4 conditions x 3 replicates) computeMatrix at TSS
  +/-3kb, `-R` set to the same full GTF used in step 2 (`$GTF`); step 4
  builds `transcript_to_gene.tsv` from that same GTF to collapse multiple
  transcripts per gene; step 5 computes the AUC and runs the paired
  Wilcoxon tests, exporting `AUC_PROMOTER_PC.PRISM_4conditions.csv` and
  `AUC_PROMOTER_PC.WILCOXON_results.csv`.
- **Figure 4E (H3K4ac vs RA013915-upregulated genes)**: this panel is a
  violin plot of promoter-proximal H3K4ac AUC at peaks overlapping
  RA013915-upregulated gene promoters, compared across all three
  conditions by paired Wilcoxon signed-rank test. The AUC/statistics step
  (`python/04b_calc_auc_rnaseq_relaxed.py`) is self-contained in this
  folder, exporting `AUC_DEGup_NEW760_3cond.csv` and
  `WILCOXON_DEGup_NEW760.csv`.
- **All violin plots (Fig. 4B, S7A, 4E) were made manually in GraphPad
  Prism from the exported CSVs above, not scripted** -- this is expected,
  not a missing step.


## Software versions

- nf-core/cutandrun v3.2.2 (Nextflow v23.10.1)
- deepTools (bigwigAverage, computeMatrix, plotHeatmap, plotProfile)
- BEDTools (intersectBed, sortBed, mergeBed, closestBed)
- pyGenomeTracks
- R: ChIPseeker, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db,
  clusterProfiler, ggplot2
- Python: pandas, scipy, matplotlib
- GraphPad Prism (violin plots: Fig. 4B, S7A, 4E, from the exported CSVs
  above)
