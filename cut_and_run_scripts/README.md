# CUT&RUN Analysis Scripts

Downstream analysis scripts for H3K4ac CUT&RUN data (nf-core/cutandrun v3.2.2
output as starting point). Cleaned up from an interactive lab-notebook
session transcript into standalone, ordered scripts, one per figure/panel.

## Contents

| Script | Figure(s) | Description |
|---|---|---|
| `01_bigwig_average_and_TSS_metaplot.sh` | 4A, 4B, S7A | Average bigWig across replicates; TSS metaplot |
| `02_consensus_peaks_union_heatmap.sh` | 4C, S7B | 4-way consensus peak union; heatmap |
| `03_feature_distribution_peaks_Fig4D.R` | 4D | Genomic feature distribution (promoter/exon/intron/intergenic) |
| `04_RNAseq_integration_Fig4E.sh` + `04b_calc_auc_rnaseq_relaxed.py` | 4E | H3K4ac signal (AUC) at RNA-seq DEG-up loci; Wilcoxon test |
| `05_overlap_and_rescued_peaks_Fig4F_G.sh` | 4F, 4G | Pairwise peak overlaps; identification of "rescued" peaks |
| `06_GO_rescued_peaks_dotplot_Fig4G.R` | 4G | GO Biological Process enrichment + dotplot of rescued peaks |
| `07_tracks_hTau_Fos_Rhob_Fig4H_I.sh` | 4H, 4I | Genome browser tracks, hTau samples |
| `08_venn_nTg_Fig_S7C.py` | S7C | nTg Veh vs nTg RA013915 peak Venn diagram |
| `09_tracks_nTg_Fos_Rhob_Fig_S7D_E.sh` | S7D, S7E | Genome browser tracks, nTg samples |
| `tracks_hTau_Fos.ini`, `tracks_hTau_Rhob.ini`, `tracks_nTg_Fos.ini`, `tracks_nTg_Rhob.ini` | 4H, 4I, S7D, S7E | pyGenomeTracks configuration files for genome browser tracks |

Run in numerical order; each script assumes the outputs of the previous ones
are present in the working directory (`$WORKDIR`, default `./cutrun_work`).

## Notes on data provenance

- **Fig. S7C Venn diagram counts** (`python/08_venn_nTg_Fig_S7C.py`): the
  peak-overlap counts (122,922 / 10,272 / 13,636) are generated with
  `intersectBed` against the consensus peak BED files; see the commands
  in the comment at the top of the script.
- **pyGenomeTracks `.ini` configs**: all four files referenced by scripts 07
  and 09 (`tracks_hTau_Fos.ini`, `tracks_hTau_Rhob.ini`, `tracks_nTg_Fos.ini`,
  `tracks_nTg_Rhob.ini`) are included in `ini_configs/`, with consistent
  track structure, colors, and signal scale (matched to the corresponding
  hTau version for each locus) across nTg and hTau variants.

## Software versions

- nf-core/cutandrun v3.2.2 (Nextflow v23.10.1)
- deepTools (bigwigAverage, computeMatrix, plotHeatmap, plotProfile)
- BEDTools (intersectBed, sortBed, mergeBed, closestBed)
- pyGenomeTracks
- R: ChIPseeker, TxDb.Mmusculus.UCSC.mm10.knownGene, org.Mm.eg.db,
  clusterProfiler, ggplot2
- Python: pandas, scipy, matplotlib
