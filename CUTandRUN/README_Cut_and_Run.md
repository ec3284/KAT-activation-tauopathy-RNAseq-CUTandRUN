# CUT&RUN Analysis Scripts

Downstream analysis of the H3K4ac CUT&RUN data, starting from the
nf-core/cutandrun v3.2.2 output. One script per panel, run in numerical order:
`01`–`04` build the files every panel reads, `05`–`10` are Figure 4, and
`11`–`14` are Extended Data Fig. 5.

## Contents

| Script | Panel | Description |
|---|---|---|
| `01_bigwig_average_and_matrices.sh` | — | Per-condition averaged bigWigs; TSS computeMatrix over all twelve libraries |
| `02_consensus_peak_union_and_heatmap_Fig4C_FigED5D.sh` | 4c, ED 5d | Three- and four-condition peak unions; a heatmap over each |
| `03_union_peak_counts.sh` | — | Raw read counts per mouse over the peak union; spike-in totals |
| `04_differential_peaks_spikein.R` | — | DESeq2 over the union with spike-in size factors; leave-one-out; the per-contrast tables |
| `05_TSS_metaplot_with_IgG_Fig4A.R` | 4a | Spike-in-normalized H3K4ac profile at the TSS, with the IgG control |
| `06_promoter_AUC_per_mouse_Fig4B_FigED5C.R` | 4b, ED 5c | Promoter AUC per animal, for both the hTau and the nTg contrast |
| `07_feature_distribution_responsive_Fig4D.R` | 4d | Where the RA013915-responding regions sit, against a coverage-matched background |
| `08_overlap_and_rescued_peaks_Fig4E.sh` | 4e | Pairwise overlap counts, both directions; the rescued peak set |
| `09_GO_rescued_peaks_Fig4F.R` | 4f | GO Biological Process enrichment of the rescued loci, against the marked genome |
| `10_tracks_hTau_Fos_Dusp1_Fig4G_H.sh` | 4g, 4h | Genome browser tracks at Fos and Dusp1, hTau conditions |
| `11_feature_distribution_four_conditions_FigED5B.R` | ED 5b | Genomic feature distribution of the four consensus sets, with a random control |
| `12_promoter_fraction_per_replicate_FigED5B.R` | ED 5b | The per-mouse promoter fraction behind that panel |
| `13_venn_nTg_FigED5E.py` | ED 5e | nTg vehicle against nTg RA013915 peak overlap |
| `14_tracks_nTg_Fos_Dusp1_FigED5F_G.sh` | ED 5f, 5g | Genome browser tracks at Fos and Dusp1, nTg conditions |
| `tracks_*.ini` | 4g, 4h, ED 5f, 5g | pyGenomeTracks configurations |

Seven scripts are R (`04`–`07`, `09`, `11`, `12`), six are shell (`01`–`03`,
`08`, `10`, `14`) and one is Python (`13`). The shell scripts regenerate the
intermediate files; with the deepTools matrices, the unions and the count
table already in place, the R scripts run on their own.

The two violin panels (Fig. 4b and Extended Data Fig. 5c) are drawn in
GraphPad Prism from the `*_PRISM_*.csv` files `06` exports. That is expected,
not a missing step.

See `../RNAseq/` for the transcriptomic analysis.

## Paths

No script contains an absolute path. Each is run from the directory where its
outputs should be written, and reads its inputs from locations given by
environment variables, all with defaults:

| Variable | Default | What it holds |
|---|---|---|
| `CUTRUN_PEAKS` | `peaks_and_consensus_peaks` | consensus and per-replicate SEACR peak BED files |
| `CUTRUN_REPLICATE_PEAKS` | `$CUTRUN_PEAKS` | the per-library peak files, if kept separately |
| `CUTRUN_BIGWIG` | `bigwig` | spike-in normalized bigWig tracks |
| `CUTRUN_BAM` | `bam` | the twelve coordinate-sorted H3K4ac BAM files |
| `CUTRUN_SPIKE_BAM` | `$CUTRUN_BAM` | the spike-in BAM files, for `03` |
| `CUTRUN_GENES_BED` | `mm10_genes_sorted.bed` | mm10 gene annotation as BED, for `09` |
| `GTF` | `Mus_musculus.GRCm38.102.gtf.gz` | the annotation used by deepTools and pyGenomeTracks |

`mm10_genes_sorted.bed` is built from the same GTF:

```bash
zcat Mus_musculus.GRCm38.102.gtf.gz \
  | awk -v OFS='\t' '$3=="gene"' \
  | sed -n 's/^\([^\t]*\)\t[^\t]*\tgene\t\([0-9]*\)\t\([0-9]*\).*gene_name "\([^"]*\)".*/chr\1\t\2\t\3\t\4/p' \
  | sortBed -i - > mm10_genes_sorted.bed
```

Two annotations are in play: the deepTools and pyGenomeTracks steps use the
Ensembl GRCm38 release 102 GTF, while the R scripts annotate against
`TxDb.Mmusculus.UCSC.mm10.knownGene`. Fig. 4a and Fig. 4b are therefore both
"at the TSS" but over different transcript sets.

## Region sets

| Set | Regions | Used by |
|---|---|---|
| Union of three conditions | 295,426 | Fig. 4c heatmap and its metaplot |
| Union of four conditions | 299,069 | Extended Data Fig. 5d heatmap |
| Same, chrM excluded | 299,068 | every quantitative comparison (`03`, `04`) |
| After the ≥120-read filter | 101,161 | the DESeq2 tests |

The one-region difference is a single peak called on chrM. `02` prints all of
these when it runs, and `11` rebuilds the four-condition union in R with
`GenomicRanges::reduce`, which is the same operation as `mergeBed`; the two
agree.

## Notes on the analysis

**The peak counts are not a result.** The four consensus sets hold 133,194 /
23,458 / 83,121 / 191,623 regions. Sequencing depth and alignment are
comparable across the twelve libraries (mean 14.3 million read pairs, range
11.8–17.1 million), and the group calling the fewest peaks — RA013915-treated
nTg, 23,458 — has the lowest duplicate rate of the four and ordinary depth. Its
libraries are the cleanest of the set. The peak-number differences therefore
report how readily each library passed the SEACR relaxed threshold, not how
much H3K4ac it carries. Every quantitative comparison here works from read
counts over a common set of regions (`03` and `04`), and the panels built from
the consensus sets are read as composition rather than amount, with the region
count printed next to each bar.

**The experimental unit is the mouse.** `06` and `12` reduce each library to
one number. An earlier version of the promoter AUC analysis collapsed the
replicates to a group median before testing, which leaves only a test across
genes; that step is gone, and the removed code is quoted in the header of `01`
so the change is on the record rather than silent.

**No hypothesis test is reported for the promoter AUC panels.** A test across
genes is invalid. A test on three per-animal medians is admissible but
underpowered: the between-animal standard deviation within groups is 0.55 log₂,
so the smallest difference detectable at 80% power with n = 3 is 1.69 log₂
(3.2-fold) — twice the 1.8-fold effect observed. The panels give the fold
change between the means of the per-animal medians. `06` still fits those tests
as a sensitivity check and writes them to
`Fig4b_perAnimal_tests_sensitivity.csv`; those values are deliberately absent
from the figures, the legends and the Statistics section.

**No region reaches FDR < 0.05 in any contrast.** The panel that asks where the
responding regions sit (`07`) therefore defines its set at a nominal threshold,
says so throughout, and tests it against a background matched on mean
normalized count rather than against all peaks, because a region with more
reads has more power and promoter peaks carry more reads than distal ones. The
non-responsive set is the internal control and moves the opposite way.

**The spike-in normalization is load-bearing, so it is checked.** The twelve
spike-in counts hard-coded in `04` match the E. coli alignment reported by the
pipeline for all twelve libraries. `04` prints them next to the size factors
DESeq2 would have estimated itself, and repeats the hTau contrast without
`hTau_915_R1`, whose spike-in count is roughly half the next lowest. A result
that does not survive that leave-one-out is a property of one library's
spike-in count.

**The GO test has a universe.** `09` uses the genes carrying any H3K4ac peak as
the background (21,900 genes; 18,362 after mapping to Entrez identifiers), not
every gene in the annotation, so the question is whether the rescued loci
differ from the rest of the marked genome. It also removed a keyword filter
that selected synaptic terms by name before plotting. The rescued set is 43% of
that universe, which caps fold enrichment at 2.31; the script prints that
fraction before anything else, and it should be read before the dotplot.

**The two-circle overlap counts are asymmetric.** `intersectBed` counts the
peaks of one set that meet the other, and that count differs by direction
because peaks merge and split differently between two calls. A Venn has one
slot for the intersection and so cannot satisfy both circle totals at once:
the nTg intersection is 10,272 from the vehicle side and 9,822 from the
RA013915 side. `08` and `13` compute both directions, print the discrepancy,
and label each circle with the consensus total it came from.

**The metaplot is drawn on a subsample.** `05` profiles a random sample of
40,000 annotated transcript TSS with a fixed seed. The mean profile is
indistinguishable from the one over every TSS; set `SUBSAMPLE <- NA` to use all
of them, and change the legend if you do.

## Software versions

- nf-core/cutandrun v3.2.2 (Nextflow v23.10.1)
- deepTools v3.5.6: `bigwigAverage`, `computeMatrix`, `plotHeatmap`
- BEDTools v2.30.0: `intersectBed`, `sortBed`, `mergeBed`, `multicov`
- SAMtools
- SEACR v1.3 (relaxed mode)
- pyGenomeTracks
- R v4.5.1: DESeq2, ChIPseeker, TxDb.Mmusculus.UCSC.mm10.knownGene,
  org.Mm.eg.db, clusterProfiler, GenomicRanges, rtracklayer, data.table,
  readr, dplyr, ggplot2, writexl
- Python: matplotlib
- GraphPad Prism 10.5.0 (the two violin panels, from the exported
  `*_PRISM_*.csv` files)
