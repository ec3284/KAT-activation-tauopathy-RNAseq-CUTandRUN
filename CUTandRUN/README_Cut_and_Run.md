# CUT&RUN Analysis Scripts

Downstream analysis of the H3K4ac CUT&RUN data, starting from the
nf-core/cutandrun v3.2.2 output. One script per panel, run in numerical order:
`01`–`04` build the files every panel reads, `05`–`09` are Figure 4, and
`10`–`13` are Extended Data Fig. 5.

## Contents

| Script | Panel | Description |
|---|---|---|
| `01_bigwig_average.sh` | — | Per-condition averaged bigWigs, three replicates per condition |
| `02_consensus_peak_union_and_heatmap_Fig4B_FigED5B.sh` | 4b, ED 5b | Three- and four-condition peak unions; a heatmap over each |
| `03_union_peak_counts.R` | — | Raw read counts per mouse over the peak union, one per fragment; spike-in totals |
| `04_differential_peaks_spikein.R` | — | DESeq2 over the union with spike-in size factors; leave-one-out; the per-contrast tables |
| `05_TSS_metaplot_with_IgG_Fig4A.R` | 4a | Spike-in-normalized H3K4ac profile at the TSS, with the IgG control |
| `06_feature_distribution_responsive_Fig4C.R` | 4c | Where the RA013915-responding regions sit, against a coverage-matched background |
| `07_overlap_and_rescued_peaks_Fig4D.sh` | 4d | Pairwise overlap counts, both directions; the rescued peak set |
| `08_GO_rescued_peaks_Fig4E.R` | 4e | GO Biological Process enrichment of the rescued loci, against the marked genome |
| `09_tracks_hTau_Fos_Dusp1_Fig4F_G.sh` | 4f, 4g | Genome browser tracks at Fos and Dusp1, hTau conditions |
| `10_feature_distribution_four_conditions_FigED5D.R` | ED 5d | Genomic feature distribution of the four consensus sets, with a random control |
| `11_promoter_fraction_per_replicate_FigED5D.R` | ED 5d | The per-mouse promoter fraction behind that panel |
| `12_venn_nTg_FigED5C.py` | ED 5c | nTg vehicle against nTg RA013915 peak overlap |
| `13_tracks_nTg_Fos_Dusp1_FigED5E_F.sh` | ED 5e, 5f | Genome browser tracks at Fos and Dusp1, nTg conditions |
| `tracks_*.ini` | 4f, 4g, ED 5e, 5f | pyGenomeTracks configurations |

Seven scripts are R (`03`–`06`, `08`, `10`, `11`), five are shell (`01`, `02`,
`07`, `09`, `13`) and one is Python (`12`). The shell scripts regenerate the
intermediate files; with the deepTools matrices, the unions and the count table
already in place, the R scripts run on their own.

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
| `CUTRUN_GENES_BED` | `mm10_genes_sorted.bed` | mm10 gene annotation as BED, for `08` |
| `GTF` | `Mus_musculus.GRCm38.102.gtf.gz` | the annotation used by pyGenomeTracks, and the source of `mm10_genes_sorted.bed` |

`mm10_genes_sorted.bed` is built from the same GTF:

```bash
zcat Mus_musculus.GRCm38.102.gtf.gz \
  | awk -v OFS='\t' '$3=="gene"' \
  | sed -n 's/^\([^\t]*\)\t[^\t]*\tgene\t\([0-9]*\)\t\([0-9]*\).*gene_name "\([^"]*\)".*/chr\1\t\2\t\3\t\4/p' \
  | sortBed -i - > mm10_genes_sorted.bed
```

Two annotations are in play: the pyGenomeTracks steps use the Ensembl GRCm38
release 102 GTF, while the R scripts annotate against
`TxDb.Mmusculus.UCSC.mm10.knownGene`. The TSS metaplot (Fig. 4a) and the
feature annotation of the peak sets therefore rest on different transcript
sets.

## Region sets

| Set | Regions | Used by |
|---|---|---|
| Union of three conditions | 295,426 | Fig. 4b heatmap and its metaplot |
| Union of four conditions | 299,069 | Extended Data Fig. 5b heatmap |
| Counted per library | 299,068 | every quantitative comparison (`03`, `04`) |
| After the ≥120-read filter | 101,161 | the DESeq2 tests |

Counting covers the first 299,068 regions of the union; the last one,
`chrX:169,929,542-169,931,142`, is absent from the count matrix and from
everything built on it. `02` prints all of these counts when it runs, and `10`
rebuilds the four-condition union in R with `GenomicRanges::reduce`, which is
the same operation as `mergeBed`; the two agree.

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

**The experimental unit is the mouse.** `11` reduces each library to one
number, the fraction of that animal's peaks that fall in a promoter. An earlier
version of the promoter analysis collapsed the replicates to a group median
before testing, which leaves only a test across genes; that step is gone, and
the removed code is quoted in the header of `11` so the change is on the record
rather than silent. No hypothesis test is reported for the promoter fraction:
with three mice per group the between-animal spread leaves the comparison
underpowered, and the panel gives the composition itself.

**No region reaches FDR < 0.05 in any contrast.** The panel that asks where the
responding regions sit (`06`) therefore defines its set at a nominal threshold,
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

**The counts were checked against the alignments.** The matrix was rebuilt from
the twelve BAM files over the same regions, counting each read pair once as a
fragment and excluding duplicate-flagged alignments. The two agree on 94.8% to
97.8% of all 299,068 cells per library, with a largest single difference of 7
reads and correlations of 0.9986 to 0.9996; library totals differ by 0.29% to
0.69%. The residual is a fragment-boundary convention, uniform along the
genome, and no region crosses the >= 120-read filter because of it.

**The GO test has a universe.** `08` uses the genes carrying any H3K4ac peak as
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
RA013915 side. `07` and `12` compute both directions, print the discrepancy,
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
