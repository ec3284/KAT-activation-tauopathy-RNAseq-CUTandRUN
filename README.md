# KAT Activation Restores Chromatin Dynamics and Resolves Tau Pathology in Tauopathy Mice — Analysis Scripts

Analysis scripts supporting the RNA-seq and CUT&RUN data reported in this
study, characterizing the hippocampal transcriptional and chromatin response
to the KAT activators RA013915 and YF2 in non-transgenic (nTg) and
hTau/Mapt-KO (hTau) mice.

## Repository structure

```
├── RNAseq/        Differential expression (limma-voom, with the sequencing
│                  run in the design and the animal as a blocking factor),
│                  gene-set enrichment, a comparison against four alternative
│                  model specifications, and the scripts for each panel of
│                  Figure 3 and its supplementary figures.
└── CUTandRUN/     H3K4ac CUT&RUN downstream analysis: spike-in-normalized
                   signal, consensus peak handling, differential analysis over
                   the peak union, genomic feature distribution, GO
                   enrichment, and genome browser tracks.
```

Each subfolder has its own `README` with a script-by-script index mapped to
the manuscript figure panels, the run order, the region sets, and software
versions.

## Running the scripts

Scripts are numbered in the order they must be run, and none of them contains
an absolute path. Run each from the directory where its outputs should be
written.

In `RNAseq/`, `01` builds the model and writes the results table every later
script reads; the inputs are expected in the working directory.

In `CUTandRUN/`, the input locations are given by environment variables
(`CUTRUN_PEAKS`, `CUTRUN_BIGWIG`, `CUTRUN_BAM`, `GTF` and a few others), each
with a sensible default; that folder's README lists them. Scripts `01`–`04`
build the files every panel reads, and the later scripts produce one panel
each.

Two panels of each figure set are drawn in GraphPad Prism from exported CSV
files rather than by the scripts. That is noted where it applies and is
expected, not a missing step.

## Raw and processed data

Raw sequencing files and processed data (count matrices, bigWig tracks, peak
calls) are deposited in the Gene Expression Omnibus (GEO):

* RNA-seq: [GSE345197](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345197)
* CUT&RUN: [GSE343636](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343636)

The RNA-seq deposition includes a sample sheet mapping each library to its
biological sample, which identifies the four hTau Vehicle animals that were
sequenced in both runs.

## A note on the two designs

The two datasets are analysed with different models for a reason worth stating
once.

The RNA-seq design is partly confounded: nTg Vehicle and nTg RA013915 were
sequenced only in the first run and hTau RA013915 only in the second, and the
four hTau Vehicle animals were sequenced in both. Those repeated libraries are
what makes the run effect estimable. Because each of the five run × group
cells holds four libraries, the design is saturated in those cells, so the
contrasts are algebraically identical to the corresponding within-run
differences and the blocking factor changes the standard errors rather than
the fold changes. `RNAseq/11_model_comparison.R` fits four alternative
specifications and shows this directly.

The CUT&RUN design has three mice per group in a single batch, so the question
there is not the model but the unit of replication: every comparison is made
with the mouse as the experimental unit, over a common set of regions, with
spike-in-derived size factors.

## Citation

If you use these scripts, please cite the associated manuscript:

*[Manuscript citation to be added upon publication]*
