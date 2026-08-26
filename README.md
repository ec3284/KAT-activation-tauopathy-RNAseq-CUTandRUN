[README main.md](https://github.com/user-attachments/files/31489901/README.main.md)
# KAT Activation Restores Chromatin Dynamics and Resolves Tau Pathology in Tauopathy Mice — Analysis Scripts

Analysis scripts supporting the RNA-seq and CUT\&RUN data reported in this
study, characterizing the hippocampal transcriptional and chromatin
response to the KAT activators RA013915 and YF2 in non-transgenic (nTg)
and hTau/Mapt-KO (hTau) mice.

## Repository structure

```
├── RNAseq/       Differential expression analysis (DESeq2) and downstream
│                 gene set definitions used for RNA-seq figures and for
│                 integration with CUT\&RUN data.
└── CUTandRUN/    H3K4ac CUT\&RUN downstream analysis scripts (peak calling
                  integration, genomic feature distribution, GO enrichment,
                  genome browser tracks, RNA-seq integration).
```

Each subfolder contains its own `README.md` with a script-by-script index
mapped to the corresponding manuscript figure panels, software versions,
and run instructions.

## Raw and processed data

Raw sequencing files and processed data (count matrices, bigWig tracks,
peak calls) are deposited in the Gene Expression Omnibus (GEO):

* RNA-seq: [GSE345197](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345197)
* CUT\&RUN: [GSE343636](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343636)

## Citation

If you use these scripts, please cite the associated manuscript:

*\[Manuscript citation to be added upon publication]*

