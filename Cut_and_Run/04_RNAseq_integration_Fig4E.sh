#!/bin/bash
# ==============================================================================
# Script 04: Integration of CUT&RUN H3K4ac signal with RNA-seq DEG-up genes
# Figure: 4E
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-./cutrun_work}"
cd "$WORKDIR"

# ---- Intersect genome-wide peak union with RNA-seq DEG-up gene regions ----
# DEG_up_regions.bed: promoter/gene regions of genes upregulated by RA013915
# in the RNA-seq differential expression analysis (see RNAseq_analysis/ scripts)
intersectBed -a union_relaxed4.bed -b DEG_up_regions.bed -u \
    > union_relaxed_on_DEGup.bed
wc -l union_relaxed_on_DEGup.bed

# ---- computeMatrix over DEG-up peaks, per-replicate bigWig tracks ----
computeMatrix reference-point --referencePoint center -b 2000 -a 2000 \
    -R union_relaxed_on_DEGup.bed \
    -S nTg-Veh-H3K4ac-rep1.bw nTg-Veh-H3K4ac-rep2.bw nTg-Veh-H3K4ac-rep3.bw \
       hTau-Veh-H3K4ac-rep1.bw hTau-Veh-H3K4ac-rep2.bw hTau-Veh-H3K4ac-rep3.bw \
       hTau-RA013915-H3K4ac-rep1.bw hTau-RA013915-H3K4ac-rep2.bw hTau-RA013915-H3K4ac-rep3.bw \
    --missingDataAsZero \
    -o matrix_relaxed_DEGup_3cond.gz -p max

# ---- AUC quantification + Wilcoxon signed-rank test ----
python3 04b_calc_auc_rnaseq_relaxed.py

# ---- Figure 4E ----
# Data source: AUC_DEGup_NEW760_3cond.csv (output of 04b above), using the
# 521 real H3K4ac consensus peaks overlapping the promoter regions
# (TSS +/-2kb) of the 505 mappable RA013915-upregulated genes (padj<0.05,
# log2FC>0; DESeq2, collapsed model). Statistics: paired Wilcoxon
# signed-rank test (see WILCOXON_DEGup_NEW760.csv).
# The violin plot itself was made manually in GraphPad Prism from
# AUC_DEGup_NEW760_3cond.csv, not scripted -- this is expected, not a
# missing step.
