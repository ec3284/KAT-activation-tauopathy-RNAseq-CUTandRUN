#!/bin/bash
# ==============================================================================
# Script 02: Consensus peak union across 4 conditions + heatmap
# Figures: 4C, S7B
#
# Requires: the four consensus peak BED files deposited in GEO
# (nTg-Veh-H3K4ac-consensus.peaks.bed, nTg-RA013915-H3K4ac-consensus.peaks.bed,
# hTau-Veh-H3K4ac-consensus.peaks.bed, hTau-RA013915-H3K4ac-consensus.peaks.bed),
# plus the per-condition averaged bigWig files from Script 01.
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-./cutrun_work}"
cd "$WORKDIR"

# ---- 3. Build 4-way union of consensus peaks (Fig. 4C, S7B) ----
cat nTg-Veh-H3K4ac-consensus.peaks.bed \
    nTg-RA013915-H3K4ac-consensus.peaks.bed \
    hTau-Veh-H3K4ac-consensus.peaks.bed \
    hTau-RA013915-H3K4ac-consensus.peaks.bed \
    > union_relaxed4_cat.bed
sortBed -i union_relaxed4_cat.bed > union_relaxed4_sorted.bed
mergeBed -i union_relaxed4_sorted.bed > union_relaxed4.bed
wc -l union_relaxed4.bed

# ---- 4. computeMatrix over the union of consensus peaks ----
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
    -R union_relaxed4.bed \
    -S nTg-Veh-H3K4ac-avg.bw nTg-RA013915-H3K4ac-avg.bw \
       hTau-Veh-H3K4ac-avg.bw hTau-RA013915-H3K4ac-avg.bw \
    --missingDataAsZero \
    -o matrix_union_relaxed4.gz -p max

# ---- 5. Heatmap + metaplot over consensus peak union (Fig. 4C) ----
plotHeatmap -m matrix_union_relaxed4.gz \
    -o heatmap_peaks_relaxed4.pdf \
    --colorMap inferno \
    --samplesLabel "nTg Veh" "nTg 915" "hTau Veh" "hTau 915" \
    --sortRegions descend --sortUsingSample 1 \
    --regionsLabel "Consensus Peaks" \
    --heatmapHeight 15 --interpolationMethod bilinear --dpi 300
