#!/bin/bash
# ==============================================================================
# Union of the consensus peak sets, and the heatmaps over it
# Figure 4c (three conditions) and Extended Data Fig. 5d (all four)
#
# The union is the common set of regions every quantitative analysis in this
# folder is built on, which is the point of building it here rather than inside
# each script: the differential test, the feature distributions and the
# heatmaps all refer to the same regions.
#
# A note that belongs next to the peak numbers rather than buried in a legend.
# The four consensus sets hold 133,194 / 23,458 / 83,121 / 191,623 regions at
# comparable sequencing depth and comparable fraction of reads in peaks. A set
# eight times larger than another is reporting how readily that library passed
# the SEACR threshold, not how much H3K4ac it carries, so the peak counts are
# not a result and the panels below are read as composition rather than amount.
#
# Requires: bedtools (sortBed, mergeBed), deepTools (computeMatrix, plotHeatmap)
#
# Input : the four consensus peak BED files, and the averaged bigWigs from 01_
# Output: union_relaxed3.bed   the three conditions of Fig. 4c (295,426 regions)
#         union_relaxed4.bed   all four conditions, the set used downstream
#         Fig4c_heatmap_union_3cond.pdf
#         FigED5d_heatmap_union_4cond.pdf
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-.}"
PEAKS="${CUTRUN_PEAKS:-peaks_and_consensus_peaks}"
THREADS="${THREADS:-max}"
cd "$WORKDIR"
[ -d "$PEAKS" ] || { echo "peak directory not found: $PEAKS  (set CUTRUN_PEAKS)" >&2; exit 1; }

NTG_VEH="$PEAKS/nTg-Veh-H3K4ac-consensus.peaks.bed"
NTG_915="$PEAKS/nTg-RA013915-H3K4ac-consensus.peaks.bed"
HTAU_VEH="$PEAKS/hTau-Veh-H3K4ac-consensus.peaks.bed"
HTAU_915="$PEAKS/hTau-RA013915-H3K4ac-consensus.peaks.bed"
for f in "$NTG_VEH" "$NTG_915" "$HTAU_VEH" "$HTAU_915"; do
    [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }
    printf "%-52s %9d regions\n" "$(basename "$f")" "$(wc -l < "$f")"
done

# ---- 1. the two unions -------------------------------------------------------
# Two region sets, because the two panels use two. Fig. 4c is drawn over the
# union of the THREE conditions it displays; Extended Data Fig. 5d adds the
# fourth condition and is drawn over the union of all four, which is also the
# set every quantitative analysis downstream uses.
cat "$NTG_VEH" "$HTAU_VEH" "$HTAU_915" \
  | sortBed -i - | mergeBed -i - > union_relaxed3.bed
cat "$NTG_VEH" "$NTG_915" "$HTAU_VEH" "$HTAU_915" \
  | sortBed -i - | mergeBed -i - > union_relaxed4.bed
printf "union_relaxed3.bed (Fig. 4c)      %9d regions\n" "$(wc -l < union_relaxed3.bed)"
printf "union_relaxed4.bed (everything else) %6d regions\n" "$(wc -l < union_relaxed4.bed)"
# Expected: 295,426 and 299,069, which are the numbers the manuscript reports
# and the sets the published panels were built on. Both include the single peak
# called on chrM; excluding it gives 295,425 and 299,068. The Methods must not
# say "after exclusion of chrM" next to these two numbers.
# 11_ rebuilds the four-way union in R with GenomicRanges::reduce, which is the
# same operation; the two agree.

# ---- 2. signal over each union ----------------------------------------------
computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
    -R union_relaxed3.bed \
    -S nTg_Veh_avg.bw hTau_Veh_avg.bw hTau_915_avg.bw \
    --missingDataAsZero \
    -o matrix_union_relaxed3.gz -p "$THREADS"

computeMatrix reference-point --referencePoint center -b 3000 -a 3000 \
    -R union_relaxed4.bed \
    -S nTg_Veh_avg.bw nTg_915_avg.bw hTau_Veh_avg.bw hTau_915_avg.bw \
    --missingDataAsZero \
    -o matrix_union_relaxed4.gz -p "$THREADS"

# ---- 3. Figure 4c ------------------------------------------------------------
# Sorted on the nTg vehicle track, so the row order is set by the reference
# condition and not by the condition whose change is being claimed.
plotHeatmap -m matrix_union_relaxed3.gz \
    --samplesLabel "nTg Veh" "hTau Veh" "hTau RA013915" \
    --regionsLabel "Consensus peaks" \
    --colorMap inferno --sortRegions descend --sortUsingSample 1 \
    --heatmapHeight 15 --interpolationMethod bilinear --dpi 300 \
    -o Fig4c_heatmap_union_3cond.pdf

# ---- 4. Extended Data Fig. 5d ------------------------------------------------
plotHeatmap -m matrix_union_relaxed4.gz \
    --samplesLabel "nTg Veh" "nTg RA013915" "hTau Veh" "hTau RA013915" \
    --regionsLabel "Consensus peaks" \
    --colorMap inferno --sortRegions descend --sortUsingSample 1 \
    --heatmapHeight 15 --interpolationMethod bilinear --dpi 300 \
    -o FigED5d_heatmap_union_4cond.pdf

echo "written: union_relaxed3.bed, union_relaxed4.bed,"
echo "         Fig4c_heatmap_union_3cond.pdf, FigED5d_heatmap_union_4cond.pdf"
