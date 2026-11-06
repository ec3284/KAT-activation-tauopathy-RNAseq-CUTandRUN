#!/bin/bash
# ==============================================================================
# Per-condition averaged tracks, and the deepTools matrix behind Fig. 4b
#
# Everything downstream reads one of the two kinds of file this script writes:
# the four averaged bigWigs (used by the metaplot, the heatmap and the genome
# browser tracks) and the twelve-library TSS matrix (used by the promoter AUC
# panel).
#
# Two steps were removed from the earlier version of this script and are not
# coming back. It built a transcript-to-gene table and then collapsed the three
# biological replicates to a group median before testing:
#
#     wide = long.groupby(["Gene", "Group"])["log2AUC"].median().unstack()
#     stat, p = wilcoxon(tmp[a], tmp[b])
#
# From that line onward the individual animals no longer exist, so the only
# test the file can support is a test across genes, which treats co-regulated
# promoters as independent replicates. The experimental unit is the mouse, and
# that analysis now lives in 06_promoter_AUC_per_mouse_Fig4B_FigED5C.R, which
# reads the matrix written here and keeps the twelve libraries separate.
#
# Requires: deepTools (bigwigAverage, computeMatrix). Spike-in normalization is
# already carried by the input bigWigs; nothing here renormalizes anything.
#
# Input : the twelve per-replicate spike-in normalized bigWigs and the mm10 GTF
#         (both as deposited in GEO)
# Output: nTg_Veh_avg.bw, nTg_915_avg.bw, hTau_Veh_avg.bw, hTau_915_avg.bw
#         AUC_PROMOTER_PC.matrix_ALL12_PC.noSkip.gz
#         AUC_PROMOTER_PC.samples.txt
# ==============================================================================
set -euo pipefail

# ---- PATHS -------------------------------------------------------------------
# Run from the directory where the outputs should be written, or set WORKDIR.
WORKDIR="${WORKDIR:-.}"
GTF="${GTF:-Mus_musculus.GRCm38.102.gtf.gz}"
THREADS="${THREADS:-max}"
cd "$WORKDIR"
[ -f "$GTF" ] || { echo "GTF not found: $GTF  (set GTF=/path/to/gtf.gz)" >&2; exit 1; }

# ---- 1. average the three replicates of each condition -----------------------
# The output names are the ones the R scripts expect. If your per-replicate
# files are named differently, change the -b lists and nothing else.
bigwigAverage -b nTg-Veh-H3K4ac-rep1.bw nTg-Veh-H3K4ac-rep2.bw nTg-Veh-H3K4ac-rep3.bw \
    -o nTg_Veh_avg.bw  -p "$THREADS"
bigwigAverage -b nTg-RA013915-H3K4ac-rep1.bw nTg-RA013915-H3K4ac-rep2.bw nTg-RA013915-H3K4ac-rep3.bw \
    -o nTg_915_avg.bw  -p "$THREADS"
bigwigAverage -b hTau-Veh-H3K4ac-rep1.bw hTau-Veh-H3K4ac-rep2.bw hTau-Veh-H3K4ac-rep3.bw \
    -o hTau_Veh_avg.bw -p "$THREADS"
bigwigAverage -b hTau-RA013915-H3K4ac-rep1.bw hTau-RA013915-H3K4ac-rep2.bw hTau-RA013915-H3K4ac-rep3.bw \
    -o hTau_915_avg.bw -p "$THREADS"

# ---- 2. TSS matrix over all twelve libraries ---------------------------------
# The ORDER of the -S list defines the column order of the matrix, and
# 06_promoter_AUC_per_mouse_Fig4B_FigED5C.R reads the columns by position
# through its SAMPLES vector. The two must agree. SAMPLES there is
#   WT_Veh_R1..R3, WT_915_R1..R3, hTau_Veh_R1..R3, hTau_915_R1..R3
# which is the order below. Do not reorder one without the other.
#
# --skipZeros drops regions that are zero in every library. The file name says
# "noSkip" for historical reasons and is kept because that is the name the
# published analysis read; the flag, not the name, is what applies.
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
    -R "$GTF" \
    -S nTg-Veh-H3K4ac-rep1.bw  nTg-Veh-H3K4ac-rep2.bw  nTg-Veh-H3K4ac-rep3.bw \
       nTg-RA013915-H3K4ac-rep1.bw nTg-RA013915-H3K4ac-rep2.bw nTg-RA013915-H3K4ac-rep3.bw \
       hTau-Veh-H3K4ac-rep1.bw hTau-Veh-H3K4ac-rep2.bw hTau-Veh-H3K4ac-rep3.bw \
       hTau-RA013915-H3K4ac-rep1.bw hTau-RA013915-H3K4ac-rep2.bw hTau-RA013915-H3K4ac-rep3.bw \
    --binSize 10 --skipZeros \
    -o AUC_PROMOTER_PC.matrix_ALL12_PC.noSkip.gz -p "$THREADS"

printf "%s\n" \
    WT_Veh_R1 WT_Veh_R2 WT_Veh_R3 \
    WT_915_R1 WT_915_R2 WT_915_R3 \
    hTau_Veh_R1 hTau_Veh_R2 hTau_Veh_R3 \
    hTau_915_R1 hTau_915_R2 hTau_915_R3 \
    > AUC_PROMOTER_PC.samples.txt

echo "written: four *_avg.bw tracks, AUC_PROMOTER_PC.matrix_ALL12_PC.noSkip.gz"
