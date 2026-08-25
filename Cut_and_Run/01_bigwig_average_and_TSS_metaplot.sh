#!/bin/bash
# ==============================================================================
# Script 01: BigWig averaging across replicates + TSS metaplot
# Figures: 4A, 4B, S7A
#
# Usage: run from a directory containing per-replicate spike-in-normalized
# bigWig files (same files deposited in GEO) and the mm10 GTF.
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-./cutrun_work}"
GTF="${GTF:-${WORKDIR}/Mus_musculus.GRCm38.102.gtf.gz}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ---- 1. Average bigWig tracks across biological replicates, per condition ----
bigwigAverage -b nTg-Veh-H3K4ac-rep1.bw nTg-Veh-H3K4ac-rep2.bw nTg-Veh-H3K4ac-rep3.bw \
    -o nTg-Veh-H3K4ac-avg.bw -p max
bigwigAverage -b nTg-RA013915-H3K4ac-rep1.bw nTg-RA013915-H3K4ac-rep2.bw nTg-RA013915-H3K4ac-rep3.bw \
    -o nTg-RA013915-H3K4ac-avg.bw -p max
bigwigAverage -b hTau-Veh-H3K4ac-rep1.bw hTau-Veh-H3K4ac-rep2.bw hTau-Veh-H3K4ac-rep3.bw \
    -o hTau-Veh-H3K4ac-avg.bw -p max
bigwigAverage -b hTau-RA013915-H3K4ac-rep1.bw hTau-RA013915-H3K4ac-rep2.bw hTau-RA013915-H3K4ac-rep3.bw \
    -o hTau-RA013915-H3K4ac-avg.bw -p max
# Output: nTg-Veh-H3K4ac-avg.bw, nTg-RA013915-H3K4ac-avg.bw,
#         hTau-Veh-H3K4ac-avg.bw, hTau-RA013915-H3K4ac-avg.bw

# ---- 2. computeMatrix at TSS +/- 3kb (Fig. 4A-B, S7A) ----
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
    -R "$GTF" \
    -S nTg-Veh-H3K4ac-avg.bw hTau-Veh-H3K4ac-avg.bw hTau-RA013915-H3K4ac-avg.bw \
    --missingDataAsZero \
    -o matrix_TSS_avg3_zeros.gz -p max

# Metaplot, Figure 4A
plotProfile -m matrix_TSS_avg3_zeros.gz \
    --colors "#CFCFCF" "#FFE0C0" "#FF6A00" \
    --plotTitle "H3K4ac at TSS" \
    --samplesLabel "nTg Veh" "hTau Veh" "hTau 915" \
    --perGroup \
    -o metaplot_TSS_3cond.pdf \
    --plotWidth 12 --plotHeight 8 --plotFileFormat pdf --dpi 300
