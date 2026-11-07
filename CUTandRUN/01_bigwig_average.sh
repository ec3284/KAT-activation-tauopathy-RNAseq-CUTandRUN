#!/bin/bash
# ==============================================================================
# Per-condition averaged tracks
#
# Everything downstream reads the four averaged bigWigs this script writes: the
# TSS metaplot, the heatmaps over the peak unions, and the genome browser
# tracks.
#
# Requires: deepTools (bigwigAverage). Spike-in normalization is already carried
# by the input bigWigs; nothing here renormalizes anything.
#
# Input : the twelve per-replicate spike-in normalized bigWigs (as deposited in
#         GEO)
# Output: nTg_Veh_avg.bw, nTg_915_avg.bw, hTau_Veh_avg.bw, hTau_915_avg.bw
# ==============================================================================
set -euo pipefail

# ---- PATHS -------------------------------------------------------------------
# Run from the directory where the outputs should be written, or set WORKDIR.
WORKDIR="${WORKDIR:-.}"
THREADS="${THREADS:-max}"
cd "$WORKDIR"

# ---- average the three replicates of each condition --------------------------
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

echo "written: nTg_Veh_avg.bw, nTg_915_avg.bw, hTau_Veh_avg.bw, hTau_915_avg.bw"
