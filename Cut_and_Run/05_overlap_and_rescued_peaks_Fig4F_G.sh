#!/bin/bash
# ==============================================================================
# Script 05: Pairwise peak overlap (Fig. 4F, Venn diagram counts) and
#            identification of "rescued" peaks (Fig. 4G)
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-./cutrun_work}"
cd "$WORKDIR"

# ---- Pairwise overlap counts, used to draw the Venn diagram (Fig. 4F) ----
echo "nTg_Veh vs hTau_Veh:"
intersectBed -a nTg-Veh-H3K4ac-consensus.peaks.bed -b hTau-Veh-H3K4ac-consensus.peaks.bed -u | wc -l

echo "nTg_Veh vs hTau_915:"
intersectBed -a nTg-Veh-H3K4ac-consensus.peaks.bed -b hTau-RA013915-H3K4ac-consensus.peaks.bed -u | wc -l

echo "hTau_Veh vs hTau_915:"
intersectBed -a hTau-Veh-H3K4ac-consensus.peaks.bed -b hTau-RA013915-H3K4ac-consensus.peaks.bed -u | wc -l

# NOTE: the Venn diagram figure itself (proportional ellipses/labels) was
# drawn manually from these counts; see 08_venn_nTg_Fig_S7C.py for a
# template plotting script (values must be filled in from the counts above).

# ---- Rescued peaks: present in nTg_Veh AND hTau_915, but absent in hTau_Veh ----
# i.e., loci lost in disease (hTau_Veh) and restored by RA013915 (hTau_915)
intersectBed -a nTg-Veh-H3K4ac-consensus.peaks.bed -b hTau-RA013915-H3K4ac-consensus.peaks.bed -u \
    > tmp_nTg_hTauVeh_relaxed.bed
wc -l < tmp_nTg_hTauVeh_relaxed.bed

echo "=== Rescued peaks ==="
intersectBed -a nTg-Veh-H3K4ac-consensus.peaks.bed -b hTau-Veh-H3K4ac-consensus.peaks.bed -v \
    > lost_in_hTauVeh_relaxed.bed
intersectBed -a lost_in_hTauVeh_relaxed.bed -b hTau-RA013915-H3K4ac-consensus.peaks.bed -u \
    > rescued_peaks_relaxed.bed
wc -l rescued_peaks_relaxed.bed
