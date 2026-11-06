#!/bin/bash
# ==============================================================================
# Pairwise peak overlap and the rescued peak set
# Figure 4e (Venn), and the input to Fig. 4f
#
# A rescued peak is one called in vehicle nTg and in RA013915-treated hTau but
# not in vehicle hTau: lost in disease and present again under treatment.
#
# What the overlap counts can and cannot say. bedtools counts the peaks of one
# set that overlap the other, and that count is not symmetric: a region split
# into two peaks in one set and merged into one in the other contributes
# differently from each side. A two-circle Venn has one number for the
# intersection and therefore cannot be consistent with both totals at once.
# Both sides are printed below, and each circle is checked against the
# consensus total it came from, so the figure cannot be drawn from a pair of
# numbers that do not add up.
#
# Requires: bedtools
#
# Input : the four consensus peak BED files
# Output: lost_in_hTau_Veh.bed, rescued_peaks.bed, Fig4e_overlap_counts.txt
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-.}"
PEAKS="${CUTRUN_PEAKS:-peaks_and_consensus_peaks}"
cd "$WORKDIR"

NTG_VEH="$PEAKS/nTg-Veh-H3K4ac-consensus.peaks.bed"
NTG_915="$PEAKS/nTg-RA013915-H3K4ac-consensus.peaks.bed"
HTAU_VEH="$PEAKS/hTau-Veh-H3K4ac-consensus.peaks.bed"
HTAU_915="$PEAKS/hTau-RA013915-H3K4ac-consensus.peaks.bed"
for f in "$NTG_VEH" "$NTG_915" "$HTAU_VEH" "$HTAU_915"; do
    [ -f "$f" ] || { echo "missing: $f" >&2; exit 1; }
done

pair () {   # pair <label A> <bed A> <label B> <bed B>
    local la=$1 a=$2 lb=$3 b=$4
    local na nb ab ba
    na=$(wc -l < "$a"); nb=$(wc -l < "$b")
    ab=$(intersectBed -a "$a" -b "$b" -u | wc -l)
    ba=$(intersectBed -a "$b" -b "$a" -u | wc -l)
    printf "%-14s %9d total   %9d overlap %-14s   %9d only\n" "$la" "$na" "$ab" "$lb" "$((na - ab))"
    printf "%-14s %9d total   %9d overlap %-14s   %9d only\n" "$lb" "$nb" "$ba" "$la" "$((nb - ba))"
    printf "   the two intersection counts differ by %d; a Venn drawn with one of\n" "$((ab > ba ? ab - ba : ba - ab))"
    printf "   them will not reproduce the other circle's total.\n\n"
}

{
  echo "=== pairwise overlap, both directions ==="
  echo
  pair "nTg Veh"  "$NTG_VEH"  "hTau Veh"      "$HTAU_VEH"
  pair "nTg Veh"  "$NTG_VEH"  "hTau RA013915" "$HTAU_915"
  pair "hTau Veh" "$HTAU_VEH" "hTau RA013915" "$HTAU_915"
  pair "nTg Veh"  "$NTG_VEH"  "nTg RA013915"  "$NTG_915"
} | tee Fig4e_overlap_counts.txt

# ---- rescued peaks -----------------------------------------------------------
intersectBed -a "$NTG_VEH" -b "$HTAU_VEH" -v > lost_in_hTau_Veh.bed
intersectBed -a lost_in_hTau_Veh.bed -b "$HTAU_915" -u > rescued_peaks.bed
printf "lost in vehicle hTau   %9d\n" "$(wc -l < lost_in_hTau_Veh.bed)"
printf "of these, rescued      %9d\n" "$(wc -l < rescued_peaks.bed)"

# 09_GO_rescued_peaks_Fig4F.R rebuilds this same set in R from the same four
# BED files, with the same definition, and writes its own rescued_peaks.bed.
# The two are expected to agree; if they do not, the BED files are not the ones
# this script read.
