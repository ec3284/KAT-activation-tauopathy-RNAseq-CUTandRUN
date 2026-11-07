#!/bin/bash
# ==============================================================================
# Raw read counts per mouse over the peak union, and the spike-in totals
#
# This is the step that turns twelve BAM files into the one table the
# differential analysis needs. It did not exist in the first version of this
# folder, because the first version compared conditions by counting called
# peaks rather than by counting reads.
#
# Counts are raw and stay raw. Normalization happens in
# 04_differential_peaks_spikein.R, where the spike-in totals enter DESeq2 as
# the size factors; anything scaled here would be normalized twice.
#
# Requires: bedtools (multicov) or featureCounts, and samtools
#
# Input : the twelve coordinate-sorted, indexed H3K4ac BAM files, and
#         union_relaxed4.bed from 02_
# Output: CUTRUN_union_peak_counts_12mice.csv.gz
#           peak_id, chrom, start, end, then one raw count column per mouse,
#           named exactly as the SPIKE vector in 04_ expects:
#           WT_Veh_R1..R3, WT_915_R1..R3, hTau_Veh_R1..R3, hTau_915_R1..R3
#         CUTRUN_spikein_reads.txt   reads mapped to the spike-in genome
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-.}"
BAMDIR="${CUTRUN_BAM:-bam}"
SPIKEDIR="${CUTRUN_SPIKE_BAM:-$BAMDIR}"
cd "$WORKDIR"
[ -f union_relaxed4.bed ] || { echo "run 02_ first: union_relaxed4.bed not found" >&2; exit 1; }
[ -d "$BAMDIR" ] || { echo "BAM directory not found: $BAMDIR  (set CUTRUN_BAM)" >&2; exit 1; }

# The twelve libraries in the column order 04_ expects, each as
# "<column name>:<bam file name>". Edit the file names, never the order.
LIBS=(
  "WT_Veh_R1:nTg-Veh-H3K4ac-rep1.target.markdup.sorted.bam"
  "WT_Veh_R2:nTg-Veh-H3K4ac-rep2.target.markdup.sorted.bam"
  "WT_Veh_R3:nTg-Veh-H3K4ac-rep3.target.markdup.sorted.bam"
  "WT_915_R1:nTg-RA013915-H3K4ac-rep1.target.markdup.sorted.bam"
  "WT_915_R2:nTg-RA013915-H3K4ac-rep2.target.markdup.sorted.bam"
  "WT_915_R3:nTg-RA013915-H3K4ac-rep3.target.markdup.sorted.bam"
  "hTau_Veh_R1:hTau-Veh-H3K4ac-rep1.target.markdup.sorted.bam"
  "hTau_Veh_R2:hTau-Veh-H3K4ac-rep2.target.markdup.sorted.bam"
  "hTau_Veh_R3:hTau-Veh-H3K4ac-rep3.target.markdup.sorted.bam"
  "hTau_915_R1:hTau-RA013915-H3K4ac-rep1.target.markdup.sorted.bam"
  "hTau_915_R2:hTau-RA013915-H3K4ac-rep2.target.markdup.sorted.bam"
  "hTau_915_R3:hTau-RA013915-H3K4ac-rep3.target.markdup.sorted.bam"
)

names=(); bams=()
for e in "${LIBS[@]}"; do
    n="${e%%:*}"; b="$BAMDIR/${e#*:}"
    [ -f "$b" ] || { echo "missing BAM: $b" >&2; exit 1; }
    [ -f "$b.bai" ] || samtools index "$b"
    names+=("$n"); bams+=("$b")
done
echo "counting ${#bams[@]} libraries over $(wc -l < union_relaxed4.bed) regions"

# ---- 1. counts ---------------------------------------------------------------
# multicov counts reads overlapping each region, in the order the BAMs are
# given, which is the order of LIBS above.
bedtools multicov -bams "${bams[@]}" -bed union_relaxed4.bed > union_counts_raw.tsv

# ---- 2. assemble the CSV -----------------------------------------------------
{
  printf "peak_id,chrom,start,end"
  for n in "${names[@]}"; do printf ",%s" "$n"; done
  printf "\n"
  awk -v OFS=, '{
      printf "peak_%d,%s,%s,%s", NR, $1, $2, $3
      for (i = 4; i <= NF; i++) printf ",%s", $i
      printf "\n"
  }' union_counts_raw.tsv
} | gzip -c > CUTRUN_union_peak_counts_12mice.csv.gz
rm -f union_counts_raw.tsv

# ---- 3. spike-in totals ------------------------------------------------------
# The numbers hard-coded as SPIKE in 04_. They come from the alignment against
# the spike-in genome, which nf-core/cutandrun reports per library; this block
# recomputes them from the spike-in BAMs so the two can be checked against each
# other. One library, hTau_915_R1, is far lower than the rest, which is why 04_
# runs a leave-one-out on it.
if compgen -G "$SPIKEDIR/*spikein*.bam" > /dev/null; then
  : > CUTRUN_spikein_reads.txt
  for e in "${LIBS[@]}"; do
      n="${e%%:*}"; stem="${e#*:}"; stem="${stem%%.*}"
      sb=$(compgen -G "$SPIKEDIR/${stem}*spikein*.bam" | head -1 || true)
      if [ -n "${sb:-}" ]; then
          printf "%s\t%s\n" "$n" "$(samtools view -c -F 260 "$sb")" >> CUTRUN_spikein_reads.txt
      else
          printf "%s\tNA\n" "$n" >> CUTRUN_spikein_reads.txt
      fi
  done
  cat CUTRUN_spikein_reads.txt
else
  echo "no spike-in BAMs found under $SPIKEDIR; take the totals from the" >&2
  echo "nf-core/cutandrun alignment report and check them against SPIKE in 04_." >&2
fi

echo "written: CUTRUN_union_peak_counts_12mice.csv.gz"
