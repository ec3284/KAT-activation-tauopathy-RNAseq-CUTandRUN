#!/bin/bash
# ==============================================================================
# Genome browser tracks at the Fos and Dusp1 loci, hTau conditions
# Figure 4f (Fos) and Figure 4g (Dusp1)
#
# These are the two loci the RNA-seq and the CUT&RUN have in common: Fos and
# Dusp1 are among the genes the compound restores in the transcriptome, and
# this panel shows the H3K4ac signal over the same promoters.
#
# The tracks are the spike-in normalized per-condition averages, and the .ini
# files fix min_value = 0 and the same max_value on every track of a panel
# (3 at Fos, 2.2 at Dusp1), IgG included. Keep it that way: a track left on its
# own autoscale shows a difference that is the autoscale.
#
# Requires: pyGenomeTracks, and the .ini files in this folder
#
# Input : tracks_hTau_Fos.ini, tracks_hTau_Dusp1.ini, which point at
#         nTg_Veh_avg.bw / hTau_Veh_avg.bw / hTau_915_avg.bw from 01_, the IgG
#         track, the consensus peak BED files and the mm10 GTF
# Output: Fig4f_Fos_hTau_track.pdf, Fig4g_Dusp1_hTau_track.pdf
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-.}"
cd "$WORKDIR"
for i in tracks_hTau_Fos.ini tracks_hTau_Dusp1.ini; do
    [ -f "$i" ] || { echo "missing config: $i" >&2; exit 1; }
done

pyGenomeTracks --tracks tracks_hTau_Fos.ini \
    --region chr12:85473000-85478000 \
    --outFileName Fig4f_Fos_hTau_track.pdf \
    --width 20 --dpi 300

pyGenomeTracks --tracks tracks_hTau_Dusp1.ini \
    --region chr17:26558000-26568000 \
    --outFileName Fig4g_Dusp1_hTau_track.pdf \
    --width 20 --dpi 300

echo "written: Fig4f_Fos_hTau_track.pdf, Fig4g_Dusp1_hTau_track.pdf"
