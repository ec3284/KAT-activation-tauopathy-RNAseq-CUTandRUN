#!/bin/bash
# ==============================================================================
# Genome browser tracks at the Fos and Dusp1 loci, nTg conditions
# Extended Data Fig. 5e (Fos) and Extended Data Fig. 5f (Dusp1)
#
# The same two loci as Fig. 4f,g, in the non-transgenic animals. This is the
# panel that shows the compound does not act on these promoters in the absence
# of the tauopathy, and it is worth reading next to one number from the peak
# calls: at the Fos promoter the nTg RA013915 group carries the highest read
# count of the four groups and is the only group with no peak called there. The
# tracks show the signal; the absence of a call is a property of the threshold.
#
# Requires: pyGenomeTracks, and the .ini files in this folder
#
# Input : tracks_nTg_Fos.ini, tracks_nTg_Dusp1.ini, which point at
#         nTg_Veh_avg.bw / nTg_915_avg.bw from 01_, the IgG track, the
#         consensus peak BED files and the mm10 GTF
# Output: FigED5e_Fos_nTg_track.pdf, FigED5f_Dusp1_nTg_track.pdf
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-.}"
cd "$WORKDIR"
for i in tracks_nTg_Fos.ini tracks_nTg_Dusp1.ini; do
    [ -f "$i" ] || { echo "missing config: $i" >&2; exit 1; }
done

pyGenomeTracks --tracks tracks_nTg_Fos.ini \
    --region chr12:85473000-85478000 \
    --outFileName FigED5e_Fos_nTg_track.pdf \
    --width 20 --dpi 300

pyGenomeTracks --tracks tracks_nTg_Dusp1.ini \
    --region chr17:26558000-26568000 \
    --outFileName FigED5f_Dusp1_nTg_track.pdf \
    --width 20 --dpi 300

echo "written: FigED5e_Fos_nTg_track.pdf, FigED5f_Dusp1_nTg_track.pdf"
