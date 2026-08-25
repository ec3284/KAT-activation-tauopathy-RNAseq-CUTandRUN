#!/bin/bash
# ==============================================================================
# Script 07: Genome browser tracks at the Fos and Dusp1 loci, hTau samples
# Figure: 4H (Fos), 4I (Dusp1)
#
# Requires pyGenomeTracks .ini config files (see tracks_hTau_Fos.ini
# and tracks_hTau_Dusp1.ini in ini_configs/) pointing to the
# nTg-Veh-H3K4ac-avg.bw / hTau-Veh-H3K4ac-avg.bw / hTau-RA013915-H3K4ac-avg.bw
# tracks, IgG control, and consensus peak BED files, plus the mm10 gene
# annotation GTF.
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-./cutrun_work}"
cd "$WORKDIR"

pyGenomeTracks --tracks tracks_hTau_Fos.ini \
    --region chr12:85473000-85478000 \
    --outFileName Fos_hTau_track.pdf \
    --width 20 --dpi 300

pyGenomeTracks --tracks tracks_hTau_Dusp1.ini \
    --region chr17:26558000-26568000 \
    --outFileName Dusp1_hTau_track.pdf \
    --width 20 --dpi 300
