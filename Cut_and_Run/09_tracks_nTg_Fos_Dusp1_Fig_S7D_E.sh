#!/bin/bash
# ==============================================================================
# Script 09: Genome browser tracks at the Fos and Dusp1 loci, nTg samples
# Figure: S7D (Fos), S7E (Dusp1)
#
# Requires pyGenomeTracks .ini config files (see tracks_nTg_Fos.ini and
# tracks_nTg_Dusp1.ini in ini_configs/) pointing to the nTg-Veh-H3K4ac-avg.bw /
# nTg-RA013915-H3K4ac-avg.bw tracks, IgG control, and consensus peak BED
# files, plus the mm10 gene annotation GTF.
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-./cutrun_work}"
cd "$WORKDIR"

pyGenomeTracks --tracks tracks_nTg_Fos.ini \
    --region chr12:85473000-85478000 \
    --outFileName Fos_nTg_track.pdf \
    --width 20 --dpi 300

pyGenomeTracks --tracks tracks_nTg_Dusp1.ini \
    --region chr17:26558000-26568000 \
    --outFileName Dusp1_nTg_track.pdf \
    --width 20 --dpi 300
