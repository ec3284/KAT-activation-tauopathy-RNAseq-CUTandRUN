#!/bin/bash
# ==============================================================================
# Script 07: Genome browser tracks at the Fos and Rhob loci, hTau samples
# Figure: 4H (Fos), 4I (Rhob)
#
# Requires pyGenomeTracks .ini config files (see tracks_hTau_Fos.ini
# and tracks_hTau_Rhob.ini in ini_configs/) pointing to the
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

pyGenomeTracks --tracks tracks_hTau_Rhob.ini \
    --region chr12:8497000-8501000 \
    --outFileName Rhob_hTau_track.pdf \
    --width 20 --dpi 300
