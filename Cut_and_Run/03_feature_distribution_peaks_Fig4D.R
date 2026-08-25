# ==============================================================================
# Script 03: Genomic feature distribution of H3K4ac peaks
# Figure: 4D
# ==============================================================================

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# ---- Input peak files ----
# nTg_veh / hTau_veh / hTau_915: consensus peak BED files deposited in GEO
# deg_up: peaks overlapping RA013915-upregulated genes ("DEG-up peaks"),
# output of 04_RNAseq_integration_Fig4E.sh (union_relaxed_on_DEGup.bed)
nTg_veh <- "nTg-Veh-H3K4ac-consensus.peaks.bed"
hTau_veh <- "hTau-Veh-H3K4ac-consensus.peaks.bed"
hTau_915 <- "hTau-RA013915-H3K4ac-consensus.peaks.bed"
deg_up_peaks <- "union_relaxed_on_DEGup.bed"
stopifnot(file.exists(deg_up_peaks))

# ---- Feature distribution: 3 conditions + DEG-up peaks, single combined plot (Fig. 4D) ----
peak_files <- list(
  "nTg Veh" = nTg_veh,
  "hTau Veh" = hTau_veh,
  "hTau 915" = hTau_915,
  "DEG-up peaks" = deg_up_peaks
)

peakAnnoList <- lapply(peak_files, function(f) {
  annotatePeak(f, TxDb = txdb, tssRegion = c(-3000, 3000),
               annoDb = "org.Mm.eg.db", verbose = FALSE)
})

pdf("feature_distribution_Fig4D.pdf", width = 10, height = 4)
plotAnnoBar(peakAnnoList)
dev.off()

for (cond in names(peakAnnoList)) {
  cat("\n===", cond, "(n=", length(peakAnnoList[[cond]]@anno), ") ===\n")
  print(peakAnnoList[[cond]]@annoStat)
}

cat("Done!\n")
