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

# ---- Feature distribution across the 3 main conditions ----
peak_files <- list(
  "nTg Veh" = nTg_veh,
  "hTau Veh" = hTau_veh,
  "hTau 915" = hTau_915
)

peakAnnoList <- lapply(peak_files, function(f) {
  annotatePeak(f, TxDb = txdb, tssRegion = c(-3000, 3000),
               annoDb = "org.Mm.eg.db", verbose = FALSE)
})

pdf("feature_distribution.pdf", width = 10, height = 4)
plotAnnoBar(peakAnnoList)
dev.off()

for (cond in names(peakAnnoList)) {
  cat("\n===", cond, "===\n")
  print(peakAnnoList[[cond]]@annoStat)
}

# ---- Feature distribution of the DEG-up ("rescued loci") peaks ----
anno_degup <- annotatePeak(deg_up_peaks, TxDb = txdb, tssRegion = c(-3000, 3000),
                            annoDb = "org.Mm.eg.db", verbose = FALSE)

n_degup <- length(anno_degup@anno)
cat(sprintf("=== DEG-up peaks (n=%d) ===\n", n_degup))
print(anno_degup@annoStat)

pdf("feature_distribution_DEGup.pdf", width = 8, height = 4)
plotAnnoBar(anno_degup)
dev.off()

# ---- Comparative barplot: genome-wide peaks vs. DEG-up (rescued) peaks (Fig. 4D) ----
anno_list <- lapply(
  list(
    "All peaks\n(genome-wide)" = hTau_915,
    "DEG-up peaks\n(rescued loci)" = deg_up_peaks
  ),
  function(f) {
    annotatePeak(f, TxDb = txdb, tssRegion = c(-3000, 3000),
                 annoDb = "org.Mm.eg.db", verbose = FALSE)
  }
)

pdf("feature_distribution_comparison.pdf", width = 10, height = 5)
plotAnnoBar(anno_list)
dev.off()

cat("Done!\n")
