#!/usr/bin/env Rscript
# ==============================================================================
# Script 02: MA plots for the two key RNA-seq comparisons
#
# Requires objects from 01_DESeq2_collapsed_rescue_analysis.R to be in memory:
#   - res_nTg          (results object, nTg 915 vs nTg Veh)
#   - res_treatment     (results object, hTau 915 vs hTau Veh)
#   - rescued_up, rescued_down  (data frames of validated rescue genes)
#
# If starting a fresh R session, re-run script 01 first, or reload the saved
# workspace / re-derive these objects from RNAseq_FINAL_DESeq2_results.xlsx.
# ==============================================================================

WORKDIR <- Sys.getenv("WORKDIR", unset = "./rnaseq_work")
setwd(WORKDIR)

# ---- nTg: Vehicle vs RA013915 ----
# Expected/confirmed result: no genes reach padj<0.05, consistent with the
# chromatin-gated model (KAT activation has minimal transcriptional effect
# in mice with intact chromatin acetylation homeostasis).
png("MAplot_nTg_VehVs915.png", width = 1400, height = 1100, res = 150)
plotMA(res_nTg, ylim = c(-3, 3), main = "nTg: Vehicle vs RA013915")
dev.off()

# ---- hTau: Vehicle vs RA013915 ----
# Blue = statistically significant DEGs (padj<0.05), matching the Figure S6C
# legend. Rescue genes (disease+treatment cross-filter) are already
# highlighted and labeled separately in the volcano plot (Figure 3A);
# re-highlighting them here would be redundant, so this MA plot shows only
# the significant/non-significant split.
ma_df <- as.data.frame(res_treatment)
ma_df$gene <- rownames(ma_df)
ma_df$isSigTreatment <- !is.na(ma_df$padj) & ma_df$padj < 0.05

png("MAplot_hTau_BLUE.png", width = 1400, height = 1100, res = 150)
plot(ma_df$baseMean, ma_df$log2FoldChange, log = "x", pch = 19, cex = 0.4, col = "grey70",
     xlab = "mean of normalized counts", ylab = "log fold change",
     main = "hTau: Vehicle vs RA013915", ylim = c(-3, 3))
abline(h = 0, col = "grey40", lwd = 2)
points(ma_df$baseMean[ma_df$isSigTreatment], ma_df$log2FoldChange[ma_df$isSigTreatment],
       col = "#3B6FA0", pch = 19, cex = 0.5)
legend("bottomright", legend = "Significant (padj < 0.05)",
       col = "#3B6FA0", pch = 19, cex = 0.8, bty = "n")
dev.off()

cat("Done! Saved in", WORKDIR, ":\n",
    "- MAplot_nTg_VehVs915.png\n",
    "- MAplot_hTau_BLUE.png\n")
