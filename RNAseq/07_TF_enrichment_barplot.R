#!/usr/bin/env Rscript
# ==============================================================================
# Script 07: Transcription factor enrichment barplot (Figure 3E)
#
# Values below were obtained from Enrichr (maayanlab.cloud/Enrichr) using the
# 22 validated rescue genes (rescued_up + rescued_down from
# 01_DESeq2_collapsed_rescue_analysis.R) as input. Three independent
# databases were queried under the Transcription tab: "ChEA ChIP-seq",
# "ARCHS4 TFs Coexp", and "TRRUST Transcription Factors 2019". A fourth
# database ("Motif Analysis") was queried but excluded from the main figure
# panel: its top hits (ELK1, SP1, AR, POU5F1, RUNX1) were weaker (all
# padj>0.04) and less biologically coherent with the rescue signature than
# the other three databases, which converge on AP-1/EP300/CREM.
#
# The top 5 terms per database (by adjusted p-value) are plotted. Note: for
# ChEA and ARCHS4, several terms beyond rank 5 share an identical p-value
# with rank 5 (e.g., KLF2/KLF4/KLF5 in ChEA; JUNB in ARCHS4) due to
# identical gene-overlap counts in the underlying ChIP-seq/coexpression
# datasets; these tied terms are omitted here for a concise 5-per-database
# panel, but are documented in the full 10-per-database table (see
# TF_barplot_4sources_ALL.pdf / source data) for transparency.
#
# This script only builds the plot from the recorded Enrichr output; it
# does not re-run the enrichment (Enrichr is a web tool, not an R package).
# To regenerate/verify: submit the 22 rescue gene symbols to Enrichr and
# export results from each of the four libraries listed above.
# ==============================================================================

library(ggplot2)
library(dplyr)

WORKDIR <- Sys.getenv("WORKDIR", unset = "./rnaseq_work")
setwd(WORKDIR)

# ---- 1. Full results, all 4 databases, top 10 terms each (source data /
#         supplementary; includes ties beyond rank 5) ----
tf_all_top10 <- data.frame(
  TF = c(
    # Motif Analysis (top 5 only; weak/incoherent signal, no ties recorded)
    "ELK1","SP1","AR","POU5F1","RUNX1",
    # ChEA ChIP-seq (top 10; ranks 7-10 tied with rank 6/7)
    "LEF1","EP300","AF4","GATA1","SMC1","FOXO1","KLF5","KLF2","KLF4_a","KLF4_b",
    # ARCHS4 TF Coexp (top 10; ranks 6-10 include ties)
    "FOSB","NR4A2","JUN","KLF6","ZFP36","JUNB","ZFP36L1","TSC22D3","NR4A1","JUND",
    # TRRUST (top 10; no ties)
    "TRP53","CREM","TP53","STAT6","CEBPB","HDAC1","ESR1","STAT3","ZBTB7A","TBP"
  ),
  Padj = c(
    0.04196, 0.04196, 0.04196, 0.04426, 0.04426,
    0.0004700, 0.002459, 0.002595, 0.002595, 0.003479, 0.003479, 0.003586, 0.003586, 0.003586, 0.003586,
    7.64e-7, 7.64e-7, 7.64e-7, 1.505e-5, 1.505e-5, 1.505e-5, 3.514e-4, 3.514e-4, 5.206e-3, 5.206e-3,
    0.00002429, 0.009182, 0.009584, 0.009584, 0.01045, 0.02291, 0.02291, 0.02291, 0.03393, 0.03393
  ),
  OddsRatio = c(
    49.90, 14.22, 13.79, 11.00, 22.49,
    11.94, 11.45, 14.97, 10.26, 9.35, 9.25, 48.54, 48.54, 48.54, 9.69,
    42.91, 42.91, 42.91, 33.60, 33.60, 33.60, 25.76, 25.76, 19.07, 19.07,
    47.06, 108.47, 24.62, 73.34, 62.32, 36.07, 33.63, 32.31, 235.02, 195.84
  ),
  Source = c(rep("Motif Analysis", 5), rep("ChEA ChIP-seq", 10),
             rep("ARCHS4 Coexp", 10), rep("TRRUST", 10))
)

tf_all_top10$Source <- factor(tf_all_top10$Source,
                                levels = c("Motif Analysis","ChEA ChIP-seq","ARCHS4 Coexp","TRRUST"))

p_barplot_tf_all <- ggplot(tf_all_top10, aes(x = reorder(TF, -log10(Padj)), y = -log10(Padj), fill = Source)) +
  geom_col() +
  facet_wrap(~Source, scales = "free", ncol = 1) +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(size = 12, face = "bold"), legend.position = "none") +
  labs(x = NULL, y = "-log10(adjusted p-value)",
       title = "Transcription factor enrichment, all databases (source data)\nRescue genes (n=22)")

print(p_barplot_tf_all)
ggsave("TF_barplot_4sources_ALL.pdf", p_barplot_tf_all, width = 8, height = 16, dpi = 300)

# ---- 2. Main figure panel: top 5 per database, 3 strongest/most coherent
#         sources (Motif Analysis excluded, see rationale above) ----
tf_3sources_top5 <- data.frame(
  TF = c("LEF1","EP300","AF4","GATA1","SMC1",
         "FOSB","NR4A2","JUN","KLF6","ZFP36",
         "TRP53","CREM","TP53","STAT6","CEBPB"),
  Padj = c(0.0004700, 0.002459, 0.002595, 0.002595, 0.003479,
           7.64e-7, 7.64e-7, 7.64e-7, 1.505e-5, 1.505e-5,
           0.00002429, 0.009182, 0.009584, 0.009584, 0.01045),
  OddsRatio = c(11.94, 11.45, 14.97, 10.26, 9.35,
                42.91, 42.91, 42.91, 33.60, 33.60,
                47.06, 108.47, 24.62, 73.34, 62.32),
  Source = rep(c("ChEA ChIP-seq","ARCHS4 Coexp","TRRUST"), each = 5)
)
tf_3sources_top5$Source <- factor(tf_3sources_top5$Source,
                                     levels = c("ChEA ChIP-seq","ARCHS4 Coexp","TRRUST"))

tf_colors_blue <- c("ChEA ChIP-seq" = "#1565C0", "ARCHS4 Coexp" = "#4FC3F7", "TRRUST" = "#0D47A1")

p_barplot_tf3 <- ggplot(tf_3sources_top5, aes(x = reorder(TF, -log10(Padj)), y = -log10(Padj), fill = Source)) +
  geom_col() +
  facet_wrap(~Source, scales = "free", ncol = 1) +
  coord_flip() +
  scale_fill_manual(values = tf_colors_blue) +
  theme_bw(base_size = 14) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        strip.text = element_text(size = 13, face = "bold"),
        legend.position = "none") +
  labs(x = NULL, y = "-log10(adjusted p-value)",
       title = "Transcription factor enrichment\nRescue genes (n=22)")

print(p_barplot_tf3)
ggsave("TF_barplot_3sources_FINAL.pdf", p_barplot_tf3, width = 6, height = 12, dpi = 300)

cat("Done! Saved in", WORKDIR, ":\n",
    "- TF_barplot_4sources_ALL.pdf (all 4 databases, top 10 each, source data)\n",
    "- TF_barplot_3sources_FINAL.pdf (top 5, 3 sources, main figure panel)\n")
