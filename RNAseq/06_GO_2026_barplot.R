#!/usr/bin/env Rscript
# ==============================================================================
# Script 06: GO Biological Process 2026 enrichment barplot (Figure 3D)
#
# Values below were obtained from Enrichr (maayanlab.cloud/Enrichr),
# "GO Biological Process 2026" library, using the 22 validated rescue genes
# (rescued_up + rescued_down from 01_DESeq2_collapsed_rescue_analysis.R) as
# input. This script only builds the plot from the recorded Enrichr output;
# it does not re-run the enrichment (Enrichr is a web tool, not an R package).
#
# Bar height = -log10(adjusted p-value), consistent with the standard
# Enrichr bar chart display.
#
# To regenerate/verify these values: submit the 22 rescue gene symbols to
# Enrichr, select "GO Biological Process 2026" under the Ontologies tab,
# and export the full results table.
# ==============================================================================

library(ggplot2)

WORKDIR <- Sys.getenv("WORKDIR", unset = "./rnaseq_work")
setwd(WORKDIR)

# ---- 1. Enrichr results (GO Biological Process 2026), top 10 terms ----
enrichr_2026_data <- data.frame(
  Term = c("Positive regulation of transcription by RNA Pol II",
           "Response to laminar fluid shear stress",
           "Neuron projection fasciculation",
           "Cellular response to fluid shear stress",
           "Axonal fasciculation",
           "Positive regulation of DNA-templated transcription",
           "Regulation of anatomical structure morphogenesis",
           "Cellular response to laminar fluid shear stress",
           "Positive regulation of nitric oxide biosynthetic process",
           "Positive regulation of macromolecule metabolic process"),
  Padj = c(0.004391, 0.004391, 0.004391, 0.004391, 0.004723,
           0.005466, 0.005720, 0.01078, 0.01078, 0.01105)
)

enrichr_2026_data$Term <- factor(enrichr_2026_data$Term,
                                    levels = rev(enrichr_2026_data$Term[order(enrichr_2026_data$Padj, decreasing = TRUE)]))

# ---- 2. Barplot, Enrichr style ----
p_go2026_bar <- ggplot(enrichr_2026_data, aes(x = -log10(Padj), y = Term)) +
  geom_col(fill = "#FF6A00") +
  theme_bw(base_size = 13) +
  theme(axis.text.y = element_text(size = 11)) +
  labs(x = "-log10(adjusted p-value)", y = NULL,
       title = "GO Biological Process 2026 enrichment\nRescue genes (n=22)")

print(p_go2026_bar)
ggsave("GO_BP_2026_rescue_barplot.pdf", p_go2026_bar, width = 9, height = 6, dpi = 300)

cat("Done! GO_BP_2026_rescue_barplot.pdf saved in", WORKDIR, "\n")
