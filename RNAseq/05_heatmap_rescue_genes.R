#!/usr/bin/env Rscript
# ==============================================================================
# Script 05: Heatmap of the 22 validated rescue genes, all 4 conditions
#
# Requires objects in memory from 01_DESeq2_collapsed_rescue_analysis.R:
#   - dds (the DESeq2 object, post-collapse, post-DESeq())
#   - rescued_up, rescued_down (validated rescue gene lists)
#
# Column annotation colors match the CUT&RUN pyGenomeTracks .ini configs
# (nTg Veh = #CFCFCF, nTg 915 = #056AE7, hTau Veh = #FFE0C0, hTau 915 = #FF6A00)
# for visual consistency across the RNA-seq and CUT&RUN figure panels.
# ==============================================================================

library(DESeq2)
library(pheatmap)
library(dplyr)

WORKDIR <- Sys.getenv("WORKDIR", unset = "./rnaseq_work")
setwd(WORKDIR)

# ---- 1. VST-transform (if vsd is not already in memory) ----
if (!exists("vsd")) {
  vsd <- vst(dds, blind = FALSE)
}

# ---- 2. Extract the 22 rescue genes, all 4 conditions ----
rescue_22 <- c(rescued_up$gene, rescued_down$gene)

vsd_mat <- assay(vsd)
heatmap_data <- vsd_mat[rescue_22, ]

# ---- 3. Z-score per gene (row) ----
heatmap_zscore <- t(scale(t(heatmap_data)))

# ---- 4. Order columns: nTg Veh, nTg 915, hTau Veh, hTau 915 ----
col_order <- c(grep("nTg_Veh", colnames(heatmap_zscore), value = TRUE),
               grep("nTg_915", colnames(heatmap_zscore), value = TRUE),
               grep("hTau_Veh", colnames(heatmap_zscore), value = TRUE),
               grep("hTau_915", colnames(heatmap_zscore), value = TRUE))
heatmap_zscore <- heatmap_zscore[, col_order]

# ---- 5. Column annotation: single "Condition" column, CUT&RUN-matched colors ----
condition_label <- case_when(
  grepl("nTg_Veh", col_order) ~ "nTg Veh",
  grepl("nTg_915", col_order) ~ "nTg 915",
  grepl("hTau_Veh", col_order) ~ "hTau Veh",
  grepl("hTau_915", col_order) ~ "hTau 915"
)

annotation_col <- data.frame(
  Condition = factor(condition_label, levels = c("nTg Veh", "nTg 915", "hTau Veh", "hTau 915")),
  row.names = col_order
)

ann_colors <- list(
  Condition = c(
    "nTg Veh" = "#CFCFCF",
    "nTg 915" = "#056AE7",
    "hTau Veh" = "#FFE0C0",
    "hTau 915" = "#FF6A00"
  )
)

# ---- 6. Custom color breaks for increased visual contrast ----
# (does not alter underlying z-score values, only how they are mapped to color)
breaks_custom <- seq(-1.5, 1.5, length.out = 101)

# ---- 7. Plot ----
pheatmap(heatmap_zscore,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = breaks_custom,
         fontsize_row = 9,
         fontsize_col = 8,
         main = "Rescue genes (n=22): all conditions",
         filename = "Heatmap_22_rescue_genes.pdf",
         width = 8, height = 7)

cat("Done! Heatmap_22_rescue_genes.pdf saved in", WORKDIR, "\n")
