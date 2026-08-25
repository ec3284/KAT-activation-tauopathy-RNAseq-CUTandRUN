#!/usr/bin/env Rscript
# ==============================================================================
# Script 04: Categorized volcano plot (hTau: RA013915 vs Vehicle)
#
# Requires objects in memory from:
#   - 01_DESeq2_collapsed_rescue_analysis.R (df_merged, rescued_up, rescued_down)
#   - 03_GO_Fisher_enrichment.R (mito_all, neuronal_all, synaptic_all)
#     -- or load("GO_Fisher_enrichment_results.RData") to restore them
#
# Points are colored by validated pathway category (Fisher padj<0.05 vs
# disease signature; see script 03). Only the 22 cross-filtered "rescue"
# genes are labeled by name, in black for readability against the
# colored/grey background points. Predicted/uncharacterized genes
# (Gm..., Mir...) are excluded from labels for readability.
# ==============================================================================

library(ggplot2)
library(ggrepel)
library(dplyr)

WORKDIR <- Sys.getenv("WORKDIR", unset = "./rnaseq_work")
setwd(WORKDIR)

# If starting fresh and only script 01 was re-run, load the saved GO/Fisher results:
# load("GO_Fisher_enrichment_results.RData")

# ---- 1. Build the annotated data frame ----
res_hTau_df <- df_merged %>%
  rename(gene_name = gene, padj = padj_treatment, log2FoldChange = log2FoldChange_treatment)

ieg_rescue <- c(rescued_up$gene, rescued_down$gene)

res_hTau_df$category <- "n.s"
res_hTau_df$category[res_hTau_df$padj < 0.05 & res_hTau_df$log2FoldChange < 0] <- "Down"
res_hTau_df$category[res_hTau_df$padj < 0.05 & res_hTau_df$log2FoldChange > 0] <- "Other_up"
res_hTau_df$category[res_hTau_df$gene_name %in% synaptic_all & res_hTau_df$padj < 0.05] <- "Synaptic"
res_hTau_df$category[res_hTau_df$gene_name %in% neuronal_all & res_hTau_df$padj < 0.05] <- "Neuronal"
res_hTau_df$category[res_hTau_df$gene_name %in% mito_all & res_hTau_df$padj < 0.05] <- "Mitochondrial"
res_hTau_df$category[res_hTau_df$gene_name %in% ieg_rescue & res_hTau_df$padj < 0.05] <- "Rescue_IEG"

res_hTau_df$category <- factor(res_hTau_df$category,
                                levels = c("Mitochondrial", "Neuronal", "Synaptic", "Rescue_IEG",
                                           "Other_up", "Down", "n.s"))

n_down <- sum(res_hTau_df$padj < 0.05 & res_hTau_df$log2FoldChange < 0, na.rm = TRUE)
n_up   <- sum(res_hTau_df$padj < 0.05 & res_hTau_df$log2FoldChange > 0, na.rm = TRUE)

cat("Category counts (padj<0.05):\n")
print(table(res_hTau_df$category[res_hTau_df$padj < 0.05]))

# ---- 2. Colors ----
cat_colors <- c("Mitochondrial" = "#1B7837", "Neuronal" = "#4FC3F7", "Synaptic" = "#1565C0",
                 "Rescue_IEG" = "#FF6A00", "Other_up" = "grey60", "Down" = "#4D4D4D", "n.s" = "grey88")

# ---- 3. Labels: only the 22 rescue genes, excluding predicted/uncharacterized genes ----
top_label <- res_hTau_df[res_hTau_df$category == "Rescue_IEG" & !is.na(res_hTau_df$padj) &
                            !grepl("^Gm[0-9]|^Mir", res_hTau_df$gene_name), ]

# ---- 4. Plot ----
p_volcano <- ggplot() +
  geom_point(data = subset(res_hTau_df, category == "n.s"),
             aes(x = log2FoldChange, y = -log10(padj), color = category), size = 0.5, alpha = 0.3) +
  geom_point(data = subset(res_hTau_df, category == "Down"),
             aes(x = log2FoldChange, y = -log10(padj), color = category), size = 1.2, alpha = 0.6) +
  geom_point(data = subset(res_hTau_df, category == "Other_up"),
             aes(x = log2FoldChange, y = -log10(padj), color = category), size = 1.2, alpha = 0.5) +
  geom_point(data = subset(res_hTau_df, category == "Synaptic"),
             aes(x = log2FoldChange, y = -log10(padj), color = category), size = 3, alpha = 0.9) +
  geom_point(data = subset(res_hTau_df, category == "Neuronal"),
             aes(x = log2FoldChange, y = -log10(padj), color = category), size = 3, alpha = 0.9) +
  geom_point(data = subset(res_hTau_df, category == "Mitochondrial"),
             aes(x = log2FoldChange, y = -log10(padj), color = category), size = 3, alpha = 0.9) +
  geom_point(data = subset(res_hTau_df, category == "Rescue_IEG"),
             aes(x = log2FoldChange, y = -log10(padj), color = category), size = 4, alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 2, color = "grey50", linewidth = 0.4) +
  geom_text_repel(data = top_label,
                   aes(x = log2FoldChange, y = -log10(padj), label = gene_name),
                   color = "black", size = 5, fontface = "bold",
                   max.overlaps = 30, box.padding = 0.6, point.padding = 0.4, force = 2) +
  scale_color_manual(values = cat_colors, name = "Category",
                      breaks = c("Mitochondrial", "Neuronal", "Synaptic", "Rescue_IEG"),
                      labels = c("Mitochondrial", "Neuronal", "Synaptic", "Rescue (IEG/vascular)")) +
  theme_bw(base_size = 13) +
  theme(text = element_text(family = "Arial")) +
  labs(x = "log2 FC (RA915/Veh)", y = "-log10 adjusted p-value",
       title = paste0(n_down, " sig down / ", n_up, " sig up")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

print(p_volcano)
ggsave("Volcano_hTau_FINAL.pdf", p_volcano, width = 10, height = 7, dpi = 300, device = cairo_pdf)

cat("\nDone! Saved Volcano_hTau_FINAL.pdf in", WORKDIR, "\n")
