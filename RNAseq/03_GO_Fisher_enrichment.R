#!/usr/bin/env Rscript
# ==============================================================================
# Script 03: GO enrichment (disease signature) and Fisher's exact test
# enrichment of disease-associated pathways among RA013915-responsive genes
#
# Requires objects from 01_DESeq2_collapsed_rescue_analysis.R to be in memory:
#   - df_merged (merged disease + treatment DESeq2 results)
#   - rescued_up, rescued_down (validated rescue gene lists)
#
# Rationale: beyond the small set of "rescue" genes (significant AND reversed
# in both comparisons), a much larger set of genes is significantly altered
# by treatment alone (no cross-filter). This script tests whether that larger
# set is enriched, at the pathway level, for genes independently found to be
# altered in the disease signature — i.e., whether RA013915 broadly modulates
# disease-relevant biological programs even without exact gene-level reversal.
# ==============================================================================

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)

WORKDIR <- Sys.getenv("WORKDIR", unset = "./rnaseq_work")
setwd(WORKDIR)
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)

# ---- 1. GO Biological Process enrichment on the disease signature ----
disease_genes_list <- df_merged %>% filter(padj_disease < 0.05) %>% pull(gene)
entrez_disease <- bitr(disease_genes_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

go_disease <- enrichGO(
  gene = entrez_disease$ENTREZID, OrgDb = org.Mm.eg.db, keyType = "ENTREZID",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE
)

cat("Disease signature GO terms found:", nrow(go_disease@result), "\n")

# ---- 2. Fisher's exact test: for each top disease GO category, is it
#         enriched among treatment-responsive genes (padj<0.05, no LFC
#         cutoff), relative to the genome-wide background? ----
treatment_sig_0 <- df_merged %>% filter(padj_treatment < 0.05)
total_genes_tested <- nrow(df_merged)
total_treatment_sig <- nrow(treatment_sig_0)

categories_to_test <- head(go_disease@result$Description, 30)
fisher_results_full <- data.frame()

for (cat_desc in categories_to_test) {
  gene_str <- go_disease@result$geneID[go_disease@result$Description == cat_desc]
  gene_list <- strsplit(gene_str, "/")[[1]]
  in_treat <- intersect(gene_list, treatment_sig_0$gene)

  m <- matrix(c(length(in_treat), length(gene_list) - length(in_treat),
                total_treatment_sig - length(in_treat),
                total_genes_tested - length(gene_list) - (total_treatment_sig - length(in_treat))),
              nrow = 2)
  ft <- fisher.test(m)

  fisher_results_full <- rbind(fisher_results_full, data.frame(
    Category = cat_desc, N_overlap = length(in_treat), N_total_disease = length(gene_list),
    Pct = round(length(in_treat) / length(gene_list) * 100, 1),
    OddsRatio = unname(ft$estimate), Pvalue = ft$p.value
  ))
}

fisher_results_full$padj <- p.adjust(fisher_results_full$Pvalue, method = "BH")
fisher_results_full <- fisher_results_full %>% arrange(Pvalue)

cat("\nFisher enrichment results (top disease GO categories vs treatment-responsive genes):\n")
print(fisher_results_full, row.names = FALSE)

# ---- 3. Define validated pathway gene sets (padj<0.05 in the Fisher test) ----
# These groupings were confirmed significant and are used for volcano plot
# coloring (see 04_volcano_plot.R). Update the Description strings here if
# go_disease results change with a re-run.
mito_all <- unique(unlist(strsplit(go_disease@result$geneID[go_disease@result$Description %in%
  c("oxidative phosphorylation", "mitochondrion organization", "cellular respiration")], "/")))

neuronal_all <- unique(unlist(strsplit(go_disease@result$geneID[go_disease@result$Description %in%
  c("forebrain development", "dendrite development", "positive regulation of cell projection organization")], "/")))

synaptic_all <- unique(unlist(strsplit(go_disease@result$geneID[go_disease@result$Description %in%
  c("regulation of synapse organization", "regulation of synapse structure or activity")], "/")))

cat("\nValidated pathway gene set sizes:\n",
    "Mitochondrial:", length(mito_all), "\n",
    "Neuronal:", length(neuronal_all), "\n",
    "Synaptic:", length(synaptic_all), "\n")

# ---- 4. Save everything ----
save(go_disease, fisher_results_full, mito_all, neuronal_all, synaptic_all,
     file = "GO_Fisher_enrichment_results.RData")

wb <- createWorkbook()
addWorksheet(wb, "Disease_GO_BP")
writeData(wb, "Disease_GO_BP", go_disease@result)
addWorksheet(wb, "Fisher_Enrichment_Treatment")
writeData(wb, "Fisher_Enrichment_Treatment", fisher_results_full)
saveWorkbook(wb, "RNAseq_GO_Fisher_enrichment.xlsx", overwrite = TRUE)

# ---- 5. Fisher enrichment dotplot (Figure 3B) ----
library(ggplot2)

fisher_plot_data <- fisher_results_full %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(15)

p_fisher_dotplot <- ggplot(fisher_plot_data, aes(x = OddsRatio, y = reorder(Category, OddsRatio))) +
  geom_point(aes(size = N_overlap, color = -log10(padj))) +
  scale_color_gradient(low = "#FFB380", high = "#CC3300", name = "-log10(padj)") +
  scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
  theme_bw(base_size = 12) +
  labs(x = "Odds Ratio", y = NULL,
       title = "Disease-associated pathways enriched\namong RA013915-responsive genes")

print(p_fisher_dotplot)
ggsave("GO_Fisher_dotplot_FINAL.pdf", p_fisher_dotplot, width = 8, height = 6, dpi = 300)

cat("\nDone! Saved in", WORKDIR, ":\n",
    "- GO_Fisher_enrichment_results.RData (load this to skip re-running GO/Fisher)\n",
    "- RNAseq_GO_Fisher_enrichment.xlsx\n",
    "- GO_Fisher_dotplot_FINAL.pdf (Figure 3B)\n",
    "\nRun 04_volcano_plot.R next to generate the categorized volcano plot.\n")
