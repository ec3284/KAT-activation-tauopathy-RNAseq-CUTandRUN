# ==============================================================================
# Script 06: GO Biological Process enrichment of rescued H3K4ac peaks + dotplot
# Figure: 4G
# Input: rescued_genes_relaxed.txt, generated from rescued_peaks_relaxed.bed
# (see 05_overlap_and_rescued_peaks_Fig4F_G.sh) via closestBed gene annotation
# ==============================================================================

# ---- Gene annotation of rescued peaks (bash step, shown here for reference) ----
# closestBed -a rescued_peaks_relaxed.bed -b mm10_genes_sorted.bed -d \
#     > rescued_relaxed_annotated.bed
# awk '$NF==0 {print $7}' rescued_relaxed_annotated.bed | sort -u \
#     > rescued_genes_relaxed.txt

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

genes <- read.table('rescued_genes_relaxed.txt', header = FALSE)$V1
cat("Input genes:", length(genes), "\n")

entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Mapped genes:", nrow(entrez), "\n")

ego_BP <- enrichGO(
  gene = entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP",
  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
  readable = TRUE
)
write.csv(ego_BP@result, 'GO_BP_rescued_relaxed.csv', row.names = FALSE)
cat("GO BP terms:", nrow(ego_BP@result[ego_BP@result$p.adjust < 0.05, ]), "\n")
cat("Done!\n")

# ---- Dotplot of top rescued-peak GO terms (Fig. 4G) ----
ego <- read.csv('GO_BP_rescued_relaxed.csv')
keywords <- 'synap|memory|learning|neuron|dendrit|potentiation|cognit|plasticity|neurotransmit|glutamat|calcium|hippoc'
selected <- ego[grepl(keywords, ego$Description, ignore.case = TRUE) & ego$p.adjust < 0.05, ]
top12 <- head(selected[order(selected$p.adjust), ], 12)
top12$Description <- factor(top12$Description, levels = rev(top12$Description))

p <- ggplot(top12, aes(x = FoldEnrichment, y = Description, size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = '#FFE0C0', high = '#FF6A00', name = '-log10(p.adj)') +
  scale_size_continuous(name = 'Gene Count', range = c(4, 12)) +
  theme_bw(base_size = 13) +
  labs(x = 'Fold Enrichment', y = '', title = 'GO Biological Process\nRescued H3K4ac Peaks') +
  theme(axis.text.y = element_text(size = 11))

pdf('GO_BP_dotplot_relaxed.pdf', width = 11, height = 7)
print(p)
dev.off()
cat("Saved!\n")
