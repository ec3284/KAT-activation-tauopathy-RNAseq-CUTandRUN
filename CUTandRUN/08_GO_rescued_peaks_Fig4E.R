# =============================================================================
# GO enrichment of the rescued H3K4ac peaks (Fig. 4e)
#
# Two changes from the first version, and the second is the one that matters.
#
# 1. A universe. enrichGO was called without the `universe` argument, so the
#    background was every gene annotated in org.Mm.eg.db, about 25,000, while
#    the foreground was roughly 12,000 genes that all carry an H3K4ac peak in
#    brain. Against that background a large brain-derived gene list returns
#    neuronal terms whatever the peaks are doing, because half of what is being
#    tested is simply "these genes have active chromatin in brain". The
#    background here is the genes that carry any H3K4ac peak, so the question
#    becomes whether the rescued loci differ from the rest of the marked
#    genome, which is a question the data can answer either way.
#
# 2. No keyword filter. The first version selected terms by name before
#    plotting:
#
#      keywords <- 'synap|memory|learning|neuron|dendrit|potentiation|...'
#      selected <- ego[grepl(keywords, ego$Description, ...), ]
#
#    and then took the top twelve of what survived. The panel therefore showed
#    synaptic terms because synaptic terms were the ones kept, and would have
#    shown them even if the strongest enrichments had been metabolic or
#    immune. Here the terms are ranked by adjusted P value and the top ones are
#    plotted whatever they are. If synaptic terms come out on top, the sentence
#    in the Results is earned rather than imposed, and it is a much stronger
#    sentence for it.
#
# Redundancy is removed the same way as in the RNA-seq figure: GO returns many
# near-copies of one concept, so a term is dropped when it shares most of its
# genes with a term already kept. Semantic simplification is not used, because
# it merges terms by name similarity rather than by what they contain.
#
# Input : the four consensus peak BED files, mm10_genes_sorted.bed
# Output: GO_BP_rescued.csv (all terms, unfiltered), Fig4e_GO_BP_rescued_top.pdf/.png,
#         rescued_genes.txt, universe_genes.txt
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Mm.eg.db)
  library(GenomicRanges); library(rtracklayer); library(ggplot2)
})

# ---- PATHS -------------------------------------------------------------------
# Run this script from the directory where the figure outputs should be written.
# The two locations below are read, never written to. Either edit them or set
# the environment variables CUTRUN_PEAKS and CUTRUN_BIGWIG.
#   CUTRUN_PEAKS   the consensus and per-replicate SEACR peak BED files
#   CUTRUN_BIGWIG  the spike-in normalized bigWig tracks
# GENES_BED is the mm10 gene annotation as BED, made from the GTF; see the README.
PEAKS  <- Sys.getenv("CUTRUN_PEAKS",  "peaks_and_consensus_peaks")
GENES_BED <- Sys.getenv("CUTRUN_GENES_BED", "mm10_genes_sorted.bed")
PEAKS_DIR <- PEAKS
stopifnot(dir.exists(PEAKS_DIR)); stopifnot(file.exists(GENES_BED))

N_SHOW   <- 12     # terms in the dotplot
MAX_OVL  <- 0.5    # drop a term sharing more than this fraction of genes with one already kept

find_one <- function(patterns, label) {
  f <- unique(unlist(lapply(patterns, function(p)
    list.files(PEAKS_DIR, pattern = p, ignore.case = TRUE, full.names = TRUE))))
  if (!length(f)) stop("no file for ", label, " in ", PEAKS_DIR, call. = FALSE)
  f[1]
}
bed <- list(
  nTg_Veh  = find_one(c("^nTg[-_]Veh.*\\.bed$", "nTg-Veh-H3K4ac-consensus"), "nTg Veh"),
  hTau_Veh = find_one(c("^hTau[-_]Veh.*\\.bed$", "hTau-Veh-H3K4ac-consensus"), "hTau Veh"),
  hTau_915 = find_one(c("^hTau[-_](RA013915|915).*\\.bed$", "hTau-RA013915-H3K4ac-consensus"), "hTau 915"),
  nTg_915  = find_one(c("^nTg[-_](RA013915|915).*\\.bed$", "nTg-RA013915-H3K4ac-consensus"), "nTg 915"))
gr <- lapply(bed, import, format = "BED")
for (k in names(gr)) cat(sprintf("%-9s %-46s %s regions\n", k, basename(bed[[k]]),
                                 format(length(gr[[k]]), big.mark = ",")))

# ---- 1. rescued peaks --------------------------------------------------------
# Present in vehicle nTg and in RA013915-treated hTau, absent in vehicle hTau.
# Same definition as the original bash, kept so the set does not change.
lost    <- gr$nTg_Veh[!overlapsAny(gr$nTg_Veh, gr$hTau_Veh)]
rescued <- lost[overlapsAny(lost, gr$hTau_915)]
cat(sprintf("\nlost in vehicle hTau: %s | of these, present in hTau RA013915: %s\n",
            format(length(lost), big.mark = ","), format(length(rescued), big.mark = ",")))
export(rescued, "rescued_peaks.bed", format = "BED")

# ---- 2. genes, and the universe ---------------------------------------------
# A gene is assigned to a peak when the peak overlaps the gene body, the same
# rule as the original awk '$NF==0'.
g <- import(GENES_BED, format = "BED")
genes_of <- function(x) unique(na.omit(g$name[subjectHits(findOverlaps(x, g))]))

rescued_genes <- genes_of(rescued)
union_peaks   <- reduce(sort(do.call(c, unname(gr))))
universe_genes <- genes_of(union_peaks)

writeLines(rescued_genes, "rescued_genes.txt")
writeLines(universe_genes, "universe_genes.txt")
cat(sprintf("rescued genes: %s | universe (genes with any H3K4ac peak): %s (%.1f%%)\n",
            format(length(rescued_genes), big.mark = ","),
            format(length(universe_genes), big.mark = ","),
            100 * length(rescued_genes) / length(universe_genes)))
# If the rescued set covers most of the universe, no enrichment test can
# separate them and the panel should not be built; the percentage above is the
# number to look at before reading anything else.

sym2ent <- function(s) {
  m <- suppressMessages(bitr(s, "SYMBOL", "ENTREZID", OrgDb = org.Mm.eg.db))
  unique(m$ENTREZID)
}
fg <- sym2ent(rescued_genes); bg <- sym2ent(universe_genes)
cat(sprintf("mapped to Entrez: %s of %s foreground, %s of %s background\n",
            format(length(fg), big.mark = ","), format(length(rescued_genes), big.mark = ","),
            format(length(bg), big.mark = ","), format(length(universe_genes), big.mark = ",")))

# ---- 3. enrichment -----------------------------------------------------------
ego <- enrichGO(gene = fg, universe = bg, OrgDb = org.Mm.eg.db, ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                readable = TRUE)
res <- as.data.frame(ego)
write.csv(res, "GO_BP_rescued.csv", row.names = FALSE)
cat(sprintf("\nterms at FDR < 0.05: %s\n", format(sum(res$p.adjust < 0.05), big.mark = ",")))
if (!nrow(res)) stop("no enriched term against this universe", call. = FALSE)

cat("\ntop 20 terms as they come out, with no filtering:\n")
print(head(res[order(res$p.adjust), c("Description", "GeneRatio", "BgRatio", "p.adjust")], 20),
      row.names = FALSE)

# ---- 4. drop near-duplicate terms -------------------------------------------
res <- res[order(res$p.adjust), ]
sets <- strsplit(res$geneID, "/")
keep <- integer(0)
for (i in seq_len(nrow(res))) {
  if (!length(keep)) { keep <- i; next }
  ov <- vapply(keep, function(j)
    length(intersect(sets[[i]], sets[[j]])) / min(length(sets[[i]]), length(sets[[j]])),
    numeric(1))
  if (max(ov) <= MAX_OVL) keep <- c(keep, i)
  if (length(keep) >= N_SHOW) break
}
top <- res[keep, ]
cat(sprintf("\n%d non-redundant terms plotted (overlap <= %.0f%% of the smaller set)\n",
            nrow(top), 100 * MAX_OVL))
print(top[, c("Description", "FoldEnrichment", "Count", "p.adjust")], row.names = FALSE)

# ---- 5. dotplot --------------------------------------------------------------
top$Description <- factor(top$Description, levels = rev(top$Description))
p <- ggplot(top, aes(FoldEnrichment, Description, size = Count, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = "#FFE0C0", high = "#FF6A00", name = expression(-log[10]*"(p.adj)")) +
  scale_size_continuous(name = "Gene count", range = c(4, 12)) +
  labs(x = "Fold enrichment", y = NULL) +
  theme_bw(base_size = 20) +
  theme(text = element_text(family = "Arial"),
        axis.text = element_text(colour = "black", size = 20),
        legend.title = element_text(size = 16), legend.text = element_text(size = 14))

ggsave("Fig4e_GO_BP_rescued_top.pdf", p, width = 13, height = 7.5, device = cairo_pdf)
ggsave("Fig4e_GO_BP_rescued_top.png", p, width = 13, height = 7.5, dpi = 400)
write.csv(top, "Fig4e_GO_BP_rescued_plotted.csv", row.names = FALSE)
cat("\nwritten to", normalizePath("."), "\n")
