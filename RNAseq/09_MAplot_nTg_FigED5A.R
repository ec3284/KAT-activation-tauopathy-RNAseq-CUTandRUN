#!/usr/bin/env Rscript
# =============================================================================
# Extended Data Figure 5a - RA013915 versus vehicle in non-transgenic mice
#
# The transcriptional counterpart to the chromatin panels of the same figure:
# what the compound does in a healthy brain, measured on the same animals'
# transcriptome. It is computed from the RNA-seq model and is therefore
# independent of how the CUT&RUN signal is normalised.
#
# Three genes pass the threshold. All three are ribosomal pseudogenes or
# predicted transcripts with nearly identical fold changes and identical
# adjusted P values, the signature of multi-mapping reads distributed across
# paralogues rather than of three genes responding. They are labelled so that
# the figure states this rather than leaving a reader to wonder.
#
# Input : RNAseq_20lib_results.xlsx (from 01_)
# Output: FigED5a_MAplot_nTg.pdf / .png, FigED5a_MAplot_nTg_source.csv
# =============================================================================

rm(list = ls())
library(readxl); library(dplyr); library(ggplot2); library(ggrepel)

# ---- PATHS -------------------------------------------------------------------
OUT <- "."
stopifnot(dir.exists(OUT)); setwd(OUT)
IN     <- "RNAseq_20lib_results.xlsx"
COUNTS <- "GEO_raw_counts_20libraries.tsv"

CONTRAST <- "drug_nTg"   # "drug_nTg", "drug_hTau" or "disease"
FDR <- 0.05

# MODE decides which analysis the panel shows, and the file names follow it so
# that the two cannot be confused.
#
#   "joint"  the model reported throughout the paper: all four groups, all 20
#            libraries, sequencing run in the design, mouse as a blocking
#            factor. Variance is estimated across every group, so the
#            non-transgenic contrast borrows information from the hTau groups.
#   "within" a sensitivity analysis restricted to the eight non-transgenic
#            libraries. These come from a single sequencing run and include no
#            repeated libraries, so no run term and no blocking are needed and
#            no variance is borrowed from another genotype. It has fewer
#            residual degrees of freedom and cannot test an interaction, which
#            is why it is not the primary analysis.
#
# Report whichever is shown as what it is. Using "joint" everywhere except the
# one panel where "within" gives a tidier number would make the choice of model
# depend on the answer.
MODE <- "within"

# ---- STYLE -------------------------------------------------------------------
# #3B6FA0 is the colour used for significant genes in the existing MA plots of
# this manuscript; kept so the panel matches its neighbours.
SIG <- "#3B6FA0"; CLOUD <- "#d9d9d6"; INK <- "#0b0b0b"; MUTED <- "#52514e"
FAMILY <- "Arial"
FONT <- 20
GENE <- 5.3     # gene labels, ggplot size units (points / 2.845)
YL <- 3         # |log2FC| axis limit
FIG_W <- 7.6; FIG_H <- 6.4
LABEL_SIGNIFICANT <- TRUE

if (requireNamespace("extrafont", quietly = TRUE))
  suppressWarnings(try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE))

# ---- 1. data -----------------------------------------------------------------
if (MODE == "joint") {
  stopifnot(file.exists(IN))
  d <- as.data.frame(read_excel(IN, sheet = "Full_Merged_Results"))
  names(d)[1] <- "gene"
  cols <- paste0(c("baseMean_", "log2FoldChange_", "padj_"), CONTRAST)
  stopifnot(all(cols %in% names(d)))
  d <- d %>% transmute(gene,
                       baseMean = .data[[cols[1]]],
                       log2FC   = .data[[cols[2]]],
                       padj     = .data[[cols[3]]])
} else if (MODE == "within") {
  stopifnot(file.exists(COUNTS), CONTRAST == "drug_nTg")
  suppressPackageStartupMessages({ library(readr); library(limma); library(edgeR) })
  x <- as.data.frame(read_tsv(COUNTS, show_col_types = FALSE))
  rownames(x) <- make.unique(ifelse(is.na(x$gene_name) | x$gene_name == "",
                                    x$gene_id, x$gene_name))
  cn <- setdiff(colnames(x), c("gene_id", "gene_name"))
  veh  <- sort(grep("^nTg_Veh_[0-9]+$", cn, value = TRUE))
  drug <- sort(grep("^nTg_RA013915_[0-9]+$", cn, value = TRUE))
  stopifnot(length(veh) == 4, length(drug) == 4)
  cat("libraries used:", paste(c(veh, drug), collapse = ", "), "\n")

  m <- round(as.matrix(x[, c(veh, drug)]))
  grp <- factor(rep(c("Veh", "RA013915"), each = 4), levels = c("Veh", "RA013915"))
  des <- model.matrix(~ grp)
  y <- DGEList(m, group = grp)
  y <- calcNormFactors(y[filterByExpr(y, design = des), , keep.lib.sizes = FALSE])
  cat(sprintf("genes retained: %d | residual d.f.: %d\n", nrow(y), ncol(m) - ncol(des)))
  fit <- eBayes(lmFit(voom(y, des, plot = FALSE), des))
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  d <- data.frame(gene = rownames(tt), baseMean = 2^tt$AveExpr,
                  log2FC = tt$logFC, padj = tt$adj.P.Val)
} else stop("MODE must be 'joint' or 'within'")

d <- d %>% filter(!is.na(baseMean), baseMean > 0, !is.na(log2FC)) %>%
  mutate(sig = !is.na(padj) & padj < FDR)

n_sig <- sum(d$sig)
cat(sprintf("\nMODE = %s | contrast %s: %d genes plotted, %d significant at FDR < %g\n",
            MODE, CONTRAST, nrow(d), n_sig, FDR))
if (n_sig) {
  cat("significant genes:\n")
  print(d %>% filter(sig) %>% arrange(padj) %>%
          transmute(gene, baseMean = round(baseMean, 1),
                    log2FC = round(log2FC, 3), padj = signif(padj, 3)),
        row.names = FALSE)
}

# Genes beyond the axis limits are drawn as triangles on the border rather than
# dropped, the same convention as the other MA plots in this manuscript.
d <- d %>% mutate(y = pmax(pmin(log2FC, YL), -YL),
                  clipped = abs(log2FC) > YL,
                  lab = if (LABEL_SIGNIFICANT) ifelse(sig, gene, NA_character_) else NA_character_)
cat(sprintf("genes drawn on the border: %d\n", sum(d$clipped)))

# ---- 2. panel ----------------------------------------------------------------
p <- ggplot(d, aes(baseMean, y)) +
  geom_hline(yintercept = 0, colour = MUTED, linewidth = 0.6) +
  geom_point(data = . %>% filter(!sig, !clipped),
             colour = CLOUD, size = 0.5, alpha = 0.5) +
  geom_point(data = . %>% filter(!sig, clipped),
             aes(shape = log2FC > 0), colour = CLOUD, size = 1.1, alpha = 0.7) +
  geom_point(data = . %>% filter(sig, !clipped),
             colour = SIG, size = 2.2) +
  geom_point(data = . %>% filter(sig, clipped),
             aes(shape = log2FC > 0), colour = SIG, fill = SIG, size = 2.6) +
  scale_shape_manual(values = c(`TRUE` = 24, `FALSE` = 25), guide = "none") +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_y_continuous(breaks = seq(-YL, YL, 1)) +
  coord_cartesian(ylim = c(-YL, YL), expand = TRUE, clip = "off") +
  # The x value is 2^(average log2-CPM), which is an average expression level,
  # not the "mean of normalized counts" that DESeq2 plots. Labelling it
  # correctly matters here because this panel sits beside MA plots produced
  # earlier with DESeq2.
  labs(x = "Average expression (CPM)",
       y = expression(log[2]*" fold change")) +
  theme_classic(base_size = FONT, base_family = FAMILY) +
  theme(text = element_text(colour = INK, family = FAMILY),
        axis.text = element_text(colour = INK, size = FONT),
        axis.title = element_text(size = FONT),
        legend.position = "none")

if (LABEL_SIGNIFICANT && n_sig)
  p <- p + ggrepel::geom_text_repel(
    data = d %>% filter(sig), aes(label = lab),
    family = FAMILY, size = GENE, colour = INK, fontface = "italic",
    min.segment.length = 0, segment.colour = MUTED, segment.size = 0.3,
    box.padding = 0.5, max.overlaps = Inf, seed = 1)

stem <- paste0("FigED5a_MAplot_nTg_", MODE)
tryCatch(ggsave(paste0(stem, ".pdf"), p, width = FIG_W, height = FIG_H, device = cairo_pdf),
  error = function(e) { message("cairo_pdf failed: ", conditionMessage(e))
    ggsave(paste0(stem, ".pdf"), p, width = FIG_W, height = FIG_H) })
ggsave(paste0(stem, ".png"), p, width = FIG_W, height = FIG_H, dpi = 400)

# ---- 3. source data ----------------------------------------------------------
write.csv(d %>% dplyr::select(gene, baseMean, log2FC, padj, sig) %>% arrange(padj),
          paste0(stem, "_source.csv"), row.names = FALSE)

cat("\nwritten to", normalizePath(OUT), ":\n ",
    paste0(stem, c(".pdf", ".png", "_source.csv"), collapse = "\n  "), "\n")
cat("\nRun this script twice, once with MODE = \"joint\" and once with\n",
    "MODE = \"within\", and state in the legend which analysis each panel shows.\n", sep = "")
