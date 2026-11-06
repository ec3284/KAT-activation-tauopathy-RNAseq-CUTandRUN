#!/usr/bin/env Rscript
# =============================================================================
# Figure 3d - Expression of the rescued genes across the four groups
#
# This panel carries the word "restores". The scatter panels show fold changes;
# only this one shows where each animal actually sits, so that the return of
# hTau RA013915 towards the nTg vehicle level is visible rather than inferred.
#
# THE BATCH PROBLEM, AND WHY IT CAN BE SOLVED HERE
# nTg vehicle and nTg RA013915 were sequenced in run C3SVT8, hTau RA013915 in
# run VJ385X. Placing the four groups side by side therefore compares across
# runs exactly where it matters most. The four hTau vehicle animals were
# sequenced in BOTH runs, so for every gene the difference between their two
# libraries is the run effect measured on paired samples. removeBatchEffect()
# uses that to put all samples on one scale, and the two libraries of each
# hTau vehicle animal are then averaged so that every column is one mouse.
#
# The correction is applied for display only. Every statistic reported in the
# manuscript comes from the model, which handles the run effect as a term in
# the design and never from these adjusted values.
#
# Input : the two expression matrices, or the combined 20-library count matrix
# Output: Fig3d_heatmap_rescued.pdf / .png, Fig3d_group_means.csv
# =============================================================================

rm(list = ls())
library(readr); library(DESeq2); library(limma); library(pheatmap)

# ---- PATHS -------------------------------------------------------------------
COUNTS <- "GEO_raw_counts_20libraries.tsv"   # gene_id, gene_name, 20 count columns
OUT    <- "."
stopifnot(file.exists(COUNTS), dir.exists(OUT)); setwd(OUT)

# ---- GENE SETS ---------------------------------------------------------------
RESCUED      <- c("Apold1", "Btg2", "Dusp1", "Fos", "Ier2", "Klf2")
SUBTHRESHOLD <- c("Ccn1", "Cldn5", "Klf4")   # same direction, above FDR threshold
GENES <- c(RESCUED, SUBTHRESHOLD)
# Both sets are shown, separated by a row annotation rather than by omission:
# dropping the three would hide that the boundary is a statistical threshold.
# Set GENES <- RESCUED for a six-row panel.

# ---- STYLE -------------------------------------------------------------------
GROUP_COLS <- c(nTg_Veh       = "#CFCFCF",
                nTg_RA013915  = "#056AE7",
                hTau_Veh      = "#FFE0C0",
                hTau_RA013915 = "#FF6A00")
# A hue used nowhere else in the figure. Reusing the group oranges here would
# make one colour carry two meanings: "treated group" along the top and
# "passes the threshold" down the side.
# Short labels on purpose: a long legend entry sets the width of the whole
# figure. What they mean belongs in the figure legend, not in the panel.
ROW_COLS <- list(Evidence = c(`FDR < 0.05` = "#5B4E8C",
                              `n.s.`       = "#DCD6EA"))
# Red-blue for expression; the four manuscript colours are reserved for the
# group bar, so that each colour carries exactly one meaning in the figure.
PALETTE <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)
CLIP <- NA     # NA = scale set by the data. A low value saturates hTau vehicle
               # into a uniform block and hides between-animal variation.
FAMILY   <- "Arial"
FONT_ROW <- 18 # gene names; matched to the other panels of the figure rather
               # than to the cell height, so that the assembled figure has one
               # type size throughout
FONT     <- 18 # legend titles, legend entries and the numbers on the colour bar
               # Raising this widens the legend block, so LEGEND_IN below has to
               # grow with it or the text is clipped.
# Cell geometry in points, which is what sets the shape of the panel. Lower
# CELL_W and higher CELL_H for narrow, tall cells; the figure size is computed
# from these below, so the panel keeps its proportions whatever you choose.
CELL_W <- 22    # width of one sample column
CELL_H <- 42    # height of one gene row
# Space to reserve beside the cells, in inches: gene names on the right, then
# the two legends. Increase LABEL_IN if long gene names are clipped.
LABEL_IN <- 1.6     # gene names on the right
LEGEND_IN <- 3.3    # the two legends; scales with FONT
LEFT_IN <- 0.5      # the Evidence bar down the left edge
MARGIN_IN <- 1.2    # group bar on top plus margins
GAP_BETWEEN_BLOCKS <- FALSE   # TRUE inserts a gap between the two evidence blocks

# ---- 1. counts and design ----------------------------------------------------
d <- as.data.frame(read_tsv(COUNTS, show_col_types = FALSE))
rownames(d) <- make.unique(ifelse(is.na(d$gene_name) | d$gene_name == "",
                                  d$gene_id, d$gene_name))

# Columns are found by pattern rather than by exact name. The count matrix has
# been written with the run suffix in two different positions at different
# times (hTau_Veh_1_runA and hTau_Veh_runA_1), and a hardcoded list fails on
# whichever one it was not written for, with an error that points at the
# symptom instead of the cause.
cn <- setdiff(colnames(d), c("gene_id", "gene_name"))
pick <- function(pat) sort(grep(pat, cn, value = TRUE))
LIB <- list(
  nTg_Veh       = pick("^nTg_Veh_[0-9]+$"),
  nTg_RA013915  = pick("^nTg_RA013915_[0-9]+$"),
  hTau_Veh_runA = pick("^hTau_Veh.*runA"),
  hTau_Veh_runB = pick("^hTau_Veh.*runB"),
  hTau_RA013915 = pick("^hTau_RA013915_[0-9]+$"))
cat("libraries found:\n")
for (k in names(LIB)) cat(sprintf("  %-14s %d: %s\n", k, length(LIB[[k]]),
                                  paste(LIB[[k]], collapse = ", ")))
if (!all(lengths(LIB) == 4))
  stop("expected four libraries in each group; found ",
       paste(sprintf("%s=%d", names(LIB), lengths(LIB)), collapse = ", "))

LIBS <- unlist(LIB, use.names = FALSE)
counts <- round(as.matrix(d[, LIBS]))

info <- data.frame(
  row.names = LIBS,
  group = factor(c(rep("nTg_Veh", 4), rep("nTg_RA013915", 4),
                   rep("hTau_Veh", 8), rep("hTau_RA013915", 4)),
                 levels = names(GROUP_COLS)),
  batch = factor(c(rep("C3SVT8", 12), rep("VJ385X", 8)), levels = c("C3SVT8", "VJ385X")))
cat(sprintf("matrix: %d genes x %d libraries\n", nrow(counts), ncol(counts)))
print(table(info$group, info$batch))

# ---- 2. variance-stabilizing transformation ----------------------------------
dds <- DESeqDataSetFromMatrix(counts, info, design = ~ group)
dds <- dds[rowSums(counts(dds)) >= 10, ]
vst <- assay(vst(dds, blind = TRUE))
cat(sprintf("genes after filtering: %d\n", nrow(vst)))

# ---- 3. run effect, estimated from the paired hTau vehicle animals -----------
runA <- LIB$hTau_Veh_runA; runB <- LIB$hTau_Veh_runB
offset <- rowMeans(vst[, runB]) - rowMeans(vst[, runA])
cat(sprintf("per-gene run effect: median %+0.3f | IQR %.3f | |offset| > 1 for %d genes\n",
            median(offset), IQR(offset), sum(abs(offset) > 1)))

vst_adj <- removeBatchEffect(vst, batch = info$batch,
                             design = model.matrix(~ group, info))

# one column per animal: the two libraries of each hTau vehicle mouse are averaged
m <- cbind(vst_adj[, c(LIB$nTg_Veh, LIB$nTg_RA013915)],
           (vst_adj[, runA] + vst_adj[, runB]) / 2,
           vst_adj[, LIB$hTau_RA013915])
colnames(m)[9:12] <- paste0("hTau_Veh_", 1:4)
grp <- factor(sub("_[0-9]+$", "", colnames(m)), levels = names(GROUP_COLS))
cat(sprintf("after averaging the paired runs: %d columns, one per animal\n", ncol(m)))

# ---- 4. matrix for the panel -------------------------------------------------
sel <- intersect(GENES, rownames(m))
if (length(sel) < length(GENES))
  warning("not found: ", paste(setdiff(GENES, sel), collapse = ", "), call. = FALSE)
M   <- m[sel, order(grp), drop = FALSE]
grp <- grp[order(grp)]

# Rows are ordered by degree of recovery rather than by clustering: with nine
# genes a dendrogram imposes an order with no interpretation, whereas this makes
# the ordering itself part of the message.
#
# Evidence takes precedence over recovery in the ordering. Sorting by recovery
# alone interleaves the two categories, and the annotation bar then reads as an
# alternating pattern instead of showing which genes meet the threshold. Within
# each block the order is still by recovery.
recovery <- rowMeans(M[, grp == "hTau_RA013915", drop = FALSE]) -
            rowMeans(M[, grp == "nTg_Veh",       drop = FALSE])
rescued_row <- rownames(M) %in% RESCUED
M <- M[order(!rescued_row, -recovery), , drop = FALSE]
cat("\nrecovery per gene (hTau RA013915 minus nTg vehicle, VST units):\n")
print(data.frame(gene = rownames(M),
                 recovery = round(recovery[rownames(M)], 3),
                 evidence = ifelse(rownames(M) %in% RESCUED,
                                   "FDR < 0.05 in both", "same direction, NS")),
      row.names = FALSE)
gaps_row <- if (GAP_BETWEEN_BLOCKS) sum(rownames(M) %in% RESCUED) else NULL

col_ann <- data.frame(Group = grp, row.names = colnames(M))
# Categories are taken from the names of ROW_COLS rather than written out
# again here: typing them twice lets the two drift apart, and pheatmap then
# silently draws empty cells instead of raising an error.
EV <- names(ROW_COLS$Evidence)   # [1] meets the threshold, [2] does not
row_ann <- data.frame(
  Evidence = factor(ifelse(rownames(M) %in% RESCUED, EV[1], EV[2]), levels = EV),
  row.names = rownames(M))
stopifnot(!anyNA(row_ann$Evidence))

z <- t(scale(t(M)))
cat(sprintf("\nrow z-scores: range %.2f to %.2f\n", min(z), max(z)))
cat("mean z-score per group:\n")
print(round(t(apply(z, 1, function(x) tapply(x, grp, mean))), 2))

# ---- 5. draw -----------------------------------------------------------------
FIG_W <- CELL_W * ncol(M) / 72 + LEFT_IN + LABEL_IN + LEGEND_IN
FIG_H <- CELL_H * nrow(M) / 72 + MARGIN_IN
cat(sprintf("\ncells %g x %g pt -> figure %.2f x %.2f inches\n",
            CELL_W, CELL_H, FIG_W, FIG_H))
cat("  If anything is clipped, raise LEGEND_IN (legend text on the right),\n",
    "  LABEL_IN (gene names) or LEFT_IN (the Evidence bar).\n", sep = "")

draw <- function()
  pheatmap(M, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
           show_colnames = FALSE, cellwidth = CELL_W, cellheight = CELL_H,
           annotation_col = col_ann, annotation_row = row_ann,
           annotation_colors = c(list(Group = GROUP_COLS), ROW_COLS),
           annotation_names_col = FALSE, annotation_names_row = FALSE,
           gaps_col = cumsum(table(grp)), gaps_row = gaps_row,
           color = PALETTE,
           breaks = if (is.na(CLIP)) NA else seq(-CLIP, CLIP, length.out = 102),
           fontsize = FONT, fontsize_row = FONT_ROW, border_color = "grey92")

png("Fig3d_heatmap_rescued.png", width = FIG_W * 200, height = FIG_H * 200, res = 200)
draw(); dev.off()

tryCatch({
  cairo_pdf("Fig3d_heatmap_rescued.pdf", width = FIG_W, height = FIG_H, family = FAMILY)
  draw(); dev.off()
}, error = function(e) {
  message("cairo_pdf failed: ", conditionMessage(e), " - retrying without the font")
  try(dev.off(), silent = TRUE)
  pdf("Fig3d_heatmap_rescued.pdf", width = FIG_W, height = FIG_H); draw(); dev.off()
})
cat(sprintf("PDF written: %s (%.0f KB)\n", file.exists("Fig3d_heatmap_rescued.pdf"),
            file.size("Fig3d_heatmap_rescued.pdf") / 1024))

# ---- 6. source data ----------------------------------------------------------
# The source data for a heatmap has to be the value behind every cell, one
# column per animal, not the group means: a reader cannot reconstruct the panel
# from four numbers per gene. Both the plotted value and the row z-score are
# given, since the panel is drawn with scale = "row".
z_out <- t(scale(t(M)))
src <- data.frame(
  gene = rownames(M),
  evidence = as.character(row_ann[rownames(M), "Evidence"]),
  round(M, 4), check.names = FALSE)
names(src)[-(1:2)] <- paste0(names(src)[-(1:2)], "_VST")
src_z <- data.frame(round(z_out, 4), check.names = FALSE)
names(src_z) <- paste0(names(src_z), "_zscore")
src <- cbind(src, src_z)
write.csv(src, "Fig3d_heatmap_source.csv", row.names = FALSE)

means <- t(apply(M, 1, function(x) tapply(x, grp, mean)))
print(round(means, 2))
write.csv(data.frame(gene = rownames(means), round(means, 4),
                     evidence = row_ann[rownames(means), "Evidence"],
                     check.names = FALSE),
          "Fig3d_group_means.csv", row.names = FALSE)

cat("\nwritten to", normalizePath(OUT), ":\n",
    " Fig3d_heatmap_rescued.pdf\n Fig3d_heatmap_rescued.png\n",
    " Fig3d_heatmap_source.csv\n Fig3d_group_means.csv\n")
cat("Values are batch-adjusted VST, for display only; all reported statistics\n",
    "come from the model, where the run effect is a term in the design.\n", sep = "")
