#!/usr/bin/env Rscript
# =============================================================================
# Supplementary Figure 1 - sequencing run effect and its removal
#
# Two things this figure has to establish, because the whole analysis rests on
# them.
#
# (a) The run effect is real and is removable. Panel A shows the samples before
#     adjustment, where the first component separates the two sequencing runs
#     as much as it separates the genotypes, and panel B after adjustment. The
#     test of a correct adjustment is not that the groups look tidy: it is that
#     the two libraries of each hTau vehicle animal, which are the same RNA,
#     come to lie on top of each other. They are drawn joined by a segment so
#     that this can be checked by eye rather than taken on trust.
#
# (b) Those repeated libraries are concordant. Panel C compares each animal
#     with itself across the two runs, against the line of identity. The
#     libraries agree closely, and the residual displacement from that line is
#     the run effect itself - which is what the term in the design removes.
#
# The justification for blocking on animal does not come from this figure. It
# comes from the consensus within-animal correlation estimated by
# duplicateCorrelation inside the model (0.62), which is net of the run effect
# and is reported in the methods. A naive comparison of within-animal against
# between-animal correlations would be confounded, because only the first
# contains a difference of run.
#
# Genes are filtered with the full design, as in 01_, so that this figure is
# computed on the same gene set as the differential expression analysis.
#
# Values shown here are adjusted for display. Every statistic in the manuscript
# comes from the model, where the run effect is a term in the design.
#
# Input : GEO_raw_counts_20libraries.tsv  (written by 01_)
# Output: FigS1_QC_batch_technical.pdf / .png, FigS1_QC_correlations.csv
# =============================================================================

rm(list = ls())
library(readr); library(limma); library(edgeR); library(ggplot2); library(dplyr)

# ---- PATHS -------------------------------------------------------------------
COUNTS <- "GEO_raw_counts_20libraries.tsv"
OUT    <- "."
stopifnot(file.exists(COUNTS), dir.exists(OUT)); setwd(OUT)

# ---- STYLE -------------------------------------------------------------------
GROUP_COLS <- c(nTg_Veh       = "#CFCFCF",
                nTg_RA013915  = "#056AE7",
                hTau_Veh      = "#FFE0C0",
                hTau_RA013915 = "#FF6A00")
INK <- "#0b0b0b"; MUTED <- "#52514e"
FAMILY <- "Arial"
FONT <- 20
N_TOP <- 2000    # most variable genes used for the ordination
FIG_W <- 15; FIG_H <- 9

if (requireNamespace("extrafont", quietly = TRUE))
  suppressWarnings(try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE))

base_theme <- theme_classic(base_size = FONT, base_family = FAMILY) +
  theme(text = element_text(colour = INK, family = FAMILY),
        axis.text = element_text(colour = INK, size = FONT),
        axis.title = element_text(size = FONT),
        plot.title = element_text(size = FONT, face = "plain", hjust = 0),
        legend.text = element_text(size = FONT - 4),
        legend.title = element_blank())

# ---- 1. data -----------------------------------------------------------------
d <- as.data.frame(read_tsv(COUNTS, show_col_types = FALSE))
rownames(d) <- make.unique(ifelse(is.na(d$gene_name) | d$gene_name == "",
                                  d$gene_id, d$gene_name))

# found by pattern, so that either spelling of the run suffix works
cn <- setdiff(colnames(d), c("gene_id", "gene_name"))
pick <- function(pat) sort(grep(pat, cn, value = TRUE))
LIB <- list(
  nTg_Veh       = pick("^nTg_Veh_[0-9]+$"),
  nTg_RA013915  = pick("^nTg_RA013915_[0-9]+$"),
  hTau_Veh_runA = pick("^hTau_Veh.*runA"),
  hTau_Veh_runB = pick("^hTau_Veh.*runB"),
  hTau_RA013915 = pick("^hTau_RA013915_[0-9]+$"))
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
  run = factor(c(rep("C3SVT8", 12), rep("VJ385X", 8)), levels = c("C3SVT8", "VJ385X")),
  animal = c(paste0("nTgV", 1:4), paste0("nTgD", 1:4),
             paste0("hTauV", 1:4), paste0("hTauV", 1:4), paste0("hTauD", 1:4)))

# the same design used in 01_, so that the gene set matches the one the
# differential expression model was fitted on
design <- model.matrix(~ run + group, info)

dge <- DGEList(counts, group = info$group)
dge <- calcNormFactors(dge[filterByExpr(dge, design = design), , keep.lib.sizes = FALSE])
lcpm <- cpm(dge, log = TRUE, prior.count = 2)
cat(sprintf("%d genes retained across %d libraries\n", nrow(lcpm), ncol(lcpm)))

lcpm_adj <- removeBatchEffect(lcpm, batch = info$run,
                              design = model.matrix(~ group, info))

RUN_A <- LIB$hTau_Veh_runA
RUN_B <- LIB$hTau_Veh_runB

# ---- 2. panels A and B: ordination before and after adjustment ---------------
ordination <- function(M, title) {
  v <- apply(M, 1, var)
  p <- prcomp(t(M[order(-v)[seq_len(min(N_TOP, nrow(M)))], ]))
  pct <- 100 * p$sdev^2 / sum(p$sdev^2)
  df <- data.frame(PC1 = p$x[, 1], PC2 = p$x[, 2], info[rownames(p$x), ],
                   library = rownames(p$x))
  # segments joining the two libraries of each repeated animal
  pairs <- data.frame(x = df[RUN_A, "PC1"], y = df[RUN_A, "PC2"],
                      xend = df[RUN_B, "PC1"], yend = df[RUN_B, "PC2"])
  gap <- sqrt((pairs$x - pairs$xend)^2 + (pairs$y - pairs$yend)^2)
  cat(sprintf("%-28s PC1 %.1f%% PC2 %.1f%% | median distance between paired libraries %.2f\n",
              title, pct[1], pct[2], median(gap)))
  g <- ggplot(df, aes(PC1, PC2)) +
    geom_segment(data = pairs, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = MUTED, linewidth = 0.4) +
    geom_point(aes(fill = group, shape = run), size = 6, stroke = 0.8, colour = "grey20") +
    scale_fill_manual(values = GROUP_COLS, drop = FALSE) +
    scale_shape_manual(values = c(C3SVT8 = 21, VJ385X = 24), drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 5)),
           shape = guide_legend(override.aes = list(fill = "grey70", size = 5))) +
    labs(title = title,
         x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])) +
    base_theme
  list(plot = g, gap = gap)
}
cat("\n")
A <- ordination(lcpm,     "Before run adjustment")
B <- ordination(lcpm_adj, "After run adjustment")
cat(sprintf("\npaired libraries move %.1f-fold closer after adjustment\n",
            median(A$gap) / median(B$gap)))

# ---- 3. panel C: each animal against itself ----------------------------------
pairs_df <- do.call(rbind, lapply(seq_along(RUN_A), function(i)
  data.frame(animal = sprintf("hTau vehicle %d", i),
             runA = lcpm[, RUN_A[i]], runB = lcpm[, RUN_B[i]])))
lab <- pairs_df %>% group_by(animal) %>%
  summarise(r = cor(runA, runB), rho = cor(runA, runB, method = "spearman"),
            median_difference = median(runB - runA),
            .groups = "drop") %>%
  mutate(text = sprintf("r = %.3f", r))
cat("\nwithin-animal agreement across the two runs:\n")
print(as.data.frame(lab %>% mutate(across(where(is.numeric), ~ round(., 4)))),
      row.names = FALSE)

pC <- ggplot(pairs_df, aes(runA, runB)) +
  geom_abline(slope = 1, intercept = 0, linetype = "42", colour = MUTED, linewidth = 0.4) +
  geom_point(colour = "#FFE0C0", size = 0.5, alpha = 0.45) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = text), inherit.aes = FALSE,
            hjust = -0.15, vjust = 1.6, family = FAMILY, size = 5.3, colour = INK) +
  facet_wrap(~ animal, nrow = 1) +
  labs(title = "Same animal, two sequencing runs",
       x = expression("log"[2]*" CPM, run C3SVT8"),
       y = expression("log"[2]*" CPM, run VJ385X")) +
  base_theme + theme(strip.background = element_blank(),
                     strip.text = element_text(size = FONT - 4, family = FAMILY))

# ---- 4. assemble -------------------------------------------------------------
save_all <- function(p, w, h, stem) {
  tryCatch(ggsave(paste0(stem, ".pdf"), p, width = w, height = h, device = cairo_pdf),
    error = function(e) { message("cairo_pdf failed: ", conditionMessage(e))
      ggsave(paste0(stem, ".pdf"), p, width = w, height = h) })
  ggsave(paste0(stem, ".png"), p, width = w, height = h, dpi = 400)
}

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  fig <- (A$plot + B$plot) / pC +
    plot_layout(heights = c(1.15, 0.85), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(family = FAMILY, face = "bold", size = FONT + 4),
          legend.position = "right")
  save_all(fig, FIG_W, FIG_H, "FigS1_QC_batch_technical")
  cat("\nwritten: FigS1_QC_batch_technical.pdf / .png\n")
} else {
  cat("\npatchwork is not installed; panels were written separately.\n")
  save_all(A$plot, 7.5, 6.5, "FigS1_QC_pca_before")
  save_all(B$plot, 7.5, 6.5, "FigS1_QC_pca_after")
  save_all(pC, FIG_W, 5, "FigS1_QC_paired_runs")
}

# ---- 5. source data ----------------------------------------------------------
# base indexing rather than dplyr::select, so that the script does not depend on
# which packages happen to be attached (AnnotationDbi masks select)
write.csv(as.data.frame(lab)[, c("animal", "r", "rho", "median_difference")],
          "FigS1_QC_correlations.csv", row.names = FALSE)

cat("\nwritten to", normalizePath(OUT), ":\n FigS1_QC_correlations.csv\n")
cat("Displayed values are run-adjusted where stated; all reported statistics\n",
    "come from the model, where the run effect is a term in the design.\n", sep = "")
