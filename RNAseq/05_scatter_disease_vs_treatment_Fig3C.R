#!/usr/bin/env Rscript
# =============================================================================
# Figure 3c - Transcriptome-wide relationship between disease and treatment
#
# Each point is one gene. The x axis is the effect of the tauopathy
# (hTau vehicle vs nTg vehicle); the y axis is the effect of RA013915 within
# hTau mice (hTau RA013915 vs hTau vehicle). Genes that are suppressed by the
# pathology and raised again by the drug fall in the upper-left quadrant, so
# rescue is a position in the plot rather than a derived score.
#
# A volcano plot cannot show this: it displays a single contrast, and the
# rescue claim is a statement about two contrasts at once.
#
# IMPORTANT CAVEAT, also stated in the figure legend
# The two axes share the hTau vehicle group. Part of any negative correlation
# between them is therefore induced by that shared denominator rather than by
# biology. The companion panel (07_scatter_treatment_specificity_Fig3E.R) does not
# have this property, because its two axes come from disjoint sets of animals.
#
# Input : RNAseq_20lib_results.xlsx, sheet "Full_Merged_Results"
#         produced by the limma-voom + duplicateCorrelation model
# Output: Fig3c_scatter_rescue.pdf / .png, Fig3c_source_data.csv
# =============================================================================

rm(list = ls())
library(ggplot2); library(dplyr); library(ggrepel); library(readxl)

# ---- PATHS -------------------------------------------------------------------
IN  <- "RNAseq_20lib_results.xlsx"
OUT <- "."
stopifnot(file.exists(IN), dir.exists(OUT)); setwd(OUT)

# ---- GENE SETS ---------------------------------------------------------------
# Rescued genes, strict definition: significant in BOTH contrasts (FDR < 0.05)
# with opposite signs, i.e. down with pathology and up with treatment.
RESCUED <- c("Apold1", "Btg2", "Dusp1", "Fos", "Ier2", "Klf2")

# Drug-responsive but not rescued: up with treatment, not altered by pathology.
TREATED_ONLY <- "Zfp36"

# Same direction as the rescued genes in both contrasts but above the FDR
# threshold under this model. Showing them makes clear that the boundary of the
# rescued set is a statistical threshold and not a discontinuity in the biology.
# Set SHOW_SUBTHRESHOLD to FALSE for a cleaner panel.
SUBTHRESHOLD <- c("Ccn1", "Cldn5", "Klf4")
SHOW_SUBTHRESHOLD <- TRUE

# ---- STYLE -------------------------------------------------------------------
# Three fills for three different kinds of evidence. A single colour for
# everything that is not rescued would place Zfp36, which is significant with
# treatment, in the same category as genes that are not, and a reader checking
# the source data would find two different things under one symbol.
ORANGE_DARK <- "#FF6A00"   # rescued: FDR < 0.05 in both contrasts
ORANGE_MID  <- "#FFA559"   # treatment-responsive, not altered by the pathology
ORANGE_PALE <- "#FFFFFF"   # same direction, above the FDR threshold (open symbol)
EDGE <- "#8a3a00"; INK <- "#0b0b0b"; MUTED <- "#52514e"; CLOUD <- "#d9d9d6"
FAMILY <- "Arial"
FONT <- 20      # axis text and axis titles, in points
GENE <- 6.3     # gene labels, in ggplot size units (points / 2.845)
                #   4.2 ~ 12 pt | 5.3 ~ 15 pt | 6.3 ~ 18 pt | 7.0 ~ 20 pt
PT <- 4.2       # diameter of highlighted points
AX <- 2         # axis limits; 99.4% of genes fall within 1.5 on both axes
FIG_W <- 7.4; FIG_H <- 7.0

if (requireNamespace("extrafont", quietly = TRUE))
  suppressWarnings(try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE))

# ---- 1. data -----------------------------------------------------------------
d <- as.data.frame(read_excel(IN, sheet = "Full_Merged_Results"))
names(d)[1] <- "gene"
NEEDED <- c("log2FoldChange_disease", "log2FoldChange_drug_hTau")
if (!all(NEEDED %in% names(d)))
  stop("missing columns: ", paste(setdiff(NEEDED, names(d)), collapse = ", "))

d <- d %>% filter(!is.na(log2FoldChange_disease), !is.na(log2FoldChange_drug_hTau))
cat(sprintf("genes with both contrasts: %d\n", nrow(d)))

highlight <- c(RESCUED, TREATED_ONLY, if (SHOW_SUBTHRESHOLD) SUBTHRESHOLD)
missing <- setdiff(highlight, d$gene)
if (length(missing)) warning("not found in results: ", paste(missing, collapse = ", "), call. = FALSE)

d <- d %>% mutate(class = factor(
  ifelse(gene %in% RESCUED, "rescued",
  ifelse(gene %in% TREATED_ONLY, "treated",
  ifelse(SHOW_SUBTHRESHOLD & gene %in% SUBTHRESHOLD, "subthreshold", "background"))),
  levels = c("background", "subthreshold", "treated", "rescued")))
print(table(d$class))

# ---- 2. quadrants and correlation, for the text ------------------------------
q <- d %>% summarise(
  rescue_quadrant = sum(log2FoldChange_disease < 0 & log2FoldChange_drug_hTau > 0),
  opposite        = sum(log2FoldChange_disease > 0 & log2FoldChange_drug_hTau < 0))
cat(sprintf("rescue quadrant: %d genes (%.1f%%) | opposite quadrant: %d (%.1f%%)\n",
            q$rescue_quadrant, 100 * q$rescue_quadrant / nrow(d),
            q$opposite, 100 * q$opposite / nrow(d)))
rho <- cor(d$log2FoldChange_disease, d$log2FoldChange_drug_hTau, method = "spearman")
cat(sprintf("Spearman (disease vs treatment): %+0.3f\n", rho))
cat("  NB: the two contrasts share the hTau vehicle group, so part of this\n",
    "  correlation is induced by the design rather than by biology.\n", sep = "")

# ---- 3. plot -----------------------------------------------------------------
# Genes outside the axis limits are drawn as triangles on the border rather
# than dropped, so nothing is hidden by the choice of scale.
d <- d %>% mutate(
  x = pmax(pmin(log2FoldChange_disease,   AX), -AX),
  y = pmax(pmin(log2FoldChange_drug_hTau, AX), -AX),
  clipped = abs(log2FoldChange_disease) > AX | abs(log2FoldChange_drug_hTau) > AX)
cat(sprintf("genes drawn on the border (|log2FC| > %g): %d\n", AX, sum(d$clipped)))

p <- ggplot(d, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.3) +
  geom_point(data = . %>% filter(class == "background", !clipped),
             colour = CLOUD, size = 0.5, alpha = 0.5) +
  geom_point(data = . %>% filter(class == "background", clipped),
             colour = CLOUD, size = 0.9, shape = 17, alpha = 0.7) +
  geom_point(data = . %>% filter(class == "subthreshold"), shape = 21,
             fill = ORANGE_PALE, colour = ORANGE_DARK, stroke = 1.0, size = PT - 0.6) +
  geom_point(data = . %>% filter(class == "treated"), shape = 21,
             fill = ORANGE_MID, colour = EDGE, stroke = 0.5, size = PT - 0.3) +
  geom_point(data = . %>% filter(class == "rescued"), shape = 21,
             fill = ORANGE_DARK, colour = EDGE, stroke = 0.5, size = PT) +
  ggrepel::geom_text_repel(
    data = . %>% filter(class != "background"), aes(label = gene),
    family = FAMILY, size = GENE, colour = INK, fontface = "italic",
    min.segment.length = 0, segment.colour = MUTED, segment.size = 0.3,
    box.padding = 0.5, point.padding = 0.25, max.overlaps = Inf, seed = 1) +
  # placed bottom-right: the upper-left corner is occupied by the rescued genes
  annotate("text", x = AX * 0.96, y = -AX * 0.94, hjust = 1, vjust = 0,
           label = sprintf("rho == %+0.3f", rho), parse = TRUE,
           family = FAMILY, size = GENE * 0.8, colour = MUTED) +
  scale_x_continuous(breaks = seq(-AX, AX, 1)) +
  scale_y_continuous(breaks = seq(-AX, AX, 1)) +
  coord_equal(xlim = c(-AX, AX), ylim = c(-AX, AX), expand = FALSE, clip = "off") +
  labs(x = expression("Disease  "*log[2]*" fold change (hTau vs nTg)"),
       y = expression("RA013915 in hTau  "*log[2]*" fold change")) +
  theme_classic(base_size = FONT, base_family = FAMILY) +
  theme(text = element_text(colour = INK, family = FAMILY),
        axis.text  = element_text(colour = INK, size = FONT),
        axis.title = element_text(size = FONT),
        legend.position = "none")

# cairo_pdf rather than pdf(): on Windows the pdf device does not know system
# fonts unless they have been registered with extrafont, and fails on "Arial".
ok <- tryCatch({ ggsave("Fig3c_scatter_rescue.pdf", p, width = FIG_W, height = FIG_H,
                        device = cairo_pdf); TRUE },
  error = function(e) { message("cairo_pdf failed: ", conditionMessage(e),
                                " - retrying with the default device")
    ggsave("Fig3c_scatter_rescue.pdf", p, width = FIG_W, height = FIG_H); TRUE })
ggsave("Fig3c_scatter_rescue.png", p, width = FIG_W, height = FIG_H, dpi = 400)

# ---- 4. source data ----------------------------------------------------------
src <- d %>% filter(class != "background") %>%
  mutate(category = ifelse(gene %in% RESCUED, "rescued",
                    ifelse(gene %in% TREATED_ONLY, "treatment-responsive",
                           "same direction, above FDR threshold"))) %>%
  dplyr::select(gene, category,
                log2FoldChange_disease,   padj_disease,
                log2FoldChange_drug_hTau, padj_drug_hTau) %>%
  arrange(category, desc(log2FoldChange_drug_hTau))
print(src %>% mutate(across(where(is.numeric), ~ signif(., 3))), row.names = FALSE)
write.csv(src, "Fig3c_source_data.csv", row.names = FALSE)

cat("\nwritten to", normalizePath(OUT), ":\n",
    " Fig3c_scatter_rescue.pdf\n Fig3c_scatter_rescue.png\n Fig3c_source_data.csv\n")
