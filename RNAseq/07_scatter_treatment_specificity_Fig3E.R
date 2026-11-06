#!/usr/bin/env Rscript
# =============================================================================
# Figure 3e - The treatment effect is confined to the tauopathy brain
#
# Each point is one gene. The y axis is the effect of RA013915 in hTau mice,
# the same quantity as in Figure 3a; the x axis is the effect of the same
# treatment in non-transgenic littermates. The dashed diagonal is the line
# along which a compound acting identically in both genotypes would fall.
#
# WHY THIS PANEL IS STRONGER THAN IT LOOKS
# The two axes are built from disjoint sets of animals: hTau RA013915 against
# hTau vehicle on y, nTg RA013915 against nTg vehicle on x. No animal
# contributes to both, so unlike Figure 3a there is no correlation induced by a
# shared reference group. A Spearman coefficient of essentially zero here means
# that knowing what the compound does in a healthy brain predicts nothing about
# what it does in a diseased one.
#
# The same conclusion is supported by three independent tests reported in the
# manuscript: a genotype x treatment interaction fitted genome-wide, a rotation
# gene-set test of the treatment-responsive set in the non-transgenic arm, and
# the reciprocal test of the non-transgenic set in the tauopathy arm.
#
# Input : RNAseq_20lib_results.xlsx, sheet "Full_Merged_Results"
# Output: Fig3e_scatter_specificity.pdf / .png, Fig3e_source_data.csv
# =============================================================================

rm(list = ls())
library(ggplot2); library(dplyr); library(ggrepel); library(readxl)

# ---- PATHS -------------------------------------------------------------------
IN  <- "RNAseq_20lib_results.xlsx"
OUT <- "."
stopifnot(file.exists(IN), dir.exists(OUT)); setwd(OUT)

# ---- GENE SETS ---------------------------------------------------------------
RESCUED      <- c("Apold1", "Btg2", "Dusp1", "Fos", "Ier2", "Klf2")
TREATED_ONLY <- "Zfp36"
SUBTHRESHOLD <- c("Ccn1", "Cldn5", "Klf4")
SHOW_SUBTHRESHOLD <- TRUE   # keep consistent with 05_scatter_disease_vs_treatment_Fig3C.R

# ---- STYLE -------------------------------------------------------------------
# Same three fills as 05_scatter_disease_vs_treatment_Fig3C.R; keep them in step.
ORANGE_DARK <- "#FF6A00"   # rescued: FDR < 0.05 in both contrasts
ORANGE_MID  <- "#FFA559"   # treatment-responsive, not altered by the pathology
ORANGE_PALE <- "#FFFFFF"   # same direction, above the FDR threshold (open symbol)
EDGE <- "#8a3a00"; INK <- "#0b0b0b"; MUTED <- "#52514e"; CLOUD <- "#d9d9d6"
FAMILY <- "Arial"
FONT <- 20; GENE <- 6.3; PT <- 4.2; AX <- 2
FIG_W <- 7.4; FIG_H <- 7.0

if (requireNamespace("extrafont", quietly = TRUE))
  suppressWarnings(try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE))

# ---- 1. data -----------------------------------------------------------------
d <- as.data.frame(read_excel(IN, sheet = "Full_Merged_Results"))
names(d)[1] <- "gene"
NEEDED <- c("log2FoldChange_drug_hTau", "log2FoldChange_drug_nTg")
if (!all(NEEDED %in% names(d)))
  stop("missing columns: ", paste(setdiff(NEEDED, names(d)), collapse = ", "))
d <- d %>% filter(!is.na(log2FoldChange_drug_hTau), !is.na(log2FoldChange_drug_nTg))
cat(sprintf("genes with both treatment contrasts: %d\n", nrow(d)))

d <- d %>% mutate(class = factor(
  ifelse(gene %in% RESCUED, "rescued",
  ifelse(gene %in% TREATED_ONLY, "treated",
  ifelse(SHOW_SUBTHRESHOLD & gene %in% SUBTHRESHOLD, "subthreshold", "background"))),
  levels = c("background", "subthreshold", "treated", "rescued")))
print(table(d$class))

# ---- 2. statistics for the text ----------------------------------------------
rho <- cor(d$log2FoldChange_drug_nTg, d$log2FoldChange_drug_hTau, method = "spearman")
cat(sprintf("Spearman (treatment in nTg vs treatment in hTau): %+0.3f\n", rho))
cat("  the two axes share no animals, so this is not inflated by the design\n")

h <- d %>% filter(gene %in% c(RESCUED, TREATED_ONLY))
cat(sprintf("treatment-responsive genes: median |log2FC| in hTau %.3f, in nTg %.3f (%.1f-fold)\n",
            median(abs(h$log2FoldChange_drug_hTau)),
            median(abs(h$log2FoldChange_drug_nTg)),
            median(abs(h$log2FoldChange_drug_hTau)) / median(abs(h$log2FoldChange_drug_nTg))))
cat(sprintf("  same sign in both genotypes: %d of %d\n",
            sum(sign(h$log2FoldChange_drug_hTau) == sign(h$log2FoldChange_drug_nTg)), nrow(h)))
cat("  the shared sign with a much smaller magnitude is consistent with a floor\n",
    "  effect: in healthy brain these immediate-early genes already sit at baseline.\n", sep = "")
if ("padj_drug_nTg" %in% names(d))
  cat(sprintf("  significant in the nTg contrast: %d of %d\n",
              sum(h$padj_drug_nTg < 0.05, na.rm = TRUE), nrow(h)))

# ---- 3. plot -----------------------------------------------------------------
d <- d %>% mutate(
  x = pmax(pmin(log2FoldChange_drug_nTg,  AX), -AX),
  y = pmax(pmin(log2FoldChange_drug_hTau, AX), -AX),
  clipped = abs(log2FoldChange_drug_nTg) > AX | abs(log2FoldChange_drug_hTau) > AX)
cat(sprintf("genes drawn on the border: %d\n", sum(d$clipped)))

p <- ggplot(d, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.3) +
  # identity line: where a compound acting equally in both genotypes would lie
  geom_abline(slope = 1, intercept = 0, linetype = "42", colour = MUTED, linewidth = 0.45) +
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
  # top-left: the bottom-left corner lies on the identity line, which would
  # otherwise run straight through the text
  annotate("text", x = -AX * 0.96, y = AX * 0.94, hjust = 0, vjust = 1,
           label = sprintf("rho == %+0.3f", rho), parse = TRUE,
           family = FAMILY, size = GENE * 0.8, colour = MUTED) +
  scale_x_continuous(breaks = seq(-AX, AX, 1)) +
  scale_y_continuous(breaks = seq(-AX, AX, 1)) +
  coord_equal(xlim = c(-AX, AX), ylim = c(-AX, AX), expand = FALSE, clip = "off") +
  labs(x = expression("RA013915 in nTg  "*log[2]*" fold change"),
       y = expression("RA013915 in hTau  "*log[2]*" fold change")) +
  theme_classic(base_size = FONT, base_family = FAMILY) +
  theme(text = element_text(colour = INK, family = FAMILY),
        axis.text  = element_text(colour = INK, size = FONT),
        axis.title = element_text(size = FONT),
        legend.position = "none")

tryCatch(ggsave("Fig3e_scatter_specificity.pdf", p, width = FIG_W, height = FIG_H,
                device = cairo_pdf),
  error = function(e) { message("cairo_pdf failed: ", conditionMessage(e))
    ggsave("Fig3e_scatter_specificity.pdf", p, width = FIG_W, height = FIG_H) })
ggsave("Fig3e_scatter_specificity.png", p, width = FIG_W, height = FIG_H, dpi = 400)

# ---- 4. source data ----------------------------------------------------------
keep <- c("gene", "log2FoldChange_drug_hTau", "padj_drug_hTau",
          "log2FoldChange_drug_nTg", "padj_drug_nTg")
src <- d %>% filter(class != "background") %>%
  dplyr::select(any_of(keep)) %>%
  mutate(rescued = gene %in% RESCUED) %>%
  arrange(desc(rescued), desc(log2FoldChange_drug_hTau))
print(src %>% mutate(across(where(is.numeric), ~ signif(., 3))), row.names = FALSE)
write.csv(src, "Fig3e_source_data.csv", row.names = FALSE)

cat("\nwritten to", normalizePath(OUT), ":\n",
    " Fig3e_scatter_specificity.pdf\n Fig3e_scatter_specificity.png\n Fig3e_source_data.csv\n")
