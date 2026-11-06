#!/usr/bin/env Rscript
# =============================================================================
# Figure 3a - Genes altered by the tauopathy
#
# Companion to the GO panel: this one names the genes, that one names the
# processes they belong to. Together they say what is changed and what it does.
#
# A volcano is the right display here and the wrong one for the treatment
# contrast. The disease contrast has 292 differentially expressed genes and a
# wide range of effect sizes, so the plot is populated and the labelled genes
# stand out from a real cloud. The treatment contrast has seven, and a volcano
# of seven points carries less than the scatter in 03_, which also shows how
# each gene behaved in the disease.
#
# Mapt is the top gene. It is the transgene, so its position is a check that
# the pipeline and the group assignment are correct rather than a finding, and
# it is labelled as such in the figure legend.
#
# Input : RNAseq_20lib_results.xlsx (from 01_)
# Output: Fig3a_volcano_disease.pdf / .png, Fig3a_volcano_source.csv
# =============================================================================

rm(list = ls())
library(readxl); library(dplyr); library(ggplot2); library(ggrepel)

# ---- PATHS -------------------------------------------------------------------
OUT <- "."
stopifnot(dir.exists(OUT)); setwd(OUT)
IN <- "RNAseq_20lib_results.xlsx"
stopifnot(file.exists(IN))

FDR <- 0.05

# ---- GENES TO NAME -----------------------------------------------------------
# Chosen for what they represent, not only for their rank: the transgene, the
# activity-dependent programme that the treatment later acts on, and the
# microglial and complement response that is the expected signature of a
# tauopathy. Every significant gene is in the source data.
LABEL_DOWN <- c("Mapt", "Fos", "Dusp1", "Btg2", "Arc", "Ier2", "Junb",
                "Nr4a1", "Egr4", "Klf2", "Apold1", "Dcx")
LABEL_UP   <- c("Trem2", "C1qa", "C1qb", "Ctss", "Fcer1g", "Il33")
LABELS <- c(LABEL_DOWN, LABEL_UP)

# ---- STYLE -------------------------------------------------------------------
# Same two colours as the GO panel and the heatmap: blue is lower in hTau, red
# is higher. One colour system across the whole figure.
COL_DOWN <- "#2166AC"; COL_UP <- "#B2182B"
CLOUD <- "#d9d9d6"; INK <- "#0b0b0b"; MUTED <- "#52514e"
FAMILY <- "Arial"
FONT <- 20
GENE <- 6.0     # gene labels, ggplot size units (points / 2.845)
PT   <- 2.6
AX   <- 3       # |log2FC| axis limit; 99.8% of genes fall within 2.5
YMAX <- 12      # -log10(FDR) limit; the maximum in these data is 11.9
FIG_W <- 7.6; FIG_H <- 7.5

if (requireNamespace("extrafont", quietly = TRUE))
  suppressWarnings(try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE))

# ---- 1. data -----------------------------------------------------------------
d <- as.data.frame(read_excel(IN, sheet = "Full_Merged_Results"))
names(d)[1] <- "gene"
d <- d %>% filter(!is.na(padj_disease), !is.na(log2FoldChange_disease)) %>%
  mutate(y_raw = -log10(padj_disease),
         class = factor(ifelse(padj_disease >= FDR, "ns",
                        ifelse(log2FoldChange_disease < 0, "down", "up")),
                        levels = c("ns", "down", "up")))
cat(sprintf("genes tested: %d\n", nrow(d)))
print(table(d$class))

missing <- setdiff(LABELS, d$gene)
if (length(missing)) warning("not found: ", paste(missing, collapse = ", "), call. = FALSE)
ns_labelled <- intersect(LABELS, d$gene[d$class == "ns"])
if (length(ns_labelled))
  cat("labelled but not significant, so shown in grey:",
      paste(ns_labelled, collapse = ", "), "\n")

# Points beyond the axis limits are drawn as triangles on the border rather
# than dropped, so the choice of scale hides nothing.
d <- d %>% mutate(
  x = pmax(pmin(log2FoldChange_disease, AX), -AX),
  y = pmin(y_raw, YMAX),
  clipped = abs(log2FoldChange_disease) > AX | y_raw > YMAX,
  lab = ifelse(gene %in% LABELS, gene, NA_character_))
cat(sprintf("genes drawn on the border: %d\n", sum(d$clipped)))

# ---- 2. panel ----------------------------------------------------------------
p <- ggplot(d, aes(x, y)) +
  geom_hline(yintercept = -log10(FDR), linetype = "22", colour = MUTED, linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.3) +
  geom_point(data = . %>% filter(class == "ns", !clipped),
             colour = CLOUD, size = 0.5, alpha = 0.5) +
  geom_point(data = . %>% filter(class == "ns", clipped),
             colour = CLOUD, size = 0.9, shape = 17, alpha = 0.7) +
  geom_point(data = . %>% filter(class != "ns", !clipped),
             aes(fill = class), shape = 21, colour = "grey25", stroke = 0.25, size = PT) +
  geom_point(data = . %>% filter(class != "ns", clipped),
             aes(fill = class), shape = 24, colour = "grey25", stroke = 0.25, size = PT) +
  ggrepel::geom_text_repel(
    data = . %>% filter(!is.na(lab)), aes(label = lab),
    family = FAMILY, size = GENE, colour = INK, fontface = "italic",
    min.segment.length = 0, segment.colour = MUTED, segment.size = 0.3,
    box.padding = 0.45, point.padding = 0.2, max.overlaps = Inf, seed = 1) +
  scale_fill_manual(values = c(down = COL_DOWN, up = COL_UP), guide = "none") +
  scale_x_continuous(breaks = seq(-AX, AX, 1)) +
  coord_cartesian(xlim = c(-AX, AX), ylim = c(0, YMAX), expand = FALSE, clip = "off") +
  labs(x = expression("Disease  "*log[2]*" fold change (hTau vs nTg)"),
       y = expression(-log[10]~"adjusted "*italic(P))) +
  theme_classic(base_size = FONT, base_family = FAMILY) +
  theme(text = element_text(colour = INK, family = FAMILY),
        axis.text = element_text(colour = INK, size = FONT),
        axis.title = element_text(size = FONT),
        legend.position = "none")

tryCatch(ggsave("Fig3a_volcano_disease.pdf", p, width = FIG_W, height = FIG_H, device = cairo_pdf),
  error = function(e) { message("cairo_pdf failed: ", conditionMessage(e))
    ggsave("Fig3a_volcano_disease.pdf", p, width = FIG_W, height = FIG_H) })
ggsave("Fig3a_volcano_disease.png", p, width = FIG_W, height = FIG_H, dpi = 400)

# ---- 3. source data ----------------------------------------------------------
src <- d %>% filter(class != "ns") %>%
  dplyr::select(gene, log2FoldChange_disease, pvalue_disease, padj_disease,
                baseMean_disease, class) %>%
  arrange(padj_disease)
write.csv(src, "Fig3a_volcano_source.csv", row.names = FALSE)
cat(sprintf("\nsignificant genes written to source data: %d\n", nrow(src)))
cat("\nlabelled genes:\n")
print(d %>% filter(!is.na(lab)) %>%
        transmute(gene, log2FC = round(log2FoldChange_disease, 3),
                  padj = signif(padj_disease, 3), class) %>%
        arrange(padj), row.names = FALSE)

cat("\nwritten to", normalizePath(OUT), ":\n",
    " Fig3a_volcano_disease.pdf\n Fig3a_volcano_disease.png\n Fig3a_volcano_source.csv\n")
