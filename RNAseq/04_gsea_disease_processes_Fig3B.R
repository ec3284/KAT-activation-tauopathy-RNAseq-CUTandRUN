#!/usr/bin/env Rscript
# =============================================================================
# Figure 3b - Biological processes altered by the tauopathy
#
# Gene set enrichment on the disease contrast, using every gene in the ranked
# transcriptome and no significance threshold.
#
# WHY NOT OVER-REPRESENTATION ON THE DIFFERENTIALLY EXPRESSED GENES
# That approach takes the 159 genes below the FDR threshold and asks which GO
# terms are over-represented among them. Those 159 are dominated by a small
# core of immediate-early genes, which are among the most studied genes in
# biology and are therefore annotated to hundreds of GO terms covering every
# context in which anyone has measured them. Three of them falling inside a
# term is enough to make it significant, and the resulting list is full of
# terms whose names describe the experiment they were annotated in rather than
# anything happening in this tissue. Ranking the whole transcriptome avoids
# both the threshold and that bias.
#
# Enrichment is fgsea throughout this paper, the same method in every panel.
#
# Input : fgsea_results.xlsx (from 02_)
# Output: Fig3b_GSEA_disease.pdf / .png, Fig3b_GSEA_disease_source.csv
# =============================================================================

rm(list = ls())
library(readxl); library(dplyr); library(ggplot2)

# ---- PATHS -------------------------------------------------------------------
OUT <- "."
stopifnot(dir.exists(OUT)); setwd(OUT)
IN <- "fgsea_results.xlsx"
stopifnot(file.exists(IN))

# ---- SETTINGS ----------------------------------------------------------------
FDR         <- 0.05
DIRECTION   <- "down"   # "down" (NES < 0), "up" (NES > 0), or "both"
N_CANDIDATE <- 60       # terms kept after removing redundant ones
N_SHOW      <- 10       # terms drawn
MAX_OVERLAP <- 0.5      # terms sharing more of their leading edge than this are dropped

# GO names are written for machines. Short display names keep the panel legible
# at this type size; the full GO identifier stays in the source data.
SHORT <- c(
  GOBP_SYNAPSE_ASSEMBLY                          = "Synapse assembly",
  GOBP_SYNAPSE_ORGANIZATION                      = "Synapse organization",
  GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING    = "Trans-synaptic signaling",
  GOBP_MODULATION_OF_CHEMICAL_SYNAPTIC_TRANSMISSION = "Synaptic transmission",
  GOBP_DENDRITE_DEVELOPMENT                      = "Dendrite development",
  GOBP_DENDRITE_MORPHOGENESIS                    = "Dendrite morphogenesis",
  GOBP_AXON_DEVELOPMENT                          = "Axon development",
  GOBP_AXONOGENESIS                              = "Axonogenesis",
  GOBP_NEURON_PROJECTION_GUIDANCE                = "Neuron projection guidance",
  GOBP_REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT = "Neuron projection development",
  GOBP_LEARNING_OR_MEMORY                        = "Learning and memory",
  GOBP_LEARNING                                  = "Learning",
  GOBP_REGULATION_OF_MEMBRANE_POTENTIAL          = "Membrane potential",
  GOBP_ACTION_POTENTIAL                          = "Action potential",
  GOBP_SYNAPTIC_VESICLE_EXOCYTOSIS               = "Synaptic vesicle exocytosis",
  GOBP_NEUROTRANSMITTER_SECRETION                = "Neurotransmitter secretion",
  GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY = "Synapse structure and activity",
  GOBP_POSTSYNAPSE_ORGANIZATION                  = "Postsynapse organization",
  GOBP_FOREBRAIN_DEVELOPMENT                     = "Forebrain development",
  GOBP_RESPONSE_TO_CAMP                          = "Response to cAMP",
  GOBP_OXIDATIVE_PHOSPHORYLATION                 = "Oxidative phosphorylation",
  GOBP_MITOCHONDRIAL_TRANSLATION                 = "Mitochondrial translation",
  GOBP_RIBOSOME_BIOGENESIS                       = "Ribosome biogenesis")

# ---- STYLE -------------------------------------------------------------------
COL_DOWN <- "#2166AC"; COL_UP <- "#B2182B"
INK <- "#0b0b0b"; MUTED <- "#52514e"
FAMILY <- "Arial"
FONT <- 20
LABEL_WRAP <- 30
# Panel size in inches. Raise FIG_H to make it taller, FIG_W to make it wider;
# the rows space themselves out over whatever height is given, so a taller
# figure also separates the ten terms more.
FIG_W <- 13.5; FIG_H <- 8.5

if (requireNamespace("extrafont", quietly = TRUE))
  suppressWarnings(try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE))
wrap_text <- function(x, n)
  vapply(x, function(s) paste(strwrap(s, width = n), collapse = "\n"), "", USE.NAMES = FALSE)

# ---- 1. candidate terms from fgsea -------------------------------------------
g <- as.data.frame(read_excel(IN, sheet = "GOBP_disease"))
cat(sprintf("GO BP terms tested: %d | significant in the disease contrast: %d\n",
            nrow(g), sum(g$padj_disease < FDR)))

cand <- g %>% filter(padj_disease < FDR)
if (DIRECTION == "down") cand <- cand %>% filter(NES_disease < 0)
if (DIRECTION == "up")   cand <- cand %>% filter(NES_disease > 0)
cand <- cand %>% arrange(padj_disease)
cat(sprintf("after keeping direction '%s': %d\n", DIRECTION, nrow(cand)))

# GO contains many near-duplicate terms. Of two terms whose leading edges, the
# genes that actually drive the enrichment, overlap by more than MAX_OVERLAP,
# only the more significant is kept.
le <- function(x) strsplit(ifelse(is.na(x), "", x), "/")[[1]]
keep <- integer(0); kept <- list()
for (i in seq_len(nrow(cand))) {
  s <- le(cand$leadingEdge_disease[i])
  if (length(s) < 3) next
  if (all(vapply(kept, function(k)
        length(intersect(s, k)) / min(length(s), length(k)) <= MAX_OVERLAP, logical(1)))) {
    kept <- c(kept, list(s)); keep <- c(keep, i)
  }
  if (length(keep) == N_CANDIDATE) break
}
cand <- cand[keep, ]
cat(sprintf("after removing overlapping leading edges: %d candidates\n", nrow(cand)))
stopifnot(nrow(cand) > 0)

d <- cand
stopifnot(nrow(d) > 0)

# ---- 2. panel ----------------------------------------------------------------
pl <- d %>% arrange(padj_disease) %>% slice_head(n = N_SHOW) %>%
  mutate(
    short = ifelse(pathway %in% names(SHORT), SHORT[pathway],
                   paste0(toupper(substr(gsub("_", " ", tolower(sub("^GOBP_", "", pathway))), 1, 1)),
                          substr(gsub("_", " ", tolower(sub("^GOBP_", "", pathway))), 2, 100))),
    lab = wrap_text(short, LABEL_WRAP)) %>%
  arrange(NES_disease) %>%
  mutate(lab = factor(make.unique(lab), levels = make.unique(lab)))

no_short <- setdiff(pl$pathway, names(SHORT))
if (length(no_short))
  cat("\nno short label for these, the GO name is used:\n  ",
      paste(no_short, collapse = "\n  "), "\n", sep = "")

cat("\nterms drawn:\n")
print(pl %>% transmute(short, NES = round(NES_disease, 2),
                       fgsea_FDR = signif(padj_disease, 3),
                       genes = size_disease),
      row.names = FALSE)
p <- ggplot(pl, aes(NES_disease, lab)) +
  geom_vline(xintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = NES_disease, yend = lab),
               colour = "#c9c9c4", linewidth = 0.5) +
  geom_point(aes(size = size_disease, fill = NES_disease > 0),
             shape = 21, colour = "grey25", stroke = 0.5) +
  scale_fill_manual(values = c(`FALSE` = COL_DOWN, `TRUE` = COL_UP), guide = "none") +
  scale_size_continuous(range = c(4, 13), name = "Genes",
                        breaks = scales::breaks_pretty(4)) +
  # extra room on the side the points sit on, so the largest dot does not touch
  # the panel edge and the axis reaches past it
  scale_x_continuous(expand = expansion(mult = if (DIRECTION == "up")
                                          c(0.04, 0.14) else c(0.14, 0.04))) +
  labs(x = "Normalized enrichment score (NES)", y = NULL) +
  theme_classic(base_size = FONT, base_family = FAMILY) +
  theme(text = element_text(colour = INK, family = FAMILY),
        axis.text = element_text(colour = INK, size = FONT),
        axis.title.x = element_text(size = FONT),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = FONT - 2),
        legend.text = element_text(size = FONT - 2))

tryCatch(ggsave("Fig3b_GSEA_disease.pdf", p, width = FIG_W, height = FIG_H, device = cairo_pdf),
  error = function(e) { message("cairo_pdf failed: ", conditionMessage(e))
    ggsave("Fig3b_GSEA_disease.pdf", p, width = FIG_W, height = FIG_H) })
ggsave("Fig3b_GSEA_disease.png", p, width = FIG_W, height = FIG_H, dpi = 400)

# ---- 3. source data ----------------------------------------------------------
write.csv(d %>% dplyr::select(pathway, size_disease, NES_disease, pval_disease,
                              padj_disease, leadingEdge_disease) %>%
            arrange(padj_disease),
          "Fig3b_GSEA_disease_source.csv", row.names = FALSE)

cat("\nwritten to", normalizePath(OUT), ":\n",
    " Fig3b_GSEA_disease.pdf\n Fig3b_GSEA_disease.png\n Fig3b_GSEA_disease_source.csv\n")
