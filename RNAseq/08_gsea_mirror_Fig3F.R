#!/usr/bin/env Rscript
# =============================================================================
# Figure 3f - Gene sets that the pathology and the treatment move in opposite
# directions
#
# This panel answers the objection that the treatment changes only a handful of
# genes. Enrichment uses the whole ranked transcriptome and no threshold, so a
# coordinated shift inside a pathway is visible even where no single gene in it
# reaches significance.
#
# Terms are kept only if they are significant in BOTH contrasts and move in
# OPPOSITE directions. Redundant terms are then removed by leading-edge
# overlap, because GO contains many near-duplicate terms and showing them all
# would inflate the panel without adding evidence.
#
# Input : fgsea_results.xlsx (from 02_fgsea_both_contrasts.R)
# Output: Fig3f_GSEA_mirror.pdf / .png, Fig3f_source_data.xlsx
# =============================================================================

rm(list = ls())
library(readxl); library(dplyr); library(ggplot2); library(openxlsx); library(tidyr)

# ---- PATHS -------------------------------------------------------------------
OUT <- "."
stopifnot(dir.exists(OUT)); setwd(OUT)
IN  <- "fgsea_results.xlsx"
stopifnot(file.exists(IN))

# ---- STYLE -------------------------------------------------------------------
# #FFE0C0 has only 1.23:1 contrast against white, so both markers take a dark
# outline rather than a white one; the fills are unchanged.
ORANGE_LIGHT <- "#FFE0C0"   # disease state
ORANGE_DARK  <- "#FF6A00"   # treatment
EDGE <- "#8a3a00"; INK <- "#0b0b0b"; MUTED <- "#52514e"
FAMILY <- "Arial"
# Matches the other panels of Figure 3. Type that differs between panels of one
# figure is visible immediately once they are assembled, so this is set to the
# same value rather than to whatever suits this panel alone. The figure is
# sized to accommodate it.
FONT <- 20
# GO term names are long, and at this type size a single line of one of them is
# wider than the plotting area. Wrapping them keeps the type large and the
# figure narrow; shrinking the font to fit would instead make this panel
# disagree with the other three.
# 32 keeps every short label on a single line. Two-line labels crowd the rows
# and start to collide with their neighbours, so prefer widening the figure or
# shortening the label over wrapping it.
LABEL_WRAP <- 34

# GO term names are written for machines, not for figures, and at this type
# size one of them is wider than the plotting area. Each term therefore gets a
# short display name here. The full GO identifier stays in the source data
# workbook, so nothing is lost and the renaming is auditable.
# A term with no entry falls back to its wrapped GO name.
SHORT_LABELS <- c(
  GOBP_PROTEIN_LOCALIZATION_TO_PRESYNAPSE                            = "Presynaptic protein targeting",
  GOBP_DENDRITE_MORPHOGENESIS                                        = "Dendrite morphogenesis",
  GOBP_NEURON_PROJECTION_GUIDANCE                                    = "Neuron projection guidance",
  GOBP_SYNAPTIC_MEMBRANE_ADHESION                                    = "Synaptic membrane adhesion",
  GOBP_SYNAPSE_ASSEMBLY                                              = "Synapse assembly",
  GOBP_REGULATION_OF_TRANS_SYNAPTIC_SIGNALING                        = "Trans-synaptic signaling",
  GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT                      = "Oxidative phosphorylation",
  GOBP_MITOCHONDRIAL_TRANSLATION                                     = "Mitochondrial translation",
  GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_MITOCHONDRIAL_MEMBRANE = "Mitochondrial protein import",
  GOBP_RIBOSOMAL_SMALL_SUBUNIT_BIOGENESIS                            = "Small subunit biogenesis",
  GOBP_POSITIVE_REGULATION_OF_DNA_TEMPLATED_TRANSCRIPTION_INITIATION = "Transcription initiation",
  GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM                        = "Vascular process",
  GOBP_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS             = "Sprouting angiogenesis",
  GOBP_REGULATION_OF_ENDOTHELIAL_CELL_DIFFERENTIATION                = "Endothelial differentiation",
  # terms that entered with the limma-based ranking
  GOBP_POSITIVE_REGULATION_OF_EXCITATORY_POSTSYNAPTIC_POTENTIAL      = "Excitatory postsynaptic potential",
  GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC                           = "Glutamatergic transmission",
  GOBP_DENDRITE_DEVELOPMENT                                          = "Dendrite development",
  GOBP_NEGATIVE_REGULATION_OF_VASCULATURE_DEVELOPMENT                = "Vasculature development",
  GOBP_VASCULAR_TRANSPORT                                            = "Vascular transport",
  GOBP_REGULATION_OF_VASCULATURE_DEVELOPMENT                         = "Vasculature development",
  GOBP_SYNAPSE_ORGANIZATION                                          = "Synapse organization",
  GOBP_REGULATION_OF_SYNAPSE_STRUCTURE_OR_ACTIVITY                   = "Synapse structure & activity",
  GOBP_MODULATION_OF_CHEMICAL_SYNAPTIC_TRANSMISSION                  = "Synaptic transmission",
  GOBP_RIBOSOME_BIOGENESIS                                           = "Ribosome biogenesis",
  GOBP_CYTOPLASMIC_TRANSLATION                                       = "Cytoplasmic translation")
# Direction is carried by the NES, so "positive regulation of" and similar
# qualifiers are dropped from the display label. The full GO name stays in the
# source data, where the qualifier can be checked.
BLOCK_WRAP <- 40   # block headings sit above each group, so they rarely wrap
# Height is set generously: panel height follows the number of terms, not the
# height of their labels, so a label that wraps onto two lines collides with its
# neighbour in a tight figure.
FIG_W <- 13.5; FIG_H <- 9.5
MAX_OVERLAP <- 0.5   # terms sharing more of their leading edge than this are dropped

wrap_text <- function(x, n)
  vapply(x, function(s) paste(strwrap(s, width = n), collapse = "\n"), "", USE.NAMES = FALSE)

if (requireNamespace("extrafont", quietly = TRUE))
  suppressWarnings(try(extrafont::loadfonts(device = "win", quiet = TRUE), silent = TRUE))

# ---- 1. mirrored terms -------------------------------------------------------
d <- read_excel(IN, sheet = "GOBP_disease")
t <- read_excel(IN, sheet = "GOBP_treatment")
m <- inner_join(d, t, by = "pathway")

mirrored <- m %>%
  filter(padj_disease < 0.05, padj_treatment < 0.05,
         sign(NES_disease) != sign(NES_treatment)) %>%
  mutate(score = abs(NES_disease) + abs(NES_treatment)) %>%
  arrange(desc(score))
cat(sprintf("terms significant in both contrasts and opposite in sign: %d\n", nrow(mirrored)))
if (nrow(mirrored) == 0)
  stop("No mirrored terms. Show the hand-built sets from the Custom_sets sheet instead.")

# ---- 2. thematic blocks ------------------------------------------------------
BLOCKS <- list(
  "Synaptic & dendritic"     = "SYNAP|DENDRIT|POSTSYNAP|PRESYNAP|NEURON_PROJECTION|AXON|NEUROTRANSMITTER",
  "Mitochondrial"            = "MITOCHOND|OXIDATIVE_PHOSPH|ELECTRON_TRANSPORT|ATP_SYNTH",
  "Ribosome & transcription" = "RIBOSOM|RRNA|TRANSLATION|TRANSCRIPTION_INITIATION|RNA_POLYMERASE|SPLICING",
  "Vascular"                 = "VASCUL|ANGIOGEN|ENDOTHELI|BLOOD_VESSEL|NITRIC_OXIDE")
# The synaptic block carries the point of the panel, so it keeps six terms and
# the room for them comes from the figure height instead.
N_PER_BLOCK <- c("Synaptic & dendritic" = 6, "Mitochondrial" = 3,
                 "Ribosome & transcription" = 2, "Vascular" = 3)

le_set <- function(x) strsplit(ifelse(is.na(x), "", x), "/")[[1]]

pick <- function(df, k) {   # drop terms sharing more than MAX_OVERLAP of their leading edge
  kept <- list(); out <- integer(0)
  for (i in seq_len(nrow(df))) {
    g <- le_set(df$leadingEdge_treatment[i])
    if (length(g) < 3) next
    ok <- all(vapply(kept, function(kg)
      length(intersect(g, kg)) / min(length(g), length(kg)) < MAX_OVERLAP, logical(1)))
    if (ok) { kept <- c(kept, list(g)); out <- c(out, i) }
    if (length(out) == k) break
  }
  df[out, ]
}

sel <- list(); used <- character(0)
for (b in names(BLOCKS)) {
  sub <- mirrored %>% filter(grepl(BLOCKS[[b]], pathway), !pathway %in% used)
  if (b == "Mitochondrial") sub <- sub %>% filter(NES_disease > 0)   # biogenesis programme
  p <- pick(sub, N_PER_BLOCK[[b]])
  if (nrow(p)) {
    used <- c(used, p$pathway)
    sel[[b]] <- p %>% mutate(block = b)
  }
  cat(sprintf("  %-26s %3d candidates -> %d selected\n", b, nrow(sub),
              if (nrow(p)) nrow(p) else 0))
}
S <- bind_rows(sel)
if (!nrow(S)) stop("No terms selected: loosen BLOCKS or raise MAX_OVERLAP.")

# ---- 3. labels and order -----------------------------------------------------
S <- S %>% mutate(
  label = pathway |> sub("^GOBP_", "", x = _) |> gsub("_", " ", x = _) |> tolower() |>
    (\(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))))() |>
    gsub("dna templated", "DNA-templated", x = _) |>
    gsub("\\batp\\b", "ATP", x = _, ignore.case = TRUE) |>
    gsub("\\brrna\\b", "rRNA", x = _, ignore.case = TRUE),
  block = factor(block, levels = names(BLOCKS)))

S <- S %>% group_by(block) %>% arrange(block, NES_treatment) %>% ungroup()

# short display name where one exists, otherwise the wrapped GO name
short <- ifelse(S$pathway %in% names(SHORT_LABELS),
                SHORT_LABELS[S$pathway], as.character(S$label))
missing_short <- S$pathway[!S$pathway %in% names(SHORT_LABELS)]
if (length(missing_short)) {
  cat("\nno short label for these terms; the GO name is wrapped instead:\n  ",
      paste(missing_short, collapse = "\n  "), "\n",
      "  Add them to SHORT_LABELS if the panel looks crowded.\n", sep = "")
}
short <- wrap_text(short, LABEL_WRAP)

S <- S %>% mutate(
  label_plot = factor(short, levels = unique(short)),
  block_plot = factor(wrap_text(as.character(block), BLOCK_WRAP),
                      levels = unique(wrap_text(names(BLOCKS), BLOCK_WRAP))),
  label = factor(label, levels = unique(label)))
print(as.data.frame(S %>% dplyr::select(block, label, NES_disease, padj_disease,
                                        NES_treatment, padj_treatment) %>%
                      mutate(across(where(is.numeric), ~ signif(., 3)))), row.names = FALSE)

# ---- 4. panel ----------------------------------------------------------------
plot_df <- S %>%
  pivot_longer(c(NES_disease, NES_treatment), names_to = "contrast", values_to = "NES") %>%
  mutate(contrast = factor(ifelse(contrast == "NES_disease",
                                  "hTau vs nTg (disease)",
                                  "RA013915 vs vehicle (hTau)"),
                           levels = c("hTau vs nTg (disease)",
                                      "RA013915 vs vehicle (hTau)")))

p <- ggplot(plot_df, aes(x = NES, y = label_plot)) +
  geom_vline(xintercept = 0, linetype = "22", colour = MUTED, linewidth = 0.3) +
  geom_line(aes(group = label_plot), colour = "#c9c9c4", linewidth = 0.8) +
  geom_point(aes(fill = contrast), shape = 21, size = 5, colour = EDGE, stroke = 0.6) +
  scale_fill_manual(values = c(ORANGE_LIGHT, ORANGE_DARK), name = NULL) +
  # one row: at this figure width both entries fit side by side
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = "Normalized enrichment score (NES)", y = NULL) +
  theme_classic(base_size = FONT, base_family = FAMILY) +
  theme(text = element_text(colour = INK, family = FAMILY),
        axis.text = element_text(colour = INK, size = FONT, family = FAMILY),
        axis.title.x = element_text(size = FONT),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.placement = "outside", strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "bold", size = FONT),
        strip.text.y.left = element_text(angle = 0, hjust = 0, face = "bold", size = FONT),
        panel.spacing.y = unit(14, "pt"),
        legend.position = "top", legend.key.size = unit(16, "pt"),
        legend.text = element_text(size = FONT))

# Block headings above each group rather than down the left margin: a column of
# headings is as wide as its longest word and that width is lost from the
# plotting area. facet_col keeps the panel heights proportional to the number of
# terms, which plain facet_wrap does not.
if (requireNamespace("ggforce", quietly = TRUE)) {
  p <- p + ggforce::facet_col(ggplot2::vars(block_plot), scales = "free_y",
                              space = "free", strip.position = "top")
} else {
  cat("ggforce is not installed; headings stay in the left margin and the figure\n",
      "stays wider. install.packages('ggforce') for the compact layout.\n", sep = "")
  p <- p + facet_grid(block_plot ~ ., scales = "free_y", space = "free_y", switch = "y")
  FIG_W <- FIG_W + 2.2
}

tryCatch(ggsave("Fig3f_GSEA_mirror.pdf", p, width = FIG_W, height = FIG_H, device = cairo_pdf),
  error = function(e) { message("cairo_pdf failed: ", conditionMessage(e))
    ggsave("Fig3f_GSEA_mirror.pdf", p, width = FIG_W, height = FIG_H) })
ggsave("Fig3f_GSEA_mirror.png", p, width = FIG_W, height = FIG_H, dpi = 400)

# ---- 5. source data ----------------------------------------------------------
wb <- createWorkbook()
addWorksheet(wb, "Selected_terms")
# the full GO name stays beside the short display label, so the renaming in
# SHORT_LABELS can be checked against the original
writeData(wb, "Selected_terms", S %>%
          mutate(display_label = gsub("\n", " ", as.character(label_plot))) %>%
          dplyr::select(block, display_label, label, pathway, size_disease,
          NES_disease, padj_disease, NES_treatment, padj_treatment,
          leadingEdge_disease, leadingEdge_treatment))
addWorksheet(wb, "All_mirrored_terms"); writeData(wb, "All_mirrored_terms", mirrored)
saveWorkbook(wb, "Fig3f_source_data.xlsx", overwrite = TRUE)

cat("\nwritten to", normalizePath(OUT), ":\n",
    " Fig3f_GSEA_mirror.pdf\n Fig3f_GSEA_mirror.png\n Fig3f_source_data.xlsx\n")
