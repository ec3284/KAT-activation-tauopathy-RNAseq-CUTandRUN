#!/usr/bin/env Rscript
# =============================================================================
# Gene-set enrichment on the two contrasts
#
# Genes are ranked by the moderated t statistic rather than by fold change.
# The t carries both the direction and the reliability of each estimate, so a
# large fold change measured with little confidence does not dominate the
# ranking. With four animals per group that distinction matters.
#
# This analysis is what the single-gene results cannot provide. Only seven
# genes pass the FDR threshold in the treatment contrast, and the p-value
# distribution of that contrast is close to uniform, so there is no large
# population of genes just below the threshold. A coordinated shift can still
# exist inside defined gene sets without producing individually significant
# genes, and that is the question asked here, using every gene and no cut-off.
#
# TWO TESTS, DELIBERATELY
#   fgsea    competitive, permutes gene labels. It does not account for
#            correlation between genes within a set, so for strongly
#            co-expressed sets (OXPHOS, ribosome) its p-values are
#            anti-conservative. Used for discovery.
#   fry      rotation test, correlation-aware, and run on the fitted model
#            itself including the animal blocking. Used as the check. A term
#            that fgsea calls significant and fry does not is a term whose
#            significance came from within-set correlation.
#
# Input : RNAseq_20lib_results.xlsx and RNAseq_20lib_model.rds (from 01_)
# Output: fgsea_results.xlsx, read by 08_gsea_mirror_Fig3F.R
# =============================================================================

rm(list = ls())
library(readxl); library(fgsea); library(msigdbr); library(dplyr)
library(limma); library(openxlsx)

# ---- PATHS -------------------------------------------------------------------
IN    <- "RNAseq_20lib_results.xlsx"
MODEL <- "RNAseq_20lib_model.rds"
OUT   <- "."
stopifnot(file.exists(IN), dir.exists(OUT)); setwd(OUT)

MIN_SIZE <- 10; MAX_SIZE <- 500
SEED <- 42

res <- as.data.frame(read_excel(IN, sheet = "Full_Merged_Results"))
cat(sprintf("genes: %d\n", nrow(res)))
stopifnot(all(c("stat_disease", "stat_drug_hTau") %in% names(res)))

# ---- 1. ranked lists ---------------------------------------------------------
mk_rank <- function(df, col) {
  r <- df[[col]]; names(r) <- df$gene
  r <- r[is.finite(r)]
  sort(r[!duplicated(names(r))], decreasing = TRUE)
}
rank_dis <- mk_rank(res, "stat_disease")
rank_trt <- mk_rank(res, "stat_drug_hTau")
cat(sprintf("ranked: disease %d | treatment %d genes\n", length(rank_dis), length(rank_trt)))

# ---- 2. published gene sets --------------------------------------------------
get_sets <- function(coll, sub) {
  m <- tryCatch(msigdbr(species = "Mus musculus", collection = coll, subcollection = sub),
                error = function(e) msigdbr(species = "Mus musculus",
                                            category = coll, subcategory = sub))
  split(m$gene_symbol, m$gs_name)
}
GOBP     <- get_sets("C5", "GO:BP")
REACTOME <- get_sets("C2", "CP:REACTOME")
cat(sprintf("GO BP sets: %d | Reactome sets: %d\n", length(GOBP), length(REACTOME)))

# ---- 3. hand-built sets ------------------------------------------------------
IEG <- c("Fos","Fosb","Jun","Junb","Jund","Egr1","Egr2","Egr3","Egr4","Arc",
         "Npas4","Nr4a1","Nr4a2","Nr4a3","Dusp1","Btg2","Ier2","Ier3","Klf2",
         "Klf4","Zfp36","Atf3","Sik1","Per1","Bdnf")
SYN <- c("Dlg4","Shank2","Shank3","Syngap1","Nlgn1","Nlgn2","Nlgn3","Nrxn1",
         "Nrxn2","Nrxn3","Grin1","Grin2a","Grin2b","Gria1","Gria2","Gria3",
         "Camk2a","Camk2b","Homer1","Dlgap1","Syn1","Syn2","Syt1","Snap25",
         "Rims1","Cplx1","Cplx2","Sv2a","Bsn","Pclo","Stx1a","Vamp2","Nrgn","Ntrk2")
# Klf2, Klf4 and Zfp36 are immediate-early genes, but they are also the
# endothelial shear-stress transcription factors and they track haemoglobin
# content across samples. The set is split so that the blood-independent
# neuronal core can be read on its own.
BLOOD <- c("Klf2", "Klf4", "Zfp36")

g <- names(rank_dis)
CUSTOM <- list(
  "IEG neuronal"          = intersect(setdiff(IEG, BLOOD), g),
  "IEG vascular"          = intersect(BLOOD, g),
  "Immediate-early (all)" = intersect(IEG, g),
  "Synaptic"              = intersect(SYN, g),
  "OXPHOS subunits" = grep("^(Nduf[abcsv]|Sdh[abcd]|Uqcr|Cyc1|Cox[4-8]|Atp5)", g,
                           value = TRUE, ignore.case = TRUE),
  "OXPHOS assembly" = grep("^(Ndufaf|Sdhaf|Coa[0-9]|Cox1[0-9]|Cox2[0-9]|Uqcc|Bcs1l|Surf1|Atpaf)",
                           g, value = TRUE, ignore.case = TRUE),
  "mt translation"  = grep("^(Mrpl|Mrps)", g, value = TRUE, ignore.case = TRUE),
  "mt import"       = grep("^(Timm|Tomm|Pam16|Dnajc19)", g, value = TRUE, ignore.case = TRUE))
cat("hand-built sets: ",
    paste(sprintf("%s (%d)", names(CUSTOM), lengths(CUSTOM)), collapse = ", "), "\n", sep = "")

# ---- 4. fgsea ----------------------------------------------------------------
run_fgsea <- function(sets, ranks, label, minSize = MIN_SIZE, maxSize = MAX_SIZE) {
  set.seed(SEED)
  r <- fgsea(pathways = sets, stats = ranks,
             minSize = minSize, maxSize = maxSize, eps = 0, nproc = 1)
  r$leadingEdge <- sapply(r$leadingEdge, paste, collapse = "/")
  r <- r[order(r$padj), ]
  names(r)[-1] <- paste0(names(r)[-1], "_", label)
  as.data.frame(r)
}

cat("\n--- GO BP ---\n")
go_dis <- run_fgsea(GOBP, rank_dis, "disease")
go_trt <- run_fgsea(GOBP, rank_trt, "treatment")
cat(sprintf("significant terms (padj < 0.05): disease %d | treatment %d\n",
            sum(go_dis$padj_disease < 0.05), sum(go_trt$padj_treatment < 0.05)))
cat("\ntop 15, disease:\n")
print(head(go_dis[, c("pathway", "NES_disease", "padj_disease", "size_disease")], 15),
      row.names = FALSE)
cat("\ntop 15, treatment:\n")
print(head(go_trt[, c("pathway", "NES_treatment", "padj_treatment", "size_treatment")], 15),
      row.names = FALSE)

# the panel in 08_ needs terms significant in both and moving in opposite
# directions; counting them here says in advance how large that panel can be
mir <- inner_join(go_dis, go_trt, by = "pathway") %>%
  filter(padj_disease < 0.05, padj_treatment < 0.05,
         sign(NES_disease) != sign(NES_treatment))
cat(sprintf("\nterms significant in BOTH contrasts and opposite in sign: %d\n", nrow(mir)))
if (nrow(mir) < 10)
  cat("  Fewer than ten. The mirror panel will be small, and it may be more\n",
      "  honest to show the hand-built sets instead of a thin GO selection.\n", sep = "")

cat("\n--- Reactome ---\n")
re_dis <- run_fgsea(REACTOME, rank_dis, "disease")
re_trt <- run_fgsea(REACTOME, rank_trt, "treatment")
cat("\ntop 10, disease:\n")
print(head(re_dis[, c("pathway", "NES_disease", "padj_disease", "size_disease")], 10),
      row.names = FALSE)
cat("\ntop 10, treatment:\n")
print(head(re_trt[, c("pathway", "NES_treatment", "padj_treatment", "size_treatment")], 10),
      row.names = FALSE)

cat("\n--- hand-built sets ---\n")
cu_dis <- run_fgsea(CUSTOM, rank_dis, "disease", minSize = 3)
cu_trt <- run_fgsea(CUSTOM, rank_trt, "treatment", minSize = 3)
cu <- merge(cu_dis[, c("pathway", "size_disease", "NES_disease", "pval_disease", "padj_disease")],
            cu_trt[, c("pathway", "NES_treatment", "pval_treatment", "padj_treatment")],
            by = "pathway", all = TRUE)
print(cu, row.names = FALSE)

# ---- 5. fry: the correlation-aware check -------------------------------------
if (!file.exists(MODEL)) {
  cat("\n", MODEL, " not found; run 01_differential_expression_model.R first to\n",
      "enable the correlation-aware check. fgsea results above are unaffected.\n", sep = "")
} else {
  M <- readRDS(MODEL)
  v <- M$voom; design <- M$design
  cat("\n\n--- fry on the same model (correlation-aware, animal as block) ---\n")

  # the hand-built sets, plus any GO term the mirror panel would use
  idx <- lapply(c(CUSTOM, GOBP[intersect(mir$pathway, names(GOBP))]),
                function(s) which(rownames(v) %in% s))
  idx <- idx[lengths(idx) >= 3]
  cat(sprintf("sets tested with fry: %d\n", length(idx)))

  fry_one <- function(contrast_name) {
    f <- fry(v, index = idx, design = design,
             contrast = M$contrasts[, contrast_name],
             block = M$block, correlation = M$correlation)
    f$set <- rownames(f)
    f[, c("set", "NGenes", "Direction", "PValue", "FDR")]
  }
  fry_dis <- fry_one("disease"); fry_trt <- fry_one("drug_hTau")
  names(fry_dis)[-1] <- paste0(names(fry_dis)[-1], "_disease")
  names(fry_trt)[-1] <- paste0(names(fry_trt)[-1], "_treatment")
  fry_both <- merge(fry_dis, fry_trt, by = "set")
  fry_both <- fry_both[order(fry_both$FDR_treatment), ]
  print(head(fry_both, 25), row.names = FALSE)

  # which fgsea calls survive the correlation-aware test
  chk <- merge(
    rbind(cu, mir[, names(cu)]) %>% dplyr::select(pathway, padj_disease, padj_treatment),
    fry_both, by.x = "pathway", by.y = "set")
  if (nrow(chk)) {
    lost <- chk$pathway[chk$padj_treatment < 0.05 & chk$FDR_treatment >= 0.05]
    cat(sprintf("\nsets significant with fgsea but not with fry (treatment): %d\n", length(lost)))
    if (length(lost))
      cat("  ", paste(head(lost, 15), collapse = "\n  "), "\n",
          "  For these, the fgsea p-value came from correlation between genes\n",
          "  inside the set. Do not report them without the fry result beside them.\n", sep = "")
  }
}

# ---- 6. workbook -------------------------------------------------------------
wb <- createWorkbook()
for (nm in c("GOBP_disease", "GOBP_treatment", "Reactome_disease", "Reactome_treatment",
             "Custom_sets", "Mirrored_terms"))
  addWorksheet(wb, nm)
writeData(wb, "GOBP_disease", go_dis);        writeData(wb, "GOBP_treatment", go_trt)
writeData(wb, "Reactome_disease", re_dis);    writeData(wb, "Reactome_treatment", re_trt)
writeData(wb, "Custom_sets", cu);             writeData(wb, "Mirrored_terms", mir)
if (exists("fry_both")) { addWorksheet(wb, "fry_check"); writeData(wb, "fry_check", fry_both) }
saveWorkbook(wb, "fgsea_results.xlsx", overwrite = TRUE)

cat("\nwritten to", normalizePath(OUT), ": fgsea_results.xlsx\n")
cat("Run 08_gsea_mirror_Fig3F.R next to draw the panel.\n")
