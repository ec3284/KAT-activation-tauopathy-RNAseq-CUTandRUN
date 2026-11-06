#!/usr/bin/env Rscript
# =============================================================================
# Model comparison for the RNA-seq design
#
# WHY THIS SCRIPT EXISTS
# The design is partly confounded: nTg Vehicle and nTg RA013915 were sequenced
# only in run C3SVT8 and hTau RA013915 only in run VJ385X, and the four hTau
# Vehicle animals were sequenced in both runs. More than one defensible model
# can be fitted to such a design, and a reader is entitled to ask whether the
# conclusions depend on which one was chosen.
#
# Five specifications are fitted here:
#
#   C_limma_dupcor    the model reported in the manuscript. All 20 libraries,
#                     ~ run + group, with the animal declared through
#                     duplicateCorrelation so that the hTau Vehicle group
#                     contributes four animals rather than eight libraries.
#
#   A_DESeq2_pseudo   all 20 libraries, ~ run + group, DESeq2, with the
#                     repeated libraries treated as independent animals.
#   Ap_limma_pseudo   the same, under limma-voom and without blocking.
#                     Both overstate the sample size of the hTau Vehicle group.
#
#   B_DESeq2_entro    each contrast fitted on its own within-run subset, so
#                     that no run term is needed and the repeated libraries
#                     are split rather than modelled: disease and drug_nTg in
#                     C3SVT8, drug_hTau in VJ385X.
#   Bp_limma_entro    the same subsets, under limma-voom.
#
# Two things are reported. The number of genes passing the threshold, which is
# what differs most between specifications, and the correlation of the effect
# sizes themselves, which is what determines whether the specifications tell
# the same biological story.
#
# A note on reading the output: the specification used in the manuscript is
# not the one that yields the most significant genes. In both contrasts that
# carry the paper's claims it is among the most conservative of the five.
#
# Input : GEO_raw_counts_20libraries.tsv  (written by 01_)
# Output: Model_comparison_DEG_counts.csv
#         Model_comparison_effect_sizes.csv
# =============================================================================

rm(list = ls())
library(readr); library(limma); library(edgeR)
stopifnot(requireNamespace("DESeq2", quietly = TRUE))

# ---- PATHS -------------------------------------------------------------------
COUNTS <- "GEO_raw_counts_20libraries.tsv"
OUT    <- "."
stopifnot(file.exists(COUNTS), dir.exists(OUT)); setwd(OUT)

FDR <- 0.05

# ---- 1. data -----------------------------------------------------------------
d <- as.data.frame(read_tsv(COUNTS, show_col_types = FALSE))
rownames(d) <- make.unique(ifelse(is.na(d$gene_name) | d$gene_name == "",
                                  d$gene_id, d$gene_name))

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
                 levels = c("nTg_Veh", "nTg_RA013915", "hTau_Veh", "hTau_RA013915")),
  run = factor(c(rep("C3SVT8", 12), rep("VJ385X", 8)), levels = c("C3SVT8", "VJ385X")),
  animal = factor(c(paste0("nTgV", 1:4), paste0("nTgD", 1:4),
                    paste0("hTauV", 1:4), paste0("hTauV", 1:4), paste0("hTauD", 1:4))))
cat(sprintf("%d libraries from %d animals\n\n", nrow(info), nlevels(info$animal)))

CONTRASTS <- c("disease", "drug_hTau", "drug_nTg")

# which libraries each contrast uses when fitted within a single run
SUBSET <- list(
  disease   = c(LIB$nTg_Veh,       LIB$hTau_Veh_runA),
  drug_hTau = c(LIB$hTau_Veh_runB, LIB$hTau_RA013915),
  drug_nTg  = c(LIB$nTg_Veh,       LIB$nTg_RA013915))

# ---- 2. the fitting routines -------------------------------------------------
# Each returns a named list of data frames, one per contrast, with log2FoldChange
# and padj. Genes are filtered with the design actually used for that fit.

fit_limma_full <- function(block) {
  design <- model.matrix(~ run + group, info)
  stopifnot(qr(design)$rank == ncol(design))
  dge <- DGEList(counts, group = info$group)
  dge <- calcNormFactors(dge[filterByExpr(dge, design = design), , keep.lib.sizes = FALSE])
  if (block) {
    v  <- voom(dge, design, plot = FALSE)
    dc <- duplicateCorrelation(v, design, block = info$animal)
    v  <- voom(dge, design, block = info$animal, correlation = dc$consensus, plot = FALSE)
    dc <- duplicateCorrelation(v, design, block = info$animal)
    cat(sprintf("  consensus within-animal correlation: %.3f\n", dc$consensus))
    fit <- lmFit(v, design, block = info$animal, correlation = dc$consensus)
  } else {
    v   <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
  }
  cm <- cbind(disease     = c(0, 0,  0,  1, 0),
              drug_hTau   = c(0, 0,  0, -1, 1),
              drug_nTg    = c(0, 0,  1,  0, 0))
  rownames(cm) <- colnames(design)
  fit <- eBayes(contrasts.fit(fit, cm))
  lapply(setNames(CONTRASTS, CONTRASTS), function(cn) {
    tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none")
    data.frame(gene = rownames(tt), log2FoldChange = tt$logFC, padj = tt$adj.P.Val)
  })
}

fit_deseq_full <- function() {
  design <- model.matrix(~ run + group, info)
  dge <- DGEList(counts, group = info$group)
  keep <- filterByExpr(dge, design = design)
  dds <- DESeq2::DESeqDataSetFromMatrix(counts[keep, ], info, ~ run + group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  pairs <- list(disease   = c("group", "hTau_Veh",      "nTg_Veh"),
                drug_hTau = c("group", "hTau_RA013915", "hTau_Veh"),
                drug_nTg  = c("group", "nTg_RA013915",  "nTg_Veh"))
  lapply(pairs, function(p) {
    r <- as.data.frame(DESeq2::results(dds, contrast = p))
    data.frame(gene = rownames(r), log2FoldChange = r$log2FoldChange, padj = r$padj)
  })
}

fit_within_run <- function(engine) {
  lapply(setNames(CONTRASTS, CONTRASTS), function(cn) {
    cols <- SUBSET[[cn]]
    sub  <- droplevels(info[cols, ])
    stopifnot(nlevels(sub$run) == 1, nlevels(sub$group) == 2)
    des  <- model.matrix(~ group, sub)
    dge  <- DGEList(counts[, cols], group = sub$group)
    keep <- filterByExpr(dge, design = des)
    if (engine == "limma") {
      dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
      fit <- eBayes(lmFit(voom(dge, des, plot = FALSE), des))
      tt  <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
      out <- data.frame(gene = rownames(tt), log2FoldChange = tt$logFC, padj = tt$adj.P.Val)
    } else {
      dds <- DESeq2::DESeqDataSetFromMatrix(counts[keep, cols], sub, ~ group)
      dds <- DESeq2::DESeq(dds, quiet = TRUE)
      r   <- as.data.frame(DESeq2::results(dds, name = DESeq2::resultsNames(dds)[2]))
      out <- data.frame(gene = rownames(r), log2FoldChange = r$log2FoldChange, padj = r$padj)
    }
    # droplevels keeps the order of the original factor, so the reference level
    # is nTg_Veh for disease and drug_nTg and hTau_Veh for drug_hTau, which is
    # the same direction as the contrasts fitted on the full design
    out
  })
}

# ---- 3. fit all five ---------------------------------------------------------
cat("C_limma_dupcor\n");  C  <- fit_limma_full(block = TRUE)
cat("Ap_limma_pseudo\n"); Ap <- fit_limma_full(block = FALSE)
cat("A_DESeq2_pseudo\n"); A  <- fit_deseq_full()
cat("Bp_limma_entro\n");  Bp <- fit_within_run("limma")
cat("B_DESeq2_entro\n");  B  <- fit_within_run("deseq")

MODELS <- list(A_DESeq2_pseudo = A, Ap_limma_pseudo = Ap,
               C_limma_dupcor  = C,
               B_DESeq2_entro  = B, Bp_limma_entro  = Bp)

# ---- 4. genes passing the threshold ------------------------------------------
counts_tab <- do.call(rbind, lapply(CONTRASTS, function(cn) {
  row <- vapply(MODELS, function(m) sum(m[[cn]]$padj < FDR, na.rm = TRUE), 0L)
  data.frame(contrasto = cn, t(row), check.names = FALSE)
}))
cat(sprintf("\ndifferentially expressed genes (FDR < %g):\n", FDR))
print(counts_tab, row.names = FALSE)
write.csv(counts_tab, "Model_comparison_DEG_counts.csv", row.names = FALSE)

# ---- 5. agreement between effect sizes ---------------------------------------
# The counts differ because the threshold is a step function. What matters for
# the conclusions is whether the estimates themselves agree, so each alternative
# is correlated against the reported model over the genes both retained.
ref <- "C_limma_dupcor"
eff <- do.call(rbind, lapply(CONTRASTS, function(cn) {
  r <- MODELS[[ref]][[cn]]
  do.call(rbind, lapply(setdiff(names(MODELS), ref), function(m) {
    a <- MODELS[[m]][[cn]]
    k <- intersect(r$gene, a$gene)
    x <- r$log2FoldChange[match(k, r$gene)]
    y <- a$log2FoldChange[match(k, a$gene)]
    ok <- is.finite(x) & is.finite(y)
    data.frame(contrasto = cn, model = m, n_genes = sum(ok),
               spearman = cor(x[ok], y[ok], method = "spearman"),
               pearson  = cor(x[ok], y[ok]))
  }))
}))
cat("\nagreement of effect sizes with the reported model:\n")
print(as.data.frame(lapply(eff, function(z) if (is.numeric(z)) round(z, 3) else z)),
      row.names = FALSE)
write.csv(eff, "Model_comparison_effect_sizes.csv", row.names = FALSE)

cat(sprintf("\nSpearman range across all comparisons: %.2f to %.2f\n",
            min(eff$spearman), max(eff$spearman)))
cat("\nwritten to", normalizePath(OUT), ":\n",
    " Model_comparison_DEG_counts.csv\n Model_comparison_effect_sizes.csv\n")
