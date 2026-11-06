#!/usr/bin/env Rscript
# =============================================================================
# Differential expression model for the RNA-seq data
#
# THE DESIGN, AND THE TWO PROBLEMS IT CONTAINS
# Sixteen animals, four per group, were sequenced in two runs. nTg vehicle and
# nTg RA013915 are present only in run C3SVT8 and hTau RA013915 only in run
# VJ385X, so run and group are confounded. The four hTau vehicle animals were
# sequenced in BOTH runs, and those repeated libraries are what makes the run
# effect estimable: the script verifies this by checking the rank of the design
# matrix before fitting anything.
#
# The same repeated libraries create the second problem. Twenty libraries come
# from sixteen animals, so hTau vehicle contributes eight libraries from four
# mice. A model that treats them as eight independent animals underestimates
# the variance of that group and inflates the significance of both contrasts
# that use it.
#
# limma-voom with duplicateCorrelation resolves both at once: the run effect is
# a term in the design, and block = animal declares that the two libraries of
# each hTau vehicle mouse are repeated measurements of the same individual.
# Pooling the two runs instead (for example with collapseReplicates) would
# destroy the very information that makes the run effect estimable, and would
# also double the sequencing depth of that group alone.
#
# Model choice was verified against four alternatives; see
# 11_model_comparison.R. Effect sizes are near-identical across all of them
# (Spearman 0.95-0.98), so the choice affects how many genes pass a threshold,
# not what the data say.
#
# Contrasts fitted
#   disease     hTau vehicle      vs nTg vehicle
#   drug_hTau   hTau RA013915     vs hTau vehicle
#   drug_nTg    nTg RA013915      vs nTg vehicle
#   interaction (drug effect in hTau) - (drug effect in nTg)
#
# The interaction is the formal test of genotype selectivity. Comparing the
# number of significant genes in two separate contrasts is not a test, because
# the difference between "significant" and "not significant" is not itself
# significant.
#
# Output: RNAseq_20lib_results.xlsx     all contrasts, one sheet
#         GEO_raw_counts_20libraries.tsv  unfiltered counts for deposition
#         GEO_sample_sheet.tsv            library to animal mapping
# =============================================================================

rm(list = ls())
library(readr); library(limma); library(edgeR)

# ---- PATHS -------------------------------------------------------------------
# Per-run expression matrices as delivered by the sequencing provider.
# Paths are relative to the working directory: place the two provider folders
# (C3SVT8_results/ and VJ385X_results/) beside this script, or set MATRIX_DIR.
MATRIX_DIR <- "."
RUN_C3 <- file.path(MATRIX_DIR, "C3SVT8_results", "C3SVT8-expression-matrix.tsv")
RUN_VJ <- file.path(MATRIX_DIR, "VJ385X_results", "VJ385X-expression-matrix.tsv")
OUT    <- "."
stopifnot(file.exists(RUN_C3), file.exists(RUN_VJ), dir.exists(OUT)); setwd(OUT)

# Column assignment within each provider matrix. Verify against the submission
# sheet before running: everything downstream depends on it.
# Library names are built as <group>_<n><suffix>, so that the animal number
# stays next to the group and the run suffix comes last: hTau_Veh_1_runA. Every
# other script in this repository reads the count matrix with these names, so
# the order of the parts is not cosmetic.
MAP <- list(
  list(run = "C3SVT8", file = RUN_C3, cols = 1:4,  group = "nTg_Veh",       suffix = ""),
  list(run = "C3SVT8", file = RUN_C3, cols = 5:8,  group = "nTg_RA013915",  suffix = ""),
  list(run = "C3SVT8", file = RUN_C3, cols = 9:12, group = "hTau_Veh",      suffix = "_runA"),
  list(run = "VJ385X", file = RUN_VJ, cols = 1:4,  group = "hTau_Veh",      suffix = "_runB"),
  list(run = "VJ385X", file = RUN_VJ, cols = 5:8,  group = "hTau_RA013915", suffix = ""))

FDR <- 0.05

# ---- 1. count matrix ---------------------------------------------------------
read_matrix <- function(f) {
  x <- as.data.frame(read_tsv(f, show_col_types = FALSE))
  rownames(x) <- make.unique(ifelse(is.na(x$gene_name) | x$gene_name == "",
                                    x$gene_id, x$gene_name))
  x
}
mats <- lapply(unique(vapply(MAP, `[[`, "", "file")), read_matrix)
names(mats) <- unique(vapply(MAP, `[[`, "", "file"))

common <- Reduce(intersect, lapply(mats, rownames))
cat(sprintf("genes annotated in both runs: %d\n", length(common)))

blocks <- lapply(MAP, function(b) {
  x <- mats[[b$file]][common, paste0(b$run, "_", b$cols, "_count"), drop = FALSE]
  colnames(x) <- paste0(b$group, "_", seq_along(b$cols), b$suffix)
  x
})
counts <- round(as.matrix(do.call(cbind, blocks)))
cat(sprintf("count matrix: %d genes x %d libraries\n", nrow(counts), ncol(counts)))
cat("library sizes (millions): ",
    paste(sprintf("%.2f", colSums(counts) / 1e6), collapse = ", "), "\n", sep = "")

# ---- 2. sample information ---------------------------------------------------
# biological_sample is shared by the two libraries of each hTau vehicle animal.
# It is both the blocking factor for the model and the declaration that those
# libraries are not independent animals, which is what any reader of the
# deposited metadata needs in order to interpret the design.
# group, run and animal are taken from MAP rather than parsed back out of the
# column names: parsing a name that the script itself just built is a way to
# introduce a mismatch that no error message points at.
info <- data.frame(
  row.names = colnames(counts),
  group = factor(unlist(lapply(MAP, function(b) rep(b$group, length(b$cols)))),
                 levels = c("nTg_Veh", "nTg_RA013915", "hTau_Veh", "hTau_RA013915")),
  run = factor(unlist(lapply(MAP, function(b) rep(b$run, length(b$cols)))),
               levels = c("C3SVT8", "VJ385X")),
  # the two libraries of each hTau vehicle animal share this identifier
  biological_sample = factor(unlist(lapply(MAP, function(b)
    paste0(b$group, "_", seq_along(b$cols))))))
stopifnot(!any(is.na(info$group)), nlevels(info$biological_sample) == 16)
print(table(info$group, info$run))
cat(sprintf("%d libraries from %d animals\n\n",
            nrow(info), nlevels(info$biological_sample)))

# ---- 3. is the run effect estimable? -----------------------------------------
design <- model.matrix(~ run + group, info)
cat(sprintf("design ~ run + group: %d columns, rank %d -> %s\n",
            ncol(design), qr(design)$rank,
            if (qr(design)$rank == ncol(design))
              "full rank, the run effect is estimable"
            else "rank deficient, run and group cannot be separated"))
if (qr(design)$rank < ncol(design))
  stop("rank-deficient design: check the column assignment in MAP")
stopifnot(identical(colnames(design),
  c("(Intercept)", "runVJ385X", "groupnTg_RA013915",
    "grouphTau_Veh", "grouphTau_RA013915")))

CONTRASTS <- cbind(
  disease     = c(0, 0,  0,  1, 0),
  drug_hTau   = c(0, 0,  0, -1, 1),
  drug_nTg    = c(0, 0,  1,  0, 0),
  interaction = c(0, 0, -1, -1, 1))
rownames(CONTRASTS) <- colnames(design)
cat("\ncontrast matrix:\n"); print(CONTRASTS)

# ---- 4. filtering ------------------------------------------------------------
# filterByExpr with the full design: the threshold is set by the smallest group
# (four libraries), not by the total number of samples.
dge <- DGEList(counts, group = info$group)
keep <- filterByExpr(dge, design = design)
cat(sprintf("\ngenes: %d total -> %d retained\n", nrow(dge), sum(keep)))
dge <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])
cat("TMM normalization factors: ",
    paste(sprintf("%.3f", dge$samples$norm.factors), collapse = ", "), "\n", sep = "")

# ---- 5. fit ------------------------------------------------------------------
# duplicateCorrelation is estimated twice, with voom weights recomputed in
# between, as recommended in the limma user guide.
v  <- voom(dge, design, plot = FALSE)
dc <- duplicateCorrelation(v, design, block = info$biological_sample)
v  <- voom(dge, design, block = info$biological_sample,
           correlation = dc$consensus, plot = FALSE)
dc <- duplicateCorrelation(v, design, block = info$biological_sample)
cat(sprintf("\nconsensus within-animal correlation: %.3f\n", dc$consensus))
cat("  A high value confirms that the two libraries of each hTau vehicle animal\n",
    "  are repeated measurements of one individual, and that a model treating\n",
    "  them as separate animals would overstate the sample size.\n", sep = "")

fit <- lmFit(v, design, block = info$biological_sample, correlation = dc$consensus)
fit <- eBayes(contrasts.fit(fit, CONTRASTS))
cat(sprintf("residual degrees of freedom: %d\n", fit$df.residual[1]))

res <- lapply(colnames(CONTRASTS), function(cn)
  topTable(fit, coef = cn, number = Inf, sort.by = "none"))
names(res) <- colnames(CONTRASTS)

# ---- 6. summary --------------------------------------------------------------
cat(sprintf("\ndifferentially expressed genes (FDR < %g):\n", FDR))
summary_tab <- do.call(rbind, lapply(names(res), function(cn) {
  r <- res[[cn]]
  data.frame(contrast = cn,
             n  = sum(r$adj.P.Val < FDR),
             up = sum(r$adj.P.Val < FDR & r$logFC > 0),
             dn = sum(r$adj.P.Val < FDR & r$logFC < 0),
             nominal_p05 = sum(r$P.Value < 0.05),
             expected_by_chance = round(0.05 * nrow(r)))
}))
print(summary_tab, row.names = FALSE)
cat("\nComparing 'nominal_p05' with 'expected_by_chance' shows whether a contrast\n",
    "carries a broad sub-threshold signal or a concentrated one.\n", sep = "")

for (cn in names(res)) {
  top <- head(rownames(res[[cn]])[order(res[[cn]]$adj.P.Val)],
              min(10, sum(res[[cn]]$adj.P.Val < FDR)))
  if (length(top)) cat(sprintf("  %-12s %s\n", cn, paste(top, collapse = ", ")))
}

# ---- 7. results table --------------------------------------------------------
# "stat" is the moderated t. It is kept because it, not the fold change, is the
# statistic that should rank genes for gene-set enrichment: it carries both the
# direction and the reliability of each estimate, so a large fold change
# measured with little confidence does not dominate the ranking.
as_block <- function(r, suffix) {
  o <- data.frame(gene = rownames(r), baseMean = 2^r$AveExpr,
                  log2FoldChange = r$logFC, stat = r$t,
                  pvalue = r$P.Value, padj = r$adj.P.Val)
  names(o)[-1] <- paste0(names(o)[-1], "_", suffix)
  o
}
merged <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE),
                 lapply(names(res), function(cn) as_block(res[[cn]], cn)))

# aliases kept for the figure scripts, which refer to the treatment contrast
for (col in c("log2FoldChange", "pvalue", "padj", "baseMean"))
  merged[[paste0(col, "_treatment")]] <- merged[[paste0(col, "_drug_hTau")]]

write.csv(merged, "RNAseq_20lib_results.csv", row.names = FALSE)
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(list(Full_Merged_Results = merged,
                           DEG_summary         = summary_tab),
                      "RNAseq_20lib_results.xlsx")
  cat("\nwritten: RNAseq_20lib_results.xlsx\n")
} else {
  cat("\nwritexl is not installed; only the CSV was written.\n",
      "install.packages('writexl') to produce the workbook read by the figure scripts.\n")
}

# ---- 8. the fitted object, for the gene-set scripts --------------------------
# fgsea works from the ranked statistics alone, but the correlation-aware check
# (fry) needs the expression values, the design and the blocking structure.
# Saving them here keeps that check on exactly the model reported above instead
# of on a second fit that could drift out of step with it.
saveRDS(list(voom = v, design = design, block = info$biological_sample,
             correlation = dc$consensus, contrasts = CONTRASTS, info = info),
        "RNAseq_20lib_model.rds")
cat("written: RNAseq_20lib_model.rds\n")

# ---- 9. files for deposition -------------------------------------------------
# Raw, unfiltered counts: the archive should hold every library generated, with
# any exclusion declared in the methods rather than by omission.
write_tsv(data.frame(gene_id = rownames(counts), gene_name = rownames(counts),
                     as.data.frame(counts), check.names = FALSE),
          "GEO_raw_counts_20libraries.tsv")
write_tsv(data.frame(library = rownames(info), info,
                     note = ifelse(duplicated(info$biological_sample) |
                                   duplicated(info$biological_sample, fromLast = TRUE),
                                   "same RNA sequenced in both runs", "")),
          "GEO_sample_sheet.tsv")

cat("\nwritten to", normalizePath(OUT), ":\n",
    " RNAseq_20lib_results.xlsx / .csv\n",
    " GEO_raw_counts_20libraries.tsv\n GEO_sample_sheet.tsv\n")
cat("\nRun 10_qc_batch_and_technical_replicates_FigS1.R next for the supplementary QC.\n")
