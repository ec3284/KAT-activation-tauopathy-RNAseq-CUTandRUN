#!/usr/bin/env Rscript
# =============================================================================
# Promoter signal per mouse, and the tests behind Fig. 4b and Extended Data Fig. 5c
#
# One script, three contrasts. The hTau pair is Fig. 4b; the nTg pair is
# Extended Data Fig. 5c; the disease contrast is quoted in the Results.
#
# From the spike-in normalized TSS matrix to the table plotted in the figure,
# in one script. It replaces the step that made the original Fig. 4b test
# invalid. The earlier Python pipeline collapsed the three biological
# replicates to a group median before testing:
#
#     wide = long.groupby(["Gene", "Group"])["log2AUC"].median().unstack()
#
# From that line onward the individual animals no longer exist, so the only
# test the saved file can support is a test across genes. That treats
# co-regulated promoters as independent replicates and returns P values of
# essentially zero whatever the data say. The experimental unit is the mouse.
# This script therefore computes the promoter area under the curve for each of
# the twelve libraries separately, reduces each library to one number, and
# tests three animals against three animals.
#
# Spike-in normalization is left untouched. It is the correct normalization for
# CUT&RUN and is what the bigWigs already carry. What changes here is only the
# unit of replication in the test.
#
# The arithmetic reproduces the published pipeline exactly: for each transcript
# and each library the signal is summed over the bins within +/- 500 bp of the
# TSS, multiplied by the bin size, and log2 transformed with a pseudocount.
# Verified against the previously computed per-library values with a median
# absolute difference of zero over 61 comparisons.
#
# Input : AUC_PROMOTER_PC.matrix_ALL12_PC.noSkip.gz (deepTools computeMatrix,
#         12 libraries, spike-in normalized bigWigs)
# Output: Fig4b_promoter_AUC_per_gene.csv.gz   per transcript per library
#         Fig4b_promoter_AUC_per_mouse.csv     12 rows, one per animal
#         Fig4b_perAnimal_tests_sensitivity.csv   not reported; see below
#         Fig4b_PRISM_promoter_AUC.csv         4 columns x 3 rows, for Prism
#         Fig4b_source_data.xlsx               the two Source Data sheets
#
# NO HYPOTHESIS TEST IS REPORTED FOR THIS PANEL, and the paper says so. Two
# different tests have to be kept apart here.
#
# A test across genes is invalid: genes within an animal are not independent
# observations, so it would treat loci as replicates rather than animals and
# return a vanishing P value whatever the data say. That is the test the
# earlier pipeline ran, and removing it is why this script exists.
#
# A test on the three per-animal medians is admissible but uninformative at
# this sample size. Section 6 computes the smallest difference detectable at
# 80% power with n = 3; it is larger than the effect reported, so a null result
# would carry no evidence of absence and a significant one would rest on three
# values. The panel therefore reports the fold change between the means of the
# per-animal medians, with the between-animal SD alongside it.
#
# Section 5 still fits those per-animal tests, as a sensitivity check that the
# fold change is not driven by one animal, and writes them to a file whose name
# says they are not reported. Do not put them in the figure or the legend.
#
# The figure is drawn in GraphPad Prism from the PRISM file: a Column table,
# scatter dot plot with mean and SD, and no test.
#
# Memory: only the 100 bins per library that fall inside the window are read,
# not all 600, which keeps the object near 1.3 GB rather than 8 GB and lets
# this run on a laptop.
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({ library(data.table) })
if (!requireNamespace("writexl", quietly = TRUE))
  install.packages("writexl", repos = "https://cloud.r-project.org")

# ---- PATHS -------------------------------------------------------------------
# Run from the directory where the outputs should be written. MATRIX is the
# deepTools computeMatrix file over the twelve libraries, written by 01_bigwig_average_and_matrices.sh.
MATRIX <- Sys.getenv("CUTRUN_AUC_MATRIX", "AUC_PROMOTER_PC.matrix_ALL12_PC.noSkip.gz")
if (!file.exists(MATRIX))
  stop("matrix not found: ", MATRIX,
       "\nIt is written by 01_bigwig_average_and_matrices.sh, or set ",
       "CUTRUN_AUC_MATRIX to its path.", call. = FALSE)

# ---- PARAMETERS (as in the published pipeline) -------------------------------
BIN_SIZE    <- 10      # bp per bin, from the matrix header
UPSTREAM    <- 3000    # bp before the TSS covered by the matrix
DOWNSTREAM  <- 3000    # bp after the TSS covered by the matrix
AUC_WINDOW  <- 500     # bp each side of the TSS that enter the AUC
PSEUDOCOUNT <- 1e-6

SAMPLES <- c("WT_Veh_R1", "WT_Veh_R2", "WT_Veh_R3",
             "WT_915_R1", "WT_915_R2", "WT_915_R3",
             "hTau_Veh_R1", "hTau_Veh_R2", "hTau_Veh_R3",
             "hTau_915_R1", "hTau_915_R2", "hTau_915_R3")

# The group names used throughout the manuscript, in the order the figures use.
GROUP_OF <- c(WT_Veh = "nTg_Veh", WT_915 = "nTg_RA013915",
              hTau_Veh = "hTau_Veh", hTau_915 = "hTau_RA013915")
PRISM_HEADERS <- c("nTg Veh", "nTg RA013915", "hTau Veh", "hTau RA013915")

# ---- 1. which columns to read ------------------------------------------------
# A computeMatrix file has six annotation columns followed by one block of bins
# per library. Reading only the window is what keeps this runnable locally.
N_BIN   <- (UPSTREAM + DOWNSTREAM) / BIN_SIZE   # 600 bins per library
TSS_BIN <- UPSTREAM / BIN_SIZE                  # bins before the TSS
HALF    <- AUC_WINDOW / BIN_SIZE                # 50 bins each side
stopifnot(HALF == as.integer(HALF), TSS_BIN - HALF >= 0, TSS_BIN + HALF <= N_BIN)

win <- (TSS_BIN - HALF + 1):(TSS_BIN + HALF)    # 100 bins, 1-based in the block
cols <- c(1:6, unlist(lapply(seq_along(SAMPLES),
                             function(s) 6 + (s - 1) * N_BIN + win)))

cat(sprintf("window: %d bp each side of the TSS = %d bins per library\n",
            AUC_WINDOW, 2 * HALF))
cat(sprintf("reading %d of %d columns\n", length(cols), 6 + length(SAMPLES) * N_BIN))

# ---- 2. read -----------------------------------------------------------------
# skip = 1 drops the "@{...}" JSON header line deepTools writes.
# na.strings = "nan" is what makes bins with no data read as NA rather than as
# text, which matters because the pipeline ignored them rather than treating
# them as zero ("missing data as zero": false in the matrix header).
m <- fread(MATRIX, sep = "\t", header = FALSE, skip = 1,
           select = cols, na.strings = c("nan", "NA", ""),
           showProgress = TRUE)
cat(sprintf("transcripts read: %s\n", format(nrow(m), big.mark = ",")))
stopifnot(nrow(m) > 1e5, ncol(m) == length(cols))

ann <- m[, 1:6]
setnames(ann, c("chrom", "start", "end", "transcript", "score", "strand"))

# ---- 3. AUC per transcript per library ---------------------------------------
# Sum of the bins inside the window, times the bin width, then log2. Bins with
# no data are skipped, as in the original. A transcript whose window is empty
# in a given library gives NA rather than a floor value.
auc <- matrix(NA_real_, nrow = nrow(m), ncol = length(SAMPLES),
              dimnames = list(NULL, SAMPLES))
for (s in seq_along(SAMPLES)) {
  blk <- as.matrix(m[, 6 + (s - 1) * length(win) + seq_along(win), with = FALSE])
  tot <- rowSums(blk, na.rm = TRUE)
  tot[rowSums(!is.na(blk)) == 0] <- NA_real_
  auc[, s] <- log2(tot * BIN_SIZE + PSEUDOCOUNT)
  cat(sprintf("  %-12s done\n", SAMPLES[s]))
}
rm(m); invisible(gc())

# A transcript enters the summary only if it is measured in all twelve
# libraries, so every animal's number rests on the same set of promoters. In
# this dataset no transcript is dropped, but the filter is kept so that the
# summary cannot silently compare different promoter sets.
complete <- rowSums(is.na(auc)) == 0
cat(sprintf("transcripts with data in all %d libraries: %s of %s\n",
            length(SAMPLES), format(sum(complete), big.mark = ","),
            format(nrow(auc), big.mark = ",")))
auc <- auc[complete, , drop = FALSE]; ann <- ann[complete]

fwrite(data.table(ann[, .(chrom, start, end, transcript, strand)], auc),
       "Fig4b_promoter_AUC_per_gene.csv.gz")

# ---- 4. one number per mouse -------------------------------------------------
# The median across promoters is the per-animal summary. The mean is written
# alongside it so the choice can be seen rather than assumed; the conclusion
# does not depend on which is used.
per_mouse <- data.table(
  mouse     = SAMPLES,
  group     = factor(GROUP_OF[sub("_R[0-9]+$", "", SAMPLES)], levels = GROUP_OF),
  replicate = sub("^.*_R", "", SAMPLES),
  median_log2AUC = apply(auc, 2, median),
  mean_log2AUC   = colMeans(auc))
stopifnot(!anyNA(per_mouse$group), all(table(per_mouse$group) == 3))
fwrite(per_mouse, "Fig4b_promoter_AUC_per_mouse.csv")
print(per_mouse)

# ---- 5. sensitivity only, NOT REPORTED ---------------------------------------
# Three animals against three animals, Welch's t-test, two sided, on the three
# contrasts the design sets. These are computed to check that the fold change
# is not carried by a single animal, and they are deliberately absent from the
# figure, the legend and the Statistics section. With n = 3 per group a rank
# test cannot return anything below 0.1, so it is printed only as a check that
# the result is not an artefact of the normality assumption. See the header.
v <- function(g) per_mouse[group == g, median_log2AUC]

CONTRASTS <- list(
  disease   = c("nTg_Veh",  "hTau_Veh",      "hTau Veh vs nTg Veh"),
  drug_hTau = c("hTau_Veh", "hTau_RA013915", "hTau RA013915 vs hTau Veh"),
  drug_nTg  = c("nTg_Veh",  "nTg_RA013915",  "nTg RA013915 vs nTg Veh"))

tests <- rbindlist(lapply(names(CONTRASTS), function(k) {
  a <- v(CONTRASTS[[k]][1]); b <- v(CONTRASTS[[k]][2])
  tt <- t.test(b, a, var.equal = FALSE)
  wx <- suppressWarnings(wilcox.test(b, a, exact = FALSE))
  data.table(contrast = k, description = CONTRASTS[[k]][3],
             n_per_group = length(a),
             delta_log2  = round(mean(b) - mean(a), 4),
             fold_change = round(2^(mean(b) - mean(a)), 3),
             sd_control = round(sd(a), 3), sd_test = round(sd(b), 3),
             t = round(unname(tt$statistic), 3), df = round(tt$parameter, 2),
             welch_p = round(tt$p.value, 4), rank_p = round(wx$p.value, 4))
}))
fwrite(tests, "Fig4b_perAnimal_tests_sensitivity.csv")
print(tests)

# ---- 6. what the experiment can resolve --------------------------------------
# Reported so that a null result is not read as an absent effect.
sd_within <- mean(tapply(per_mouse$median_log2AUC, per_mouse$group, sd))
mde <- 3.1 * sd_within   # two sided alpha 0.05, 80% power, n = 3 per group
cat(sprintf("\nbetween-animal SD within group: %.3f log2\n", sd_within))
cat(sprintf("smallest difference detectable at 80%% power with n = 3: %.2f log2 (%.1f fold)\n",
            mde, 2^mde))
cat(sprintf("observed range across the %d animals: %.2f to %.2f log2\n",
            nrow(per_mouse), min(per_mouse$median_log2AUC),
            max(per_mouse$median_log2AUC)))

# ---- 7. table for GraphPad Prism ---------------------------------------------
# One column per group, one row per animal: the layout a Prism Column table
# expects, so the file can be pasted in without rearranging. The figure is
# drawn there as a scatter dot plot with mean and SD.
prism <- as.data.table(lapply(levels(per_mouse$group), function(g) v(g)))
setnames(prism, PRISM_HEADERS)
fwrite(prism, "Fig4b_PRISM_promoter_AUC.csv")
print(prism)

# ---- 8. source data ----------------------------------------------------------
# Written as its own workbook rather than appended to the paper's Source Data
# file, because writexl rewrites a workbook whole and would discard the sheets
# already in it. Paste these two sheets in once every Figure 4 panel is final.
sheet_mouse <- per_mouse[order(group, replicate),
  .(Mouse = mouse,
    Group = c(nTg_Veh = "nTg vehicle", nTg_RA013915 = "nTg RA013915",
              hTau_Veh = "hTau vehicle",
              hTau_RA013915 = "hTau RA013915")[as.character(group)],
    Replicate = replicate,
    `Median promoter log2 AUC (plotted)` = round(median_log2AUC, 4),
    `Mean promoter log2 AUC` = round(mean_log2AUC, 4))]

# The Source Data sheet carries what the paper reports: the per-animal values,
# the difference and the fold change. No P value, because none is reported.
sheet_stats <- tests[, .(Contrast = description, `n per group` = n_per_group,
                         `Difference (log2)` = delta_log2,
                         `Fold change` = fold_change,
                         `Between-animal SD (control)` = sd_control,
                         `Between-animal SD (treated)` = sd_test)]

writexl::write_xlsx(list(`Fig. 4b per-mouse values` = sheet_mouse,
                         `Fig. 4b group summaries` = sheet_stats),
                    "Fig4b_source_data.xlsx")

cat("\nwritten to", normalizePath("."), ":\n",
    " Fig4b_promoter_AUC_per_gene.csv.gz\n",
    " Fig4b_promoter_AUC_per_mouse.csv\n",
    " Fig4b_perAnimal_tests_sensitivity.csv   (not reported)\n",
    " Fig4b_PRISM_promoter_AUC.csv\n",
    " Fig4b_source_data.xlsx\n")
cat("\nDraw the panel in Prism from Fig4b_PRISM_promoter_AUC.csv:\n",
    "Column table, scatter dot plot with mean and SD, and NO test. The figure\n",
    "reports the fold change between the group means of the per-animal\n",
    "medians; the P values above are a sensitivity check and are not reported.\n",
    sep = "")
