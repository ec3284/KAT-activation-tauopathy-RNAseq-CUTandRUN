#!/usr/bin/env Rscript
# =============================================================================
# Differential H3K4ac peaks, spike-in normalized
# Source of the responding-region set used by Fig. 4c; no panel of its own
#
# This is the quantitative analysis the figure needs, and it replaces three
# things that could not support one.
#
#   The count of called peaks per condition. The four consensus sets hold
#   133,194 / 23,458 / 83,121 / 191,623 regions at equal sequencing depth
#   (8.1 to 9.3 million reads in peaks). A set that is eight times larger than
#   another is reporting how readily that library passed the caller's
#   threshold. At the Fos promoter, nTg RA013915 carries the highest read
#   count of the four groups and is the only group with no call there.
#
#   The paired test across peaks. A signed-rank test over the peaks of a gene
#   set asks whether their median change differs from zero. In this dataset
#   96.9% of all peaks rise with the compound in hTau mice, so six peaks all
#   rising has probability 0.83 by chance. The test returns a small P because
#   the genome-wide shift is present at every peak, and the P value scales
#   with how many peaks are placed in the set: 11 peaks give P = 0.004 and 514
#   peaks give P = 3e-85 for the same per-animal effect.
#
#   Peak sets selected by position. A set built by intersecting peaks with
#   promoter windows is 100% promoter by construction.
#
# What this script does instead: one common set of regions, raw counts per
# animal, spike-in size factors, and a model that knows there are three mice
# per group. Nothing here requires the BAM files. Spike-in normalization is
# applied in full: it enters DESeq2 as the size factors, which is what
# spike-in normalization is, and the twelve spike-in read counts are set
# below.
#
# The leave-one-out section is not optional. hTau_915_R1 carries 19,363
# spike-in reads against a minimum of 35,143 in the other eleven libraries,
# so it receives a 2.66-fold multiplier while its group mates receive 0.78 and
# 1.46. Across the three treated animals the size factors alone differ from
# the vehicles by 0.63 log2 before any data are read, and the global shift in
# the data is 0.79 log2. A result that disappears when that one library is
# dropped is a property of its spike-in count, not of the compound. Report
# what the leave-one-out shows either way.
#
# Input : CUTRUN_union_peak_counts_12mice.csv.gz
#           peak_id, chrom, start, end, then one raw count column per mouse,
#           over the four-way consensus union (299,069 regions, chrM included;
#           see 02_).
# Output: DiffPeaks_results.xlsx         all three contrasts
#         DiffPeaks_up_hTau_FDR05.bed    for the feature distribution panel
#         DiffPeaks_MAplot_hTau.pdf
#         DiffPeaks_sizefactors.csv
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({ library(DESeq2); library(readr); library(dplyr) })
if (!requireNamespace("writexl", quietly = TRUE))
  install.packages("writexl", repos = "https://cloud.r-project.org")

# ---- PATHS -------------------------------------------------------------------
# Run from the directory where the outputs should be written. COUNTS is the
# per-mouse count matrix over the peak union, written by 03_union_peak_counts.sh.
COUNTS <- Sys.getenv("CUTRUN_COUNTS", "CUTRUN_union_peak_counts_12mice.csv.gz")
stopifnot(file.exists(COUNTS))

# ---- SPIKE-IN ----------------------------------------------------------------
# Reads mapped to the spike-in genome, from the bowtie2 logs. These are the
# normalization: DESeq2 divides each library's counts by its size factor, so
# setting the size factors to the spike-in counts is spike-in normalization,
# the same operation DiffBind performs with normalize = DBA_NORM_SPIKEIN.
SPIKE <- c(
  WT_Veh_R1 = 74184, WT_Veh_R2 = 35143, WT_Veh_R3 = 54326,
  WT_915_R1 = 75966, WT_915_R2 = 36713, WT_915_R3 = 54923,
  hTau_Veh_R1 = 56405, hTau_Veh_R2 = 57732, hTau_Veh_R3 = 51374,
  hTau_915_R1 = 19363, hTau_915_R2 = 66207, hTau_915_R3 = 35154)

MIN_COUNT <- 120   # total reads across the twelve libraries for a peak to be tested
FDR <- 0.05

# ---- 1. counts ---------------------------------------------------------------
d <- read_csv(COUNTS, show_col_types = FALSE)
stopifnot(all(c("peak_id", "chrom", "start", "end") %in% names(d)))
regions <- d[, c("peak_id", "chrom", "start", "end")]
m <- as.matrix(d[, names(SPIKE)])          # columns reordered to match SPIKE
rownames(m) <- d$peak_id
cat(sprintf("regions: %s | libraries: %d\n", format(nrow(m), big.mark = ","), ncol(m)))
cat(sprintf("reads in peaks per library: %s\n",
            paste(format(colSums(m), big.mark = ","), collapse = ", ")))

keep <- rowSums(m) >= MIN_COUNT
cat(sprintf("regions tested (>= %d reads in total): %s\n",
            MIN_COUNT, format(sum(keep), big.mark = ",")))
m <- m[keep, ]; regions <- regions[keep, ]

group <- factor(sub("_R[0-9]+$", "", colnames(m)),
                levels = c("WT_Veh", "WT_915", "hTau_Veh", "hTau_915"))
coldata <- data.frame(row.names = colnames(m), group = group)
print(table(group))

# ---- 2. the model ------------------------------------------------------------
fit <- function(mat, cd, label) {
  dds <- DESeqDataSetFromMatrix(mat, cd, design = ~ group)
  dds <- estimateSizeFactors(dds)
  own <- sizeFactors(dds)                         # what DESeq2 would use
  sf  <- SPIKE[colnames(mat)] / mean(SPIKE[colnames(mat)])
  sizeFactors(dds) <- sf                          # spike-in normalization
  cat(sprintf("\n--- %s ---\n", label))
  print(round(rbind(`spike-in` = sf, `DESeq2 internal` = own), 3))
  # Printing both matters: where they disagree, every difference the analysis
  # reports is the disagreement. They are not interchangeable and the reader
  # of the Methods should be able to see which was used and by how much it
  # differs from the internal estimate.
  DESeq(dds, quiet = TRUE)
}

dds <- fit(m, coldata, "all twelve libraries")

get <- function(dds, a, b, name) {
  r <- as.data.frame(results(dds, contrast = c("group", b, a), alpha = FDR))
  r$peak_id <- rownames(r)
  up <- sum(r$padj < FDR & r$log2FoldChange > 0, na.rm = TRUE)
  dn <- sum(r$padj < FDR & r$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("%-26s FDR<%.2f:  up %7s   down %7s   (%.1f%% of regions)   median log2FC %+.3f\n",
              name, FDR, format(up, big.mark = ","), format(dn, big.mark = ","),
              100 * (up + dn) / nrow(r), median(r$log2FoldChange, na.rm = TRUE)))
  r %>% left_join(regions, by = "peak_id") %>%
    select(peak_id, chrom, start, end, baseMean, log2FoldChange, lfcSE, pvalue, padj)
}

cat("\n")
res_drug_hTau <- get(dds, "hTau_Veh", "hTau_915", "RA013915 in hTau")
res_drug_nTg  <- get(dds, "WT_Veh",   "WT_915",   "RA013915 in nTg")
res_disease   <- get(dds, "WT_Veh",   "hTau_Veh", "hTau vs nTg (vehicle)")

# A median log2FC far from zero is the genome-wide shift, not a set of
# responding regions. If most regions move together, the size factors and not
# the chromatin are what differ between the groups.

# ---- 3. leave-one-out --------------------------------------------------------
# The same contrast without the library whose spike-in count is an outlier.
drop <- "hTau_915_R1"
m2 <- m[, setdiff(colnames(m), drop)]
cd2 <- data.frame(row.names = colnames(m2),
                  group = factor(sub("_R[0-9]+$", "", colnames(m2)),
                                 levels = levels(group)))
dds2 <- fit(m2, cd2, paste("leave-one-out: without", drop))
cat("\n")
res_loo <- get(dds2, "hTau_Veh", "hTau_915", "RA013915 in hTau (LOO)")

shared <- res_drug_hTau$peak_id[which(res_drug_hTau$padj < FDR)]
kept   <- res_loo$peak_id[which(res_loo$padj < FDR)]
cat(sprintf("\nregions significant with all three treated mice: %s\n",
            format(length(shared), big.mark = ",")))
cat(sprintf("still significant after dropping %s: %s (%.1f%%)\n", drop,
            format(length(intersect(shared, kept)), big.mark = ","),
            100 * length(intersect(shared, kept)) / max(1, length(shared))))

# ---- 4. sensitivity: what the internal size factors give ---------------------
# Not an alternative to report in place of the spike-in result. It is the
# comparison that says how much of the result is carried by the normalization:
# if the spike-in analysis finds tens of thousands of regions and this one
# finds almost none, the difference between the groups lives in the twelve
# size factors rather than in the counts.
dds3 <- DESeq(estimateSizeFactors(
  DESeqDataSetFromMatrix(m, coldata, design = ~ group)), quiet = TRUE)
cat("\n--- sensitivity, DESeq2 internal size factors ---\n")
invisible(get(dds3, "hTau_Veh", "hTau_915", "RA013915 in hTau (internal)"))

# ---- 5. outputs --------------------------------------------------------------
sig <- res_drug_hTau %>% filter(!is.na(padj), padj < FDR, log2FoldChange > 0) %>%
  arrange(padj)
write.table(sig %>% select(chrom, start, end), "DiffPeaks_up_hTau_FDR05.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("\nDiffPeaks_up_hTau_FDR05.bed: %s regions\n", format(nrow(sig), big.mark = ",")))
# This file is written for completeness and is expected to be empty or nearly
# so: no region reaches FDR < 0.05 in any contrast of this dataset. The panel
# that asks where the responding regions sit is therefore built in
# 06_feature_distribution_responsive_Fig4C.R from a nominal threshold, named as
# such, against a coverage-matched background. Do not pass this BED file to the
# four-condition feature-distribution script as if it were a result.

pdf("DiffPeaks_MAplot_hTau.pdf", width = 7, height = 6)
plotMA(results(dds, contrast = c("group", "hTau_915", "hTau_Veh")), ylim = c(-4, 4),
       main = "RA013915 vs vehicle, hTau (spike-in normalized)")
abline(h = median(res_drug_hTau$log2FoldChange, na.rm = TRUE), col = "red", lty = 2)
dev.off()
# The dashed line is the median fold change. A cloud centred on it rather than
# on zero is a global shift.

writexl::write_xlsx(list(
  `RA013915 in hTau`      = res_drug_hTau,
  `RA013915 in nTg`       = res_drug_nTg,
  `Disease (vehicle)`     = res_disease,
  `RA013915 hTau leave1out` = res_loo),
  "DiffPeaks_results.xlsx")

write.csv(data.frame(library = names(SPIKE), spikein_reads = as.integer(SPIKE),
                     size_factor = round(SPIKE / mean(SPIKE), 4)),
          "DiffPeaks_sizefactors.csv", row.names = FALSE)

# The per-contrast tables are also written as plain CSV, because
# 06_feature_distribution_responsive_Fig4C.R reads the hTau contrast from one.
write.csv(res_drug_hTau, "DiffPeaks_RA013915_in_hTau.csv", row.names = FALSE)
write.csv(res_drug_nTg,  "DiffPeaks_RA013915_in_nTg.csv",  row.names = FALSE)
write.csv(res_disease,   "DiffPeaks_disease_vehicle.csv",  row.names = FALSE)

cat("\nwritten to", normalizePath("."), ":\n",
    " DiffPeaks_results.xlsx\n",
    " DiffPeaks_RA013915_in_hTau.csv\n DiffPeaks_RA013915_in_nTg.csv\n",
    " DiffPeaks_disease_vehicle.csv\n",
    " DiffPeaks_up_hTau_FDR05.bed\n DiffPeaks_MAplot_hTau.pdf\n",
    " DiffPeaks_sizefactors.csv\n")
