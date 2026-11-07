#!/usr/bin/env Rscript
# =============================================================================
# Promoter-proximal fraction of H3K4ac peaks, one value per mouse
# Extended Data Fig. 5d, the per-mouse version of that panel
#
# Why this script exists. The feature-distribution panel is built from the four
# consensus peak sets, and a consensus set is an aggregate: one set per group,
# built by requiring a peak in 2 of 3 replicates. A chi-squared or Fisher test
# between two consensus sets therefore has the peak as its unit and n = 1 mouse
# per group. With 133,194 against 23,458 regions such a test returns a vanishing
# P value that measures how many peaks were called, not how the mark is
# distributed.
#
# Here each replicate is annotated on its own, which gives one number per
# mouse, the proportion of that animal's peaks that fall in a promoter window.
# Groups are then compared on three values each, with the mouse as the
# experimental unit. No P value from this script is reported: with three mice
# per group the comparison is underpowered, and the panel gives the composition
# itself. The contrasts are written out as a sensitivity check and are
# deliberately absent from the figure, the legend and the Statistics section.
#
# The panel itself stays as it is, drawn from the consensus sets, because that
# is the set the rest of the figure uses. This script is the per-animal check
# behind it, not a replacement for it.
#
# Input : per-replicate SEACR peak BED files (nf-core/cutandrun, the called
#         peaks directory, before the consensus step)
# Output: FigED5d_PromoterFraction_per_mouse.csv  one row per library
#         FigED5d_PromoterFraction_tests.csv      the contrasts
#         FigED5d_PromoterFraction_PRISM.csv      wide, for Prism
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(ChIPseeker); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db); library(GenomicRanges); library(rtracklayer)
})
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# ---- PATHS -------------------------------------------------------------------
# Run this script from the directory where the figure outputs should be written.
# The two locations below are read, never written to. Either edit them or set
# the environment variables CUTRUN_PEAKS and CUTRUN_BIGWIG.
#   CUTRUN_PEAKS   the consensus and per-replicate SEACR peak BED files
#   CUTRUN_BIGWIG  the spike-in normalized bigWig tracks
# REP_DIR must hold the PER-REPLICATE SEACR peak files, one per library, not
# the four consensus sets. In an nf-core/cutandrun output they are under
# 03_peak_calling/04_called_peaks. Files whose name contains "consensus" or
# "IgG" are dropped below, and the classification is printed before anything is
# computed, so a wrong directory is visible immediately.
PEAKS  <- Sys.getenv("CUTRUN_PEAKS",  "peaks_and_consensus_peaks")
REP_DIR <- Sys.getenv("CUTRUN_REPLICATE_PEAKS", PEAKS)
stopifnot(dir.exists(REP_DIR))

TSS_REGION <- c(-3000, 3000)   # same window as the panel
# The two definitions reported. The panel's bars use the full promoter window;
# the 1 kb figure is the one quoted in the Results, so both are computed and
# the legend can state which it means.

# ---- 1. find the per-replicate files ----------------------------------------
# Every .bed in REP_DIR is listed, then files whose name contains "consensus"
# or "IgG" are dropped, and the rest are classified by genotype, treatment and
# replicate from the file name. Adjust the three patterns if your names differ;
# the classification is printed before anything is computed, so a wrong pattern
# is visible immediately rather than silently producing the wrong groups.
all_bed <- list.files(REP_DIR, pattern = "\\.bed$", full.names = TRUE)
all_bed <- all_bed[!grepl("consensus|IgG", basename(all_bed), ignore.case = TRUE)]
if (!length(all_bed))
  stop("no per-replicate .bed files in ", REP_DIR, call. = FALSE)

bn       <- basename(all_bed)
genotype <- ifelse(grepl("hTau",         bn, ignore.case = TRUE), "hTau",
            ifelse(grepl("nTg|^WT|_WT",  bn, ignore.case = TRUE), "nTg", NA))
treat    <- ifelse(grepl("RA013915|_915|-915", bn, ignore.case = TRUE), "RA013915",
            ifelse(grepl("Veh",                bn, ignore.case = TRUE), "Vehicle", NA))
rep_id   <- sub(".*?[Rr]ep?([0-9]+).*", "\\1", bn)
rep_id[rep_id == bn] <- NA

info <- data.frame(file = all_bed, name = bn, genotype, treat,
                   replicate = rep_id, stringsAsFactors = FALSE)
info$group <- ifelse(is.na(info$genotype) | is.na(info$treat), NA,
                     paste(info$genotype, info$treat))

cat("files found in", REP_DIR, ":\n")
print(info[, c("name", "group", "replicate")], row.names = FALSE)

bad <- is.na(info$group) | is.na(info$replicate)
if (any(bad)) {
  cat("\nnot classified, and therefore not used:\n")
  print(info$name[bad], quote = FALSE)
  info <- info[!bad, ]
}
tb <- table(info$group)
cat("\nlibraries per group:\n"); print(tb)
if (any(tb != 3))
  warning("not three libraries in every group. Check the patterns above before ",
          "reading the tests: the n of each group is what the test uses.",
          call. = FALSE)
if (!nrow(info)) stop("nothing left to annotate", call. = FALSE)

# ---- 2. one value per mouse --------------------------------------------------
share <- function(a, pattern) {
  s <- a@annoStat
  sum(s$Frequency[grepl(pattern, s$Feature)])
}

res <- do.call(rbind, lapply(seq_len(nrow(info)), function(i) {
  a <- annotatePeak(info$file[i], TxDb = txdb, tssRegion = TSS_REGION,
                    annoDb = "org.Mm.eg.db", verbose = FALSE)
  out <- data.frame(
    library    = info$name[i],
    group      = info$group[i],
    replicate  = info$replicate[i],
    n_peaks    = length(a@anno),
    prom_1kb   = share(a, "^Promoter \\(<=1kb\\)"),
    prom_3kb   = share(a, "^Promoter"),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-40s %8s peaks   <=1kb %5.2f%%   <=3kb %5.2f%%\n",
              out$library, format(out$n_peaks, big.mark = ","),
              out$prom_1kb, out$prom_3kb))
  out
}))
cat("\n")

res$group <- factor(res$group,
  levels = c("nTg Vehicle", "nTg RA013915", "hTau Vehicle", "hTau RA013915"))
res <- res[order(res$group, res$replicate), ]
write.csv(res, "FigED5d_PromoterFraction_per_mouse.csv", row.names = FALSE)

cat("per-group means (one value per mouse, n = number of mice):\n")
for (g in levels(res$group)) {
  x <- res[res$group == g, ]
  if (!nrow(x)) next
  cat(sprintf("  %-14s n = %d   <=1kb %5.2f%% +/- %4.2f   <=3kb %5.2f%% +/- %4.2f   peaks %s\n",
              g, nrow(x), mean(x$prom_1kb), sd(x$prom_1kb),
              mean(x$prom_3kb), sd(x$prom_3kb),
              paste(format(x$n_peaks, big.mark = ","), collapse = ", ")))
}

# ---- 3. the contrasts --------------------------------------------------------
# The test is a pooled two-sample t test, chosen in advance and not after
# looking at the variances. Welch is printed next to it as a sensitivity check
# only: selecting between the two on the basis of an F test on three values per
# group distorts the error rate of whichever test is then reported, and with
# df = 2,2 the F test has a critical value of 39, so it almost never rejects.
# If the two columns disagree, say so rather than picking the smaller one.
contrast <- function(a, b, measure) {
  x <- res[[measure]][res$group == a]; y <- res[[measure]][res$group == b]
  if (length(x) < 2 || length(y) < 2) return(NULL)
  tp <- t.test(y, x, var.equal = TRUE)
  tw <- t.test(y, x, var.equal = FALSE)
  data.frame(measure = measure, reference = a, comparison = b,
             mean_ref = mean(x), mean_comp = mean(y),
             fold = mean(y) / mean(x),
             diff_pp = mean(y) - mean(x),
             n_ref = length(x), n_comp = length(y),
             t = unname(tp$statistic), df = unname(tp$parameter),
             P_pooled = tp$p.value, P_welch = tw$p.value,
             stringsAsFactors = FALSE)
}

pairs <- list(c("nTg Vehicle",  "nTg RA013915"),    # the Extended Data comparison
              c("hTau Vehicle", "hTau RA013915"),   # the same treatment in disease
              c("nTg Vehicle",  "hTau Vehicle"))    # the disease effect
tests <- do.call(rbind, c(
  lapply(pairs, function(p) contrast(p[1], p[2], "prom_1kb")),
  lapply(pairs, function(p) contrast(p[1], p[2], "prom_3kb"))))

cat("\n--- contrasts, mouse as the experimental unit ---\n")
for (i in seq_len(nrow(tests))) with(tests[i, ],
  cat(sprintf("%-9s %-14s -> %-14s %5.2f%% vs %5.2f%%  %.2f-fold  (%+.2f pp)  t(%g) = %.2f  P = %s   [Welch P = %s]\n",
              sub("prom_", "<=", measure), reference, comparison,
              mean_comp, mean_ref, fold, diff_pp, df, t,
              format.pval(P_pooled, digits = 3), format.pval(P_welch, digits = 3))))
write.csv(tests, "FigED5d_PromoterFraction_tests.csv", row.names = FALSE)

# A note worth keeping in view when reading the nTg contrast. The treated
# libraries call far fewer peaks, and a less sensitive call retains the
# strongest regions, which are promoter-biased. A shift toward promoters and a
# smaller set are not independent observations, so the legend should carry the
# peak numbers next to the percentages, and the claim is a redistribution of
# the called set rather than a change in how many promoters carry the mark.
cat("\nper-library peak numbers are in the table above and in the CSV; report\n",
    "them next to the percentages.\n", sep = "")

# ---- 4. Prism ----------------------------------------------------------------
# One column per group, one row per mouse, which is the layout Prism expects
# for a grouped scatter with mean and SD.
wide <- function(measure) {
  m <- lapply(levels(res$group), function(g) res[[measure]][res$group == g])
  names(m) <- levels(res$group)
  n <- max(vapply(m, length, integer(1)))
  as.data.frame(lapply(m, function(v) c(v, rep(NA_real_, n - length(v)))),
                check.names = FALSE)
}
p1 <- wide("prom_1kb"); p3 <- wide("prom_3kb")
names(p1) <- paste(names(p1), "<=1kb"); names(p3) <- paste(names(p3), "<=3kb")
write.csv(cbind(p1, p3), "FigED5d_PromoterFraction_PRISM.csv", row.names = FALSE, na = "")

cat("\nwritten to", normalizePath("."), ":\n",
    " FigED5d_PromoterFraction_per_mouse.csv\n",
    " FigED5d_PromoterFraction_tests.csv\n",
    " FigED5d_PromoterFraction_PRISM.csv\n")
