# =============================================================================
# Where the regions that respond to RA013915 sit in the genome (Fig. 4c)
#
# The panel asks whether the regions whose H3K4ac signal responds to the
# compound are distributed like the rest of the peaks, and it is built so that
# the answer can come out either way.
#
# The responsive set is defined by signal, never by position: it is the peaks
# with the smallest P values in the DESeq2 comparison of RA013915 against
# vehicle in hTau mice, over the union of consensus peaks, with spike-in size
# factors and the mouse as the experimental unit. An earlier version of this
# panel defined its set by intersecting peaks with promoter windows, which
# made the set 100% promoter by construction; nothing here uses a TSS to
# choose a region, only to classify one afterwards.
#
# Two thresholds matter and both are reported. No region reaches FDR < 0.05 in
# any contrast, so the set is defined at a nominal threshold and named as such
# throughout. At P < 0.01, about 1,012 of the 3,587 regions are expected by
# chance; that admixture dilutes any real enrichment rather than creating one,
# so the value measured here is a lower bound on the enrichment among regions
# that genuinely respond.
#
# The non-responsive set (P > 0.5) is the internal control. If the enrichment
# in the responsive set were an artefact of how peaks are distributed, the
# control would move the same way. It moves the opposite way.
#
# The comparison is made against a coverage-matched background, because a
# region with more reads has more power and so a smaller P, and promoter peaks
# carry more reads than distal ones. Without that matching the panel would
# measure sequencing depth. With it, mean baseMean is 19.5 in both sets.
#
# Input : DiffPeaks_RA013915_in_hTau.csv (from 04_differential_peaks_spikein.R)
# Output: Fig4c_FeatureDist_responsive.pdf / .png / .csv
#         peaks_responsive_p01.bed, peaks_nonresponsive_p50.bed
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(ChIPseeker); library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db); library(GenomicRanges); library(rtracklayer); library(ggplot2)
})
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# ---- PATHS -------------------------------------------------------------------
# Run from the same directory as 05_, which writes the table read here.
RES <- "DiffPeaks_RA013915_in_hTau.csv"
if (!file.exists(RES))
  stop("run 04_differential_peaks_spikein.R first; it writes ", RES, call. = FALSE)

P_RESPONSIVE <- 0.01    # nominal
P_CONTROL    <- 0.5     # nominal, the least responsive regions
N_MATCH      <- 5       # background regions matched per responsive region

# ---- 1. the three sets -------------------------------------------------------
d <- read.csv(RES, stringsAsFactors = FALSE)
d <- d[!is.na(d$pvalue), ]
cat(sprintf("regions with a P value: %s\n", format(nrow(d), big.mark = ",")))
cat(sprintf("expected below P < %.2f by chance: %s\n", P_RESPONSIVE,
            format(round(P_RESPONSIVE * nrow(d)), big.mark = ",")))

resp <- which(d$pvalue < P_RESPONSIVE)
ctrl <- which(d$pvalue > P_CONTROL)
cat(sprintf("responsive (nominal P < %.2f): %s\n", P_RESPONSIVE,
            format(length(resp), big.mark = ",")))
cat(sprintf("non-responsive (P > %.1f)   : %s\n", P_CONTROL,
            format(length(ctrl), big.mark = ",")))

gr <- function(i) GRanges(d$chrom[i], IRanges(d$start[i] + 1, d$end[i]))
export(gr(resp), "peaks_responsive_p01.bed", format = "BED")
export(gr(ctrl), "peaks_nonresponsive_p50.bed", format = "BED")
export(gr(seq_len(nrow(d))), "peaks_all_tested.bed", format = "BED")

# ---- 2. coverage-matched background -----------------------------------------
# For each responsive region, the N_MATCH regions with the closest baseMean
# among those not in the set, each used once.
pool <- setdiff(seq_len(nrow(d)), resp)
pool <- pool[order(d$baseMean[pool])]
pbm  <- d$baseMean[pool]
taken <- logical(length(pool))
for (i in resp) {
  j <- findInterval(d$baseMean[i], pbm); lo <- j; hi <- j + 1L; got <- 0L
  while (got < N_MATCH && (lo >= 1L || hi <= length(pool))) {
    pick <- if (hi > length(pool)) lo else if (lo < 1L) hi else
      if (abs(pbm[lo] - d$baseMean[i]) <= abs(pbm[hi] - d$baseMean[i])) lo else hi
    if (!taken[pick]) { taken[pick] <- TRUE; got <- got + 1L }
    if (pick == lo) lo <- lo - 1L else hi <- hi + 1L
  }
}
match_idx <- pool[taken]
cat(sprintf("coverage-matched background: %s regions | mean baseMean %.1f vs %.1f in the responsive set\n",
            format(length(match_idx), big.mark = ","),
            mean(d$baseMean[match_idx]), mean(d$baseMean[resp])))
export(gr(match_idx), "peaks_coverage_matched.bed", format = "BED")

# ---- 3. annotate -------------------------------------------------------------
sets <- list("All peaks" = "peaks_all_tested.bed",
             "RA013915-responsive" = "peaks_responsive_p01.bed",
             "Coverage-matched" = "peaks_coverage_matched.bed",
             "Non-responsive" = "peaks_nonresponsive_p50.bed")
anno <- lapply(sets, function(f)
  annotatePeak(f, TxDb = txdb, tssRegion = c(-3000, 3000),
               annoDb = "org.Mm.eg.db", verbose = FALSE))
names(anno) <- sprintf("%s\n(n = %s)", names(anno),
  format(vapply(anno, function(a) length(a@anno), integer(1)), big.mark = ","))

# ---- 4. the test -------------------------------------------------------------
# Proximal is within 1 kb of a TSS, the definition used in the published
# panels this one follows. The test is the responsive set against the
# coverage-matched background, which is the comparison the panel is for.
prox <- function(a) {
  s <- a@annoStat
  sum(s$Frequency[s$Feature == "Promoter (<=1kb)"]) / 100
}
n_of <- function(a) length(a@anno)
i_r <- grep("^RA013915-responsive", names(anno))
i_m <- grep("^Coverage-matched", names(anno))
i_c <- grep("^Non-responsive", names(anno))
i_a <- grep("^All peaks", names(anno))

tab <- function(x, y) {
  m <- matrix(round(c(n_of(anno[[x]]) * prox(anno[[x]]), n_of(anno[[x]]) * (1 - prox(anno[[x]])),
                      n_of(anno[[y]]) * prox(anno[[y]]), n_of(anno[[y]]) * (1 - prox(anno[[y]])))),
              nrow = 2, byrow = TRUE,
              dimnames = list(c(sub("\n.*", "", names(anno)[x]), sub("\n.*", "", names(anno)[y])),
                              c("proximal", "distal")))
  ft <- fisher.test(m)
  cat("\n"); print(m)
  cat(sprintf("%.1f%% vs %.1f%%   OR = %.2f   P = %s\n",
              100 * prox(anno[[x]]), 100 * prox(anno[[y]]),
              unname(ft$estimate), format.pval(ft$p.value, digits = 3)))
  invisible(ft)
}
cat("\n--- responsive against coverage-matched background (the panel's test) ---")
tab(i_r, i_m)
cat("\n--- non-responsive against all peaks (internal control) ---")
tab(i_c, i_a)

# ---- 5. panel ----------------------------------------------------------------
p <- plotAnnoBar(anno) +
  theme(text = element_text(size = 20, family = "Arial"),
        axis.text = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 16), legend.title = element_text(size = 18))
ggsave("Fig4c_FeatureDist_responsive.pdf", p, width = 13, height = 5.5, device = cairo_pdf)
ggsave("Fig4c_FeatureDist_responsive.png", p, width = 13, height = 5.5, dpi = 400)

out <- do.call(rbind, lapply(names(anno), function(k) {
  s <- anno[[k]]@annoStat
  data.frame(Set = sub("\n", " ", k), Feature = s$Feature, Percentage = round(s$Frequency, 3))
}))
write.csv(out, "Fig4c_FeatureDist_responsive.csv", row.names = FALSE)
cat("\nwritten to", normalizePath("."), "\n")
