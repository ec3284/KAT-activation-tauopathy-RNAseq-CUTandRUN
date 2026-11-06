# =============================================================================
# Genomic feature distribution of H3K4ac peaks, four conditions
# Extended Data Fig. 5b
#
# Entirely in R. No bedtools, no Python.
#
# Three changes from the first version.
#
# 1. A row of randomly placed regions. Saying that 16 to 17% of peaks fall in
#    promoters states a proportion, not an enrichment: randomly placed regions
#    also land in promoters at some rate, and without that rate the number
#    cannot support the phrase "promoter-associated". The shuffled set has the
#    same number of regions and the same width distribution as the real peak
#    union, placed at random positions on the same chromosomes, so the
#    difference between the two rows is the enrichment. Quote the ratio in the
#    text rather than the raw percentage.
#
#    Note what the panel will also show: the large majority of peaks are not
#    in promoters. The claim the data support is that the mark is enriched at
#    promoters, not that it is confined to them, and a reader comparing the
#    two rows can see which claim is being made.
#
# 2. The "DEG-up peaks" row is gone. It was built by intersecting the peak
#    union with promoter windows of TSS +/- 2 kb, so every region in it is a
#    promoter region by construction and the bar reproduced the selection
#    rather than measuring anything. Set DIFF_PEAKS to a set defined by
#    differential signal if one becomes available; leave it NULL and no such
#    row is drawn.
#
# 3. Peak counts appear in the labels. The four condition sets hold 133,194 /
#    23,458 / 83,121 / 191,623 regions at comparable sequencing depth and
#    comparable fraction of reads in peaks (about 10% in every library), so
#    they differ in how readily each library passed the SEACR relaxed
#    threshold rather than in library quality. Composition is a proportion and
#    so is more robust to that than a count is, but not immune: the smallest
#    set is about three times more promoter-proximal than the others, because
#    a less sensitive call keeps only the strongest peaks and those sit at
#    promoters. Printing n next to each bar lets a reader see this.
#
# Output: FigED5b_feature_distribution.pdf / .png / .csv
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
})

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# ---- PATHS -------------------------------------------------------------------
# Run this script from the directory where the figure outputs should be written.
# The two locations below are read, never written to. Either edit them or set
# the environment variables CUTRUN_PEAKS and CUTRUN_BIGWIG.
#   CUTRUN_PEAKS   the consensus and per-replicate SEACR peak BED files
#   CUTRUN_BIGWIG  the spike-in normalized bigWig tracks
PEAKS  <- Sys.getenv("CUTRUN_PEAKS",  "peaks_and_consensus_peaks")
PEAKS_DIR <- PEAKS
stopifnot(dir.exists(PEAKS_DIR))

DIFF_PEAKS <- NULL   # a set defined by differential signal, when one exists

# The peak files have travelled under two naming schemes, the GEO names and
# the pipeline's own "<group>_relaxed.bed". Rather than fixing one, each group
# is located by pattern, and the script reports what it matched so a silent
# mismatch cannot pass unnoticed.
find_one <- function(patterns, label) {
  f <- unlist(lapply(patterns, function(p)
    list.files(PEAKS_DIR, pattern = p, ignore.case = TRUE, full.names = TRUE)))
  f <- unique(f[file.exists(f)])
  if (length(f) == 0)
    stop(sprintf("no file found for %s in\n  %s\nExpected something matching: %s",
                 label, PEAKS_DIR, paste(patterns, collapse = " or ")), call. = FALSE)
  if (length(f) > 1)
    message(sprintf("several candidates for %s (%s); using the first",
                    label, paste(basename(f), collapse = ", ")))
  f[1]
}

groups <- list(
  "nTg Veh"       = find_one(c("^nTg[-_]Veh.*\\.bed$", "^nTg_Veh_relaxed\\.bed$",
                               "nTg-Veh-H3K4ac-consensus\\.peaks\\.bed$"), "nTg Veh"),
  "nTg RA013915"  = find_one(c("^nTg[-_](RA013915|915).*\\.bed$",
                               "nTg-RA013915-H3K4ac-consensus\\.peaks\\.bed$"), "nTg RA013915"),
  "hTau Veh"      = find_one(c("^hTau[-_]Veh.*\\.bed$",
                               "hTau-Veh-H3K4ac-consensus\\.peaks\\.bed$"), "hTau Veh"),
  "hTau RA013915" = find_one(c("^hTau[-_](RA013915|915).*\\.bed$",
                               "hTau-RA013915-H3K4ac-consensus\\.peaks\\.bed$"), "hTau RA013915"))

cat("peak files used (read from", PEAKS_DIR, "):\n")
for (k in names(groups))
  cat(sprintf("  %-14s %-48s %s regions\n", k, basename(groups[[k]]),
              format(length(import(groups[[k]], format = "BED")), big.mark = ",")))

# Which condition rows the panel shows: all four.
#
# An earlier version drew three and left nTg RA013915 for the Extended Data
# figure. That was a mistake, because the fourth set is the only one that
# differs. Measured as the distance from each peak to the nearest TSS in the
# same annotation, the proportion of peaks within 1 kb of a TSS is
#
#   nTg Veh        8.79%   (133,194 regions)
#   nTg RA013915  15.17%   ( 23,458 regions)
#   hTau Veh       8.21%   ( 83,121 regions)
#   hTau RA013915  8.89%   (191,623 regions)
#
# so three conditions sit within 0.7 percentage points of each other and the
# fourth is 1.73-fold higher. A three-row panel showed only the part of the
# result where nothing happens.
#
# The n of each row has to be printed next to it, and here is why it is not a
# caveat but part of the finding. nTg RA013915 has 23,458 called regions
# against 83,000 to 192,000 for the others. Fewer regions pass the caller's
# threshold in that library, and the regions that do pass are the
# highest-signal ones, which are disproportionately promoters. The shift
# toward promoters and the smaller set are therefore two views of the same
# thing, and the panel is read as a change in the composition of the called
# set, not as a change in the number of promoters that carry H3K4ac.
SHOW_CONDITIONS <- c("nTg Veh", "nTg RA013915", "hTau Veh", "hTau RA013915")
peak_files <- groups[intersect(SHOW_CONDITIONS, names(groups))]

# The union of the four call sets, built here so the script does not depend on
# a file produced elsewhere. mergeBed and GenomicRanges::reduce do the same
# thing: overlapping intervals are collapsed into one.
UNION <- "union_relaxed4.bed"
if (!file.exists(UNION)) {
  u <- reduce(sort(do.call(c, unname(lapply(groups, import, format = "BED")))))
  export(u, UNION, format = "BED")
  cat(sprintf("built %s: %s regions\n", UNION, format(length(u), big.mark = ",")))
}
union_gr0 <- import(UNION, format = "BED")
cat(sprintf("union: %s regions\n", format(length(union_gr0), big.mark = ",")))

# ---- gene sets ---------------------------------------------------------------
# Peaks are selected by the whole gene locus, body plus FLANK on each side,
# never by a promoter window. A set chosen inside promoter windows is 100%
# promoter whatever the data say, which is what the earlier "DEG-up peaks" row
# measured. Selecting by locus leaves the position of the peak free, so the
# panel can come out either way.
#
# The background row is peaks at every other gene's locus, not random genomic
# regions, so that gene loci are compared with gene loci and the comparison is
# not carried by intergenic space.
#
# One caveat to keep in mind when reading the result: promoter share falls
# with gene length, because a long gene has one promoter and many intragenic
# peaks. The leading-edge genes are synaptic and long (median about 48 kb);
# the RNA-seq up genes are shorter. The script prints the median length of
# each set so that a difference in promoter share can be checked against it
# before being read as biology.
FLANK <- 5000
# The gene-set rows are off. They are not on the same footing as the condition
# rows: a condition row holds every peak called in that library, 38% of them
# distal intergenic, while a gene-set row holds only peaks inside gene loci,
# where that share falls to about 7%. Drawn together, the two blocks invite a
# comparison that is carried by that restriction rather than by the genes.
# Set SHOW_GENESETS to TRUE to draw them, and say in the legend which rows are
# restricted to gene loci.
SHOW_GENESETS <- FALSE

GENESETS <- list(
  # The seven genes the new model calls significant for the compound in hTau
  # mice. They carry 34 peaks between them, so each peak moves the bar by 2.9
  # percentage points and a category at 20% has a 95% interval from 7% to 33%.
  # Their median length is also 3,469 bp against about 17,000 for genes in
  # general, and promoter share rises as genes get shorter, so a high promoter
  # bar here would be the length of the genes rather than a property of the
  # compound. The row is drawn because it is worth seeing; it is not a
  # comparison the panel can carry, and the n in the label says why.
  "Rescued genes (7)"  = "rescued_7_genes.txt",
  "Leading-edge genes" = "LE_up_treatment.txt")
# Both files come from the RNA-seq analysis and are only read when
# SHOW_GENESETS is TRUE, which it is not for the published panel. A missing
# file skips its row rather than stopping the script.

g <- genes(txdb, single.strand.genes.only = TRUE)
sym <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db, names(g),
                                              "SYMBOL", "ENTREZID"))
g$symbol <- unname(sym)
g <- g[!is.na(g$symbol)]

peaks_at <- function(symbols) {
  gg <- g[g$symbol %in% symbols]
  if (length(gg) == 0) return(list(gr = GRanges(), n = 0, len = NA))
  loci <- trim(resize(gg, width(gg) + 2 * FLANK, fix = "center"))
  list(gr = union_gr0[queryHits(findOverlaps(union_gr0, loci))] |> unique(),
       n = length(gg), len = median(width(gg)))
}

used <- character(0)
for (k in if (SHOW_GENESETS) names(GENESETS) else character(0)) {
  f <- GENESETS[[k]]
  if (!file.exists(f)) { message("missing gene list, row skipped: ", f); next }
  s <- unique(trimws(readLines(f))); s <- s[nzchar(s)]
  r <- peaks_at(s); used <- union(used, s)
  if (length(r$gr) == 0) { message("no peaks for ", k); next }
  fn <- paste0("peaks_", gsub("[^A-Za-z0-9]+", "_", k), ".bed")
  export(r$gr, fn, format = "BED")
  peak_files[[k]] <- fn
  cat(sprintf("%-20s %s genes matched, median gene length %s bp, %s peaks in locus +/- %s kb\n",
              k, format(r$n, big.mark = ","), format(r$len, big.mark = ","),
              format(length(r$gr), big.mark = ","), FLANK / 1000))
}

# The "all other genes" row is not drawn. It is still computed, because the
# promoter share of the gene-set rows can only be read against it: peaks at
# genes in general sit at about 17% promoter, and a set above or below that is
# what the rows are saying. The number is printed rather than plotted. Set
# SHOW_OTHER_GENES to TRUE to put the bar back.
SHOW_OTHER_GENES <- FALSE
other <- if (SHOW_GENESETS) setdiff(g$symbol, used) else character(0)
r <- peaks_at(other)
if (length(r$gr)) export(r$gr, "peaks_other_genes.bed", format = "BED")
if (SHOW_OTHER_GENES) peak_files[["All other genes"]] <- "peaks_other_genes.bed"
if (length(r$gr))
  cat(sprintf("%-20s %s genes, median gene length %s bp, %s peaks%s\n",
              "All other genes", format(r$n, big.mark = ","),
              format(r$len, big.mark = ","), format(length(r$gr), big.mark = ","),
              if (SHOW_OTHER_GENES) "" else "  [computed, not plotted]"))

# The union itself is not drawn as a row. Every peak in it belongs to one of
# the rows above, so the bar would be the weighted average of them and would
# state nothing the others do not. It stays as the source the rows are drawn
# from, and as the set the random control is matched to.

# ---- 1. randomly placed regions ----------------------------------------------
# Same count and same widths as the union, placed uniformly at random on the
# standard chromosomes, with the chromosome drawn in proportion to its length
# so that random placement is genuinely uniform across the genome. Regions are
# kept inside the chromosome. The seed makes the row reproducible; the result
# is stable because the set is large.
set.seed(1)
union_gr <- union_gr0
sl <- seqlengths(txdb)
sl <- sl[names(sl) %in% paste0("chr", c(1:19, "X", "Y"))]
sl <- sl[!is.na(sl)]

w <- width(union_gr)
chr <- sample(names(sl), length(w), replace = TRUE, prob = sl / sum(sl))
maxstart <- pmax(as.numeric(sl[chr]) - w, 1)
start_r <- floor(runif(length(w)) * maxstart) + 1
shuffled <- GRanges(seqnames = chr,
                    ranges = IRanges(start = start_r, width = w))
shuffled <- sort(shuffled)
export(shuffled, "union_shuffled.bed", format = "BED")
cat(sprintf("shuffled control: %s regions, median width %s bp\n",
            format(length(shuffled), big.mark = ","),
            format(median(width(shuffled)), big.mark = ",")))

# The random row is left out: once each gene set is compared with peaks at all
# other genes, the genomic expectation adds a reference the panel no longer
# needs. Set SHOW_RANDOM to TRUE to draw it.
SHOW_RANDOM <- TRUE
if (SHOW_RANDOM) peak_files[["Random regions"]] <- "union_shuffled.bed"
if (!is.null(DIFF_PEAKS) && file.exists(DIFF_PEAKS))
  peak_files[["RA013915-responsive peaks"]] <- DIFF_PEAKS

# ---- 2. annotate -------------------------------------------------------------
peakAnnoList <- lapply(peak_files, function(f)
  annotatePeak(f, TxDb = txdb, tssRegion = c(-3000, 3000),
               annoDb = "org.Mm.eg.db", verbose = FALSE))

names(peakAnnoList) <- sprintf("%s\n(n = %s)", names(peakAnnoList),
  format(vapply(peakAnnoList, function(a) length(a@anno), integer(1)),
         big.mark = ","))

# ---- 3. the enrichment the panel is for --------------------------------------
# One number for the text: how much more often a real peak falls in a promoter
# than a randomly placed region of the same size does.
prom <- function(a) {
  s <- a@annoStat
  sum(s$Frequency[grepl("^Promoter", s$Feature)])
}
pr <- vapply(peakAnnoList, prom, numeric(1))
names(pr) <- sub("\n.*", "", names(pr))

ref  <- if ("All other genes" %in% names(pr)) pr[["All other genes"]] else
        if (file.exists("peaks_other_genes.bed"))
          prom(annotatePeak("peaks_other_genes.bed", TxDb = txdb,
                            tssRegion = c(-3000, 3000), annoDb = "org.Mm.eg.db",
                            verbose = FALSE)) else NA_real_
rand <- if ("Random regions" %in% names(pr)) pr[["Random regions"]] else NA_real_
cat("\npromoter share (TSS +/- 3 kb)\n")
for (k in names(pr))
  cat(sprintf("  %-22s %5.1f%%\n", k, pr[[k]]))
if (!is.na(ref)) {
  cat("\neach set against peaks at all other genes:\n")
  for (k in setdiff(names(pr), c("All other genes", "Random regions")))
    cat(sprintf("  %-22s %5.1f%% vs %.1f%%  ->  %.2f-fold\n", k, pr[[k]], ref, pr[[k]] / ref))
}
# ---- the enrichment the random row is there for --------------------------
# Peaks against randomly placed regions of the same number and width. This is
# the comparison that supports calling the mark promoter-associated: 16% in
# promoters means nothing on its own, because randomly placed regions also
# land there at some rate, and that rate is what the test is against.
#
# The conditions are deliberately not tested against each other. With more
# than 100,000 regions per set, a difference of one percentage point returns
# P < 1e-15 while the odds ratio is 1.10, and the message of this panel is
# that the conditions do not differ. A small P there would state the opposite
# of what the panel shows.
if ("Random regions" %in% names(pr)) {
  # The union is annotated here for the test only; it is not drawn, because
  # every peak in it already appears in one of the condition rows.
  au <- annotatePeak(UNION, TxDb = txdb, tssRegion = c(-3000, 3000),
                     annoDb = "org.Mm.eg.db", verbose = FALSE)
  nm_r <- grep("^Random regions", names(peakAnnoList), value = TRUE)
  n_u <- length(au@anno); n_r <- length(peakAnnoList[[nm_r]]@anno)
  f_u <- prom(au) / 100; f_r <- pr[["Random regions"]] / 100
  tab <- matrix(round(c(n_u * f_u, n_u * (1 - f_u),
                        n_r * f_r, n_r * (1 - f_r))), nrow = 2, byrow = TRUE,
                dimnames = list(c("peaks", "random"), c("promoter", "other")))
  print(tab)
  ct <- chisq.test(tab)
  orr <- (tab[1,1] / tab[1,2]) / (tab[2,1] / tab[2,2])
  cat(sprintf("promoter enrichment: %.1f%% vs %.1f%%, %.2f-fold, chi-squared = %.0f, P = %s\n",
              100 * f_u, 100 * f_r, orr, unname(ct$statistic),
              format.pval(ct$p.value, digits = 3)))
}

if (!is.na(rand)) {
  cat(sprintf("\nfor scale, randomly placed regions sit at %.1f%%, so peaks at genes in\n", rand))
  cat(sprintf("general are %.2f-fold promoter-enriched over the genomic expectation.\n", ref / rand))
}
cat("Read any difference between gene sets against the median gene lengths printed\n",
    "above: promoter share falls as genes get longer, so two sets of different\n",
    "length are not directly comparable on this axis.\n", sep = "")

# ---- 4. panel ----------------------------------------------------------------
p <- plotAnnoBar(peakAnnoList) +
  theme(text = element_text(size = 20, family = "Arial"),
        axis.text = element_text(size = 20, colour = "black"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18))

ggsave("FigED5b_feature_distribution.pdf", p, width = 13, height = 6.5, device = cairo_pdf)
ggsave("FigED5b_feature_distribution.png", p, width = 13, height = 6.5, dpi = 400)

# ---- 5. source data ----------------------------------------------------------
tab <- do.call(rbind, lapply(names(peakAnnoList), function(k) {
  s <- peakAnnoList[[k]]@annoStat
  data.frame(Set = k, Feature = s$Feature, Percentage = round(s$Frequency, 3))
}))
write.csv(tab, "FigED5b_feature_distribution.csv", row.names = FALSE)

for (k in names(peakAnnoList)) {
  cat("\n===", k, "===\n"); print(peakAnnoList[[k]]@annoStat)
}
cat("\nwritten to", normalizePath("."), "\n")
