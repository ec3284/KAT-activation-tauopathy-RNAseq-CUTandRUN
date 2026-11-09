#!/usr/bin/env Rscript
# =============================================================================
# Raw read counts per mouse over the peak union
#
# This is the step that turns twelve BAM files into the one table the
# differential analysis needs. Each read pair is counted once, over the
# interval it spans, and is credited to every region it overlaps;
# duplicate-flagged alignments are excluded. Counts are raw and stay raw:
# normalization happens in 04_, where the spike-in totals enter DESeq2 as the
# size factors, so anything scaled here would be normalized twice.
#
# Counting runs one chromosome at a time, which keeps memory bounded and never
# splits a read pair across chunks.
#
# Requires: Rsamtools, GenomicAlignments, GenomicRanges, data.table
#
# Input : the twelve coordinate-sorted, indexed H3K4ac BAM files
#         union_relaxed4.bed from 02_
# Output: CUTRUN_union_peak_counts_12mice.csv.gz
#           peak_id, chrom, start, end, then one count column per library,
#           named as the GEO submission names them
#         CUTRUN_spikein_reads.txt   reads mapped to the spike-in genome
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(Rsamtools); library(GenomicAlignments)
  library(GenomicRanges); library(data.table)
})

# ---- PATHS -------------------------------------------------------------------
# Run from the directory where the outputs should be written, or set WORKDIR.
WORKDIR  <- Sys.getenv("WORKDIR", ".")
BED      <- Sys.getenv("CUTRUN_UNION", "union_relaxed4.bed")
BAMDIR   <- Sys.getenv("CUTRUN_BAM", "bam")
SPIKEDIR <- Sys.getenv("CUTRUN_SPIKE_BAM", BAMDIR)
setwd(WORKDIR)
stopifnot(file.exists(BED), dir.exists(BAMDIR))

MIN_COUNT <- 120   # the filter 04_ applies; printed here only as a check

# The twelve libraries, in the order the outputs use. The left-hand name is the
# BAM file stem, the right-hand one is the name the GEO submission uses for that
# library and the name the count column carries. Edit the stems if the pipeline
# wrote different file names; never edit the pairing.
LIBS <- c(
  WT_Veh_R1    = "nTg-Veh-H3K4ac-rep1",
  WT_Veh_R2    = "nTg-Veh-H3K4ac-rep2",
  WT_Veh_R3    = "nTg-Veh-H3K4ac-rep3",
  WT_915_R1    = "nTg-RA013915-H3K4ac-rep1",
  WT_915_R2    = "nTg-RA013915-H3K4ac-rep2",
  WT_915_R3    = "nTg-RA013915-H3K4ac-rep3",
  hTau_Veh_R1  = "hTau-Veh-H3K4ac-rep1",
  hTau_Veh_R2  = "hTau-Veh-H3K4ac-rep2",
  hTau_Veh_R3  = "hTau-Veh-H3K4ac-rep3",
  hTau_915_R1  = "hTau-RA013915-H3K4ac-rep1",
  hTau_915_R2  = "hTau-RA013915-H3K4ac-rep2",
  hTau_915_R3  = "hTau-RA013915-H3K4ac-rep3")

# ---- 1. the regions ----------------------------------------------------------
# BED is 0-based half-open and GRanges is 1-based inclusive, so the start is
# shifted by one for counting; the table written at the end carries the BED's
# own coordinates.
bed <- read.table(BED, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
regions <- GRanges(bed$V1, IRanges(bed$V2 + 1L, bed$V3))
cat(sprintf("regions: %s\n", format(length(regions), big.mark = ",")))

# ---- 2. the BAM files --------------------------------------------------------
# Each library's BAM is located by its stem, allowing for the suffixes
# nf-core/cutandrun appends. The mapping is printed before anything is counted,
# so a wrong stem is visible immediately rather than producing silently
# mislabelled columns.
locate <- function(stem) {
  hits <- list.files(BAMDIR, pattern = paste0("^", stem, "\\..*\\.bam$|^", stem, "\\.bam$"),
                     full.names = TRUE)
  if (!length(hits)) stop("no BAM for ", stem, " in ", normalizePath(BAMDIR), call. = FALSE)
  hits[1]
}
bams <- vapply(names(LIBS), locate, character(1))
print(data.frame(column = unname(LIBS), file = basename(bams)), row.names = FALSE)

if (!all(file.exists(paste0(bams, ".bai")) |
         file.exists(sub("\\.bam$", ".bai", bams))))
  stop("some BAM files are not indexed. Run samtools index first.", call. = FALSE)

# ---- 3. counting -------------------------------------------------------------
count_one <- function(bam, gr) {
  flag <- scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE,
                      isNotPassingQualityControls = FALSE)
  len  <- seqlengths(seqinfo(BamFile(bam)))
  out  <- integer(length(gr))
  chroms <- unique(as.character(seqnames(gr)))
  missing <- setdiff(chroms, names(len))
  if (length(missing))
    warning("not in the header of ", basename(bam), ": ",
            paste(missing, collapse = ", "), call. = FALSE)
  for (ch in intersect(chroms, names(len))) {
    sel <- as.character(seqnames(gr)) == ch
    win <- GRanges(ch, IRanges(1, min(max(end(gr[sel])) + 10000L, len[[ch]])))
    frg <- granges(readGAlignmentPairs(bam, param = ScanBamParam(which = win, flag = flag)),
                   on.discordant.seqnames = "drop")
    if (length(frg)) out[sel] <- countOverlaps(gr[sel], frg, ignore.strand = TRUE)
  }
  out
}

m <- matrix(0L, nrow = length(regions), ncol = length(LIBS),
            dimnames = list(NULL, unname(LIBS)))
for (i in seq_along(bams)) {
  m[, i] <- count_one(bams[i], regions)
  cat(sprintf("  %-28s %s fragments in regions\n", LIBS[i],
              format(sum(m[, i]), big.mark = ",")))
}

counts <- data.table(peak_id = sprintf("pk%d", seq_len(nrow(bed))),
                     chrom   = bed$V1,
                     start   = bed$V2,
                     end     = bed$V3)
counts <- cbind(counts, as.data.table(m))
fwrite(counts, "CUTRUN_union_peak_counts_12mice.csv.gz")

cat(sprintf("\nregions with >= %d reads across the twelve libraries: %s\n",
            MIN_COUNT, format(sum(rowSums(m) >= MIN_COUNT), big.mark = ",")))

# ---- 4. spike-in totals ------------------------------------------------------
# The numbers hard-coded as SPIKE in 04_. They come from the alignment against
# the spike-in genome, which nf-core/cutandrun reports per library; recomputing
# them here lets the two be checked against each other.
spike <- list.files(SPIKEDIR, pattern = "spikein.*\\.bam$", full.names = TRUE)
if (length(spike)) {
  tot <- vapply(names(LIBS), function(st) {
    f <- spike[grepl(paste0("^", st, "\\."), basename(spike))]
    if (!length(f)) return(NA_integer_)
    as.integer(countBam(f[1], param = ScanBamParam(
      flag = scanBamFlag(isUnmappedQuery = FALSE,
                         isSecondaryAlignment = FALSE)))$records)
  }, integer(1))
  fwrite(data.table(library = unname(LIBS), spikein_reads = unname(tot)),
         "CUTRUN_spikein_reads.txt", sep = "\t")
  print(data.table(library = unname(LIBS), spikein_reads = unname(tot)))
} else {
  message("no spike-in BAMs under ", SPIKEDIR,
          "; take the totals from the nf-core/cutandrun alignment report and ",
          "check them against SPIKE in 04_.")
}

cat("\nwritten: CUTRUN_union_peak_counts_12mice.csv.gz\n")
