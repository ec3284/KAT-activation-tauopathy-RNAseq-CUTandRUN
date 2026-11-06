#!/usr/bin/env Rscript
# =============================================================================
# TSS metaplot of H3K4ac, with the IgG control (Fig. 4a)
#
# The R replacement for computeMatrix + plotProfile. It reads the spike-in
# normalized bigWig files directly, so nothing here renormalizes anything: the
# scaling is already in the tracks.
#
# Second version. The first used genomation::ScoreMatrixBin, which failed, and
# because the failure happened inside an lapply the object `prof` was never
# created; every later line then referred to `df`, found the stats::df function
# instead, and produced a cascade of "closure is not subsettable". The binning
# is now done by rtracklayer's own bigWig summary, which is the same operation
# computeMatrix performs and needs no extra package. Each track is wrapped so
# that a failure prints its real error and stops, instead of being discovered
# forty lines later.
#
# Strand is handled: profiles over minus-strand transcripts are reversed, so
# upstream is on the left for every gene.
#
# Input : the per-condition averaged bigWigs from bigwigAverage, plus the IgG
#         bigWigs (one library per genotype)
# Output: Fig4a_TSS_metaplot.pdf / .png
#         Fig4a_TSS_metaplot_profiles.csv   the plotted values, for source data
# =============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(GenomicRanges); library(rtracklayer)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(ggplot2)
})
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# ---- PATHS -------------------------------------------------------------------
# Run from the directory where the figure should be written. BIGWIG is the
# folder holding the spike-in normalized tracks: the four per-condition
# averages from 01_, the twelve per-replicate tracks, and the IgG tracks.
# Set CUTRUN_BIGWIG to override it. Nothing is written to that folder.
BIGWIG <- Sys.getenv("CUTRUN_BIGWIG", "bigwig")
if (!dir.exists(BIGWIG)) BIGWIG <- "."

# The tracks to draw, in plotting order, with the palette. IgG is dashed in a
# darker grey rather than #CFCFCF, which is the nTg vehicle colour: two greys
# in one panel cannot be told apart, and the dash says "control" without the
# legend having to.
#
# The averaged file names are the ones that actually exist in the GEO folder:
# nTg_Veh_avg.bw, nTg_915_avg.bw, hTau_Veh_avg.bw, hTau_915_avg.bw.
#
# The IgG track is the same file the pyGenomeTracks configs in this folder use,
# hTau-IgG-rep1.bw. An earlier run of this script did not find it, which was a
# path problem rather than a missing file: locate() below searches both the .bw
# and the .bigWig spelling, but only under BIGWIG and the directories listed in
# ROOTS. If it still is not found, point CUTRUN_BIGWIG at the folder holding it
# before concluding it does not exist, and do not let the legend mention a
# control that the panel does not draw.
TRACKS <- list(
  list(label = "nTg Veh",       file = "nTg_Veh_avg.bw",        color = "#CFCFCF", dash = FALSE),
  list(label = "hTau Veh",      file = "hTau_Veh_avg.bw",       color = "#FFE0C0", dash = FALSE),
  list(label = "hTau RA013915", file = "hTau_915_avg.bw",       color = "#FF6A00", dash = FALSE),
  list(label = "IgG",           file = "hTau-IgG-rep1.bw",      color = "#4D4D4D", dash = TRUE))

# Four-condition version for Extended Data: uncomment instead of the above.
# TRACKS <- list(
#   list(label = "nTg Veh",       file = "nTg_Veh_avg.bw",   color = "#CFCFCF", dash = FALSE),
#   list(label = "nTg RA013915",  file = "nTg_915_avg.bw",   color = "#056AE7", dash = FALSE),
#   list(label = "hTau Veh",      file = "hTau_Veh_avg.bw",  color = "#FFE0C0", dash = FALSE),
#   list(label = "hTau RA013915", file = "hTau_915_avg.bw",  color = "#FF6A00", dash = FALSE))
#
# Per-replicate version, in the "bigwig" subfolder, if you would rather average
# the profiles than the tracks. The mean of the three profiles equals the
# profile of the averaged track, so either route gives the same panel:
#   nTg-Veh-H3K4ac-rep1.bw ... rep3.bw, and the same pattern for the other
#   three groups.

UPSTREAM   <- 3000
DOWNSTREAM <- 3000
N_BIN      <- 300      # 20 bp per bin. computeMatrix used 10 bp; the mean
                       # profile is identical and this runs twice as fast.
SUBSAMPLE  <- 40000    # TSS drawn at random, with the seed below. This is the
                       # value used for the published panel, so the legend must
                       # say a random sample of 40,000 annotated transcripts
                       # and not "all annotated transcripts". Set it to NA to
                       # use every TSS; the mean profile is visually
                       # indistinguishable, but then the legend changes too.
SEED       <- 1
CHUNK      <- 4000     # windows per call, to keep memory flat

OUT <- "Fig4a_TSS_metaplot"

# ---- 1. the TSS set ----------------------------------------------------------
# Transcripts, not genes, which is what computeMatrix does when given a GTF:
# every transcript contributes its own TSS.
tx  <- transcripts(txdb)
tx  <- tx[as.character(seqnames(tx)) %in% paste0("chr", c(1:19, "X"))]
tss <- promoters(tx, upstream = 0, downstream = 1)
cat(sprintf("transcripts: %s\n", format(length(tss), big.mark = ",")))

if (!is.na(SUBSAMPLE) && length(tss) > SUBSAMPLE) {
  set.seed(SEED)
  tss <- tss[sort(sample(length(tss), SUBSAMPLE))]
  cat(sprintf("subsampled to %s TSS (seed %d)\n",
              format(length(tss), big.mark = ","), SEED))
}

win <- resize(tss, width = UPSTREAM + DOWNSTREAM, fix = "center")
sl  <- seqlengths(txdb)[as.character(seqnames(win))]
win <- win[start(win) >= 1 & !is.na(sl) & end(win) <= sl]
cat(sprintf("windows: %s of width %d\n", format(length(win), big.mark = ","),
            UPSTREAM + DOWNSTREAM))

# ---- 2. locate the files -----------------------------------------------------
# The bigWigs are looked for in BIGWIG, then in the working directory, then by
# searching ROOTS recursively, because the averaged tracks are not always in
# the folder that holds the peaks. The search is by file name, so a track found
# in an unexpected place is still the right track, and the directory each file
# came from is printed below for checking.
ROOTS <- unique(c(BIGWIG, ".", dirname(BIGWIG), file.path(".", "bigwig")))
ROOTS <- ROOTS[dir.exists(ROOTS)]

found <- character(0)
locate <- function(nm) {
  for (d in c(BIGWIG, ".")) {
    p <- file.path(d, nm); if (file.exists(p)) return(p)
  }
  if (!length(found)) {
    cat("searching for the bigWig files under:\n")
    for (d in ROOTS) cat("  ", d, "\n")
    # Both extensions, because the H3K4ac tracks are .bw while the IgG tracks
    # were written as .bigWig.
    found <<- unlist(lapply(ROOTS, list.files, pattern = "\\.(bw|bigWig|bigwig)$",
                            recursive = TRUE, full.names = TRUE))
    found <<- unique(found)
    cat(sprintf("  %d .bw files found\n", length(found)))
  }
  hit <- found[basename(found) == nm]
  if (length(hit)) return(hit[1])
  NA_character_
}
files <- vapply(TRACKS, function(t) locate(t$file), character(1))

if (anyNA(files)) {
  cat("\nnot found:\n")
  for (i in which(is.na(files))) cat("  ", TRACKS[[i]]$file, "\n")
  if (length(found)) {
    cat("\n.bw files that do exist, with the directory holding them:\n")
    for (f in utils::head(found, 40))
      cat(sprintf("  %-42s %s\n", basename(f), dirname(f)))
    if (length(found) > 40) cat(sprintf("  ... and %d more\n", length(found) - 40))
  }
  stop("set CUTRUN_BIGWIG to the directory listed above, or correct the names in ",
       "TRACKS. If the averaged tracks were never written, run bigwigAverage ",
       "first, or point TRACKS at the per-replicate bigWigs and the script ",
       "will average them by taking the mean profile.", call. = FALSE)
}

cat("bigWig files:\n")
for (f in files)
  cat(sprintf("  %-38s %7.0f MB  %s\n", basename(f), file.size(f) / 1e6,
              dirname(normalizePath(f, winslash = "/"))))
# A size of 0 MB, or a few KB for a track that should be tens of megabytes, is
# a cloud-storage placeholder rather than the file itself.

# ---- 3. one profile per track ------------------------------------------------
# rtracklayer's summary() on a BigWigFile returns, for each window, `size`
# binned means. That is computeMatrix's reference-point mode. defaultValue = 0
# is deeptools' --missingDataAsZero: a bin with no coverage counts as zero
# rather than dropping the window from the average.
# Opening a bigWig can fail with "UCSC library operation failed / Couldn't
# open" even when file.exists() is TRUE. Two causes, and both are handled here:
# the UCSC code wants an absolute path and chokes on a relative one, and a file
# held by cloud storage as an on-demand placeholder exists for R but has no
# bytes to read until something copies it. The copy below forces it to
# materialize and gives the library a short ASCII path at the same time.
open_bw <- function(f) {
  a <- normalizePath(f, winslash = "/", mustWork = TRUE)
  bwf <- tryCatch({ x <- BigWigFile(a); invisible(seqinfo(x)); x },
                  error = function(e) NULL, warning = function(w) NULL)
  if (!is.null(bwf)) return(bwf)
  tmp <- file.path(normalizePath(tempdir(), winslash = "/"), basename(f))
  cat("\n    opening in place failed, copying to a local temporary file ... ")
  if (!file.copy(a, tmp, overwrite = TRUE))
    stop("could not copy ", basename(f), " to ", tmp, call. = FALSE)
  cat(sprintf("%.0f MB\n    ", file.size(tmp) / 1e6))
  x <- BigWigFile(tmp); invisible(seqinfo(x)); x
}

# summary() on a BigWigFile returns a different class depending on the version
# of rtracklayer: a matrix, a NumericList, or an RleList. unlist() has no
# method for some of them, which is where "no method for coercing this S4 class
# to a vector" comes from. This turns any of the three into a plain matrix.
as_bin_matrix <- function(s, nb) {
  if (is.matrix(s)) return(s)
  v <- tryCatch(as.numeric(unlist(s, use.names = FALSE)),
                error = function(e) NULL)
  if (is.null(v))
    v <- unlist(lapply(as.list(s), function(x) as.numeric(as.vector(x))),
                use.names = FALSE)
  if (length(v) %% nb != 0)
    stop("summary() returned ", length(v), " values, not a multiple of ", nb,
         call. = FALSE)
  matrix(v, ncol = nb, byrow = TRUE)
}

# Newer rtracklayer accepts as = "matrix", which skips the conversion
# altogether. Older versions reject the argument, so it is tried once.
bw_summary <- function(bwf, w, nb) {
  s <- tryCatch(summary(bwf, which = w, size = nb, type = "mean",
                        defaultValue = 0, as = "matrix"),
                error = function(e)
                  summary(bwf, which = w, size = nb, type = "mean",
                          defaultValue = 0))
  as_bin_matrix(s, nb)
}

profile_of <- function(f, label) {
  bwf <- open_bw(f)
  si  <- seqinfo(bwf)
  keep <- as.character(seqnames(win)) %in% seqnames(si)
  if (!any(keep))
    stop("no shared chromosome names between ", basename(f),
         " and the annotation. The bigWig uses: ",
         paste(utils::head(seqnames(si), 5), collapse = ", "),
         call. = FALSE)
  if (!all(keep))
    cat(sprintf("    %s of %s windows are on chromosomes absent from this file\n",
                format(sum(!keep), big.mark = ","),
                format(length(win), big.mark = ",")))
  w <- win[keep]
  neg <- as.character(strand(w)) == "-"

  total <- numeric(N_BIN); n <- 0L
  idx <- split(seq_along(w), ceiling(seq_along(w) / CHUNK))
  for (k in idx) {
    m <- bw_summary(bwf, w[k], N_BIN)
    m[is.na(m)] <- 0
    rev_k <- neg[k]
    if (any(rev_k)) m[rev_k, ] <- m[rev_k, N_BIN:1, drop = FALSE]
    total <- total + colSums(m); n <- n + nrow(m)
  }
  total / n
}

# A probe on three windows before the real work, so that a problem with the
# return class shows up here, with the class named, rather than inside the loop.
{
  pb <- open_bw(files[1])
  s  <- tryCatch(summary(pb, which = win[1:3], size = 10, type = "mean",
                         defaultValue = 0, as = "matrix"),
                 error = function(e)
                   summary(pb, which = win[1:3], size = 10, type = "mean",
                           defaultValue = 0))
  cat(sprintf("probe on %s: summary() returned class %s\n",
              basename(files[1]), paste(class(s), collapse = "/")))
  pm <- as_bin_matrix(s, 10)
  cat(sprintf("  converted to a %d x %d matrix, values %.4g to %.4g\n",
              nrow(pm), ncol(pm), min(pm, na.rm = TRUE), max(pm, na.rm = TRUE)))
  if (all(pm == 0, na.rm = TRUE))
    cat("  every probe value is zero. If this persists for whole tracks, check\n",
        "  that the chromosome naming in the bigWig matches the annotation.\n",
        sep = "")
}

prof <- vector("list", length(TRACKS))
for (i in seq_along(TRACKS)) {
  cat(sprintf("  %-38s ", basename(files[i])))
  prof[[i]] <- tryCatch(profile_of(files[i], TRACKS[[i]]$label),
    error = function(e) stop("failed on ", basename(files[i]), ": ",
                             conditionMessage(e), call. = FALSE))
  pos0 <- seq(-UPSTREAM, DOWNSTREAM, length.out = N_BIN)
  cat(sprintf("peak %.3f at %+d bp\n", max(prof[[i]]),
              round(pos0[which.max(prof[[i]])])))
}

pos <- seq(-UPSTREAM, DOWNSTREAM, length.out = N_BIN)
dat <- do.call(rbind, lapply(seq_along(TRACKS), function(i)
  data.frame(position = pos, signal = prof[[i]],
             track = TRACKS[[i]]$label, stringsAsFactors = FALSE)))
dat$track <- factor(dat$track, levels = vapply(TRACKS, `[[`, character(1), "label"))
write.csv(dat, paste0(OUT, "_profiles.csv"), row.names = FALSE)
# The data frame is called `dat` and not `df` on purpose: `df` is the name of a
# base R function, so a script that fails before creating it carries on with
# the function in its place and reports errors that point nowhere near the
# real one.

# ---- 4. the panel ------------------------------------------------------------
cols <- setNames(vapply(TRACKS, `[[`, character(1), "color"), levels(dat$track))
ltys <- setNames(ifelse(vapply(TRACKS, `[[`, logical(1), "dash"), "22", "solid"),
                 levels(dat$track))

# A line in #CFCFCF on white is faint. OUTLINE draws a dark line underneath
# every curve so the pale ones stay visible in print.
OUTLINE <- FALSE

p <- ggplot(dat, aes(position, signal)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4)
if (OUTLINE)
  p <- p + geom_line(aes(group = track), linewidth = 2.0, colour = "grey25")
p <- p +
  geom_line(aes(colour = track, linetype = track), linewidth = 1.4) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_linetype_manual(values = ltys, name = NULL) +
  scale_x_continuous(breaks = c(-3000, -1500, 0, 1500, 3000),
                     labels = c("-3.0", "-1.5", "TSS", "+1.5", "+3.0"),
                     expand = expansion(0)) +
  labs(x = "Distance from TSS (kb)", y = "Spike-in-normalized H3K4ac") +
  theme_classic(base_size = 20) +
  theme(text = element_text(family = "Arial"),
        axis.text  = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 20),
        axis.line  = element_line(colour = "black", linewidth = 0.6),
        axis.ticks = element_line(colour = "black"),
        legend.text = element_text(size = 18),
        legend.key.width = unit(1.6, "cm"),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_blank())

ggsave(paste0(OUT, ".pdf"), p, width = 9, height = 6.5, device = cairo_pdf)
ggsave(paste0(OUT, ".png"), p, width = 9, height = 6.5, dpi = 400)

cat("\nwritten to", normalizePath("."), ":\n ", paste0(OUT, ".pdf"),
    "\n ", paste0(OUT, ".png"), "\n ", paste0(OUT, "_profiles.csv"), "\n")

# Read the per-track line printed above before trusting the panel. If the IgG
# maximum is not low and its profile not flat relative to the three H3K4ac
# tracks, either the wrong file was picked up or that library is not behaving
# as a control, and the sentence in the Results that depends on it should not
# be written yet.
