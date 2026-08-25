#!/usr/bin/env Rscript
# ==============================================================================
# RNA-seq analysis: final pipeline with technical-replicate collapsing
# and cross-filtered "rescue" gene identification
#
# Design: 20 sequenced samples (nTg Veh/915 from batch C3SVT8; hTau Veh from
# BOTH batch C3SVT8 and VJ385X, same biological RNA re-sequenced; hTau 915
# from batch VJ385X only). The 8 hTau_Veh samples (4 from each batch) are
# collapsed into 4 samples using DESeq2::collapseReplicates(), summing raw
# counts for the same biological sample sequenced twice. This avoids both
# (a) discarding half the hTau_Veh data, and (b) the need to statistically
# estimate a large batch effect (batch_VJ_vs_C3 was found to affect ~3,900
# genes when modeled explicitly, an assumption too strong to rely on with
# only 4 anchor samples).
#
# "Rescue" genes are defined by a cross-filter: significant (padj<0.05) in
# opposite directions between the disease comparison (hTau_Veh vs nTg_Veh)
# and the treatment comparison (hTau_915 vs hTau_Veh) — i.e., genes altered
# by disease and reversed by RA013915. This definition inherently controls
# for batch/technical artifacts, since spurious batch-driven genes are not
# expected to show this specific opposite-direction, dual-significance
# pattern by chance.
# ==============================================================================

library(readr)
library(DESeq2)
library(dplyr)
library(tibble)
library(openxlsx)

# ---- 0. Setup ----
WORKDIR <- Sys.getenv("WORKDIR", unset = "./rnaseq_work")
dir.create(WORKDIR, showWarnings = FALSE, recursive = TRUE)
setwd(WORKDIR)
dir.create(tempdir(), recursive = TRUE, showWarnings = FALSE)  # guards against broken R temp dir

# ---- 1. Load raw count matrices from both Plasmidsaurus batches ----
# C3SVT8: samples 1-4 = nTg_Veh, 5-8 = nTg_915, 9-12 = hTau_Veh
# VJ385X: samples 1-4 = hTau_Veh (same biological RNA as C3SVT8 9-12,
#         re-sequenced), 5-8 = hTau_915
data_c3 <- read_tsv("C3SVT8-expression-matrix.tsv")
data_c3_df <- as.data.frame(data_c3)
gene_names_c3 <- ifelse(is.na(data_c3_df$gene_name), data_c3_df$gene_id, data_c3_df$gene_name)

data_vj <- read_tsv("VJ385X-expression-matrix.tsv")
data_vj_df <- as.data.frame(data_vj)
gene_names_vj <- ifelse(is.na(data_vj_df$gene_name), data_vj_df$gene_id, data_vj_df$gene_name)

# ---- 2. Build per-group count matrices ----
ntg_veh <- data_c3_df[, paste0("C3SVT8_", 1:4, "_count")]
rownames(ntg_veh) <- make.unique(gene_names_c3)
colnames(ntg_veh) <- c("nTg_Veh_1", "nTg_Veh_2", "nTg_Veh_3", "nTg_Veh_4")

ntg_915 <- data_c3_df[, paste0("C3SVT8_", 5:8, "_count")]
rownames(ntg_915) <- make.unique(gene_names_c3)
colnames(ntg_915) <- c("nTg_915_1", "nTg_915_2", "nTg_915_3", "nTg_915_4")

htau_veh_c3 <- data_c3_df[, paste0("C3SVT8_", 9:12, "_count")]
rownames(htau_veh_c3) <- make.unique(gene_names_c3)
colnames(htau_veh_c3) <- c("hTauVeh_C3_1", "hTauVeh_C3_2", "hTauVeh_C3_3", "hTauVeh_C3_4")

htau_veh_vj <- data_vj_df[, paste0("VJ385X_", 1:4, "_count")]
rownames(htau_veh_vj) <- make.unique(gene_names_vj)
colnames(htau_veh_vj) <- c("hTauVeh_VJ_1", "hTauVeh_VJ_2", "hTauVeh_VJ_3", "hTauVeh_VJ_4")

htau_915 <- data_vj_df[, paste0("VJ385X_", 5:8, "_count")]
rownames(htau_915) <- make.unique(gene_names_vj)
colnames(htau_915) <- c("hTau915_VJ_1", "hTau915_VJ_2", "hTau915_VJ_3", "hTau915_VJ_4")

# ---- 3. Merge on common genes, build the full 20-sample matrix ----
common_genes <- Reduce(intersect, list(rownames(ntg_veh), rownames(ntg_915),
                                         rownames(htau_veh_c3), rownames(htau_veh_vj),
                                         rownames(htau_915)))

counts_matrix <- cbind(
  ntg_veh[common_genes, ], ntg_915[common_genes, ],
  htau_veh_c3[common_genes, ], htau_veh_vj[common_genes, ],
  htau_915[common_genes, ]
)
counts_matrix <- round(counts_matrix)

# ---- 4. Sample metadata: assign shared sample_id to the two hTau_Veh
#         sequencing runs of the same biological sample, for collapsing ----
sample_info <- data.frame(
  row.names = colnames(counts_matrix),
  sample_id = factor(c(
    "nTg_Veh_1", "nTg_Veh_2", "nTg_Veh_3", "nTg_Veh_4",
    "nTg_915_1", "nTg_915_2", "nTg_915_3", "nTg_915_4",
    "hTau_Veh_1", "hTau_Veh_2", "hTau_Veh_3", "hTau_Veh_4",   # C3 run
    "hTau_Veh_1", "hTau_Veh_2", "hTau_Veh_3", "hTau_Veh_4",   # VJ run, same RNA
    "hTau_915_1", "hTau_915_2", "hTau_915_3", "hTau_915_4"
  )),
  condition = factor(c(
    rep("nTg_Veh", 4), rep("nTg_915", 4),
    rep("hTau_Veh", 4), rep("hTau_Veh", 4),
    rep("hTau_915", 4)
  ), levels = c("nTg_Veh", "nTg_915", "hTau_Veh", "hTau_915")),
  batch = factor(c(
    rep("C3SVT8", 8), rep("C3SVT8", 4), rep("VJ385X", 4), rep("VJ385X", 4)
  ))
)

# ---- 5. DESeq2: build, collapse technical replicates, filter, run ----
dds_raw <- DESeqDataSetFromMatrix(countData = as.matrix(counts_matrix),
                                    colData = sample_info, design = ~condition)

dds_collapsed <- collapseReplicates(dds_raw, groupby = dds_raw$sample_id)
stopifnot(ncol(dds_collapsed) == 16)  # sanity check: 4 groups x 4 replicates

dds_filtered <- dds_collapsed[rowSums(counts(dds_collapsed)) >= 10, ]
dds <- DESeq(dds_filtered)

# ---- 6. QC: PCA on VST-transformed data ----
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "condition") +
  theme_bw() + ggtitle("PCA post-collapse (n=4 per group)")
ggsave("PCA_post_collapse.pdf", pca_plot, width = 7, height = 6, dpi = 300)

# ---- 7. Extract the two key contrasts ----
res_disease   <- results(dds, contrast = c("condition", "hTau_Veh", "nTg_Veh"))
res_treatment <- results(dds, contrast = c("condition", "hTau_915", "hTau_Veh"))

df_disease <- as.data.frame(res_disease) %>%
  rename_with(~ paste0(., "_disease")) %>% rownames_to_column("gene")
df_treatment <- as.data.frame(res_treatment) %>%
  rename_with(~ paste0(., "_treatment")) %>% rownames_to_column("gene")

df_merged <- inner_join(df_disease, df_treatment, by = "gene")

# ---- 8. Define rescue genes: significant (padj<0.05), opposite direction
#         between disease and treatment ----
rescued_up   <- df_merged %>%  # UP in disease, DOWN with treatment
  filter(padj_disease < 0.05 & log2FoldChange_disease > 0,
         padj_treatment < 0.05 & log2FoldChange_treatment < 0) %>%
  arrange(padj_treatment)

rescued_down <- df_merged %>%  # DOWN in disease, UP with treatment (the main rescue direction)
  filter(padj_disease < 0.05 & log2FoldChange_disease < 0,
         padj_treatment < 0.05 & log2FoldChange_treatment > 0) %>%
  arrange(padj_treatment)

cat("Rescued (UP in disease, reversed DOWN by treatment):", nrow(rescued_up), "\n")
cat("Rescued (DOWN in disease, reversed UP by treatment):", nrow(rescued_down), "\n")
cat("Total rescue genes:", nrow(rescued_up) + nrow(rescued_down), "\n")

# ---- 9. Also report the two comparisons independently (all significant genes) ----
disease_all <- df_merged %>% filter(padj_disease < 0.05) %>% arrange(padj_disease)
treatment_all <- df_merged %>% filter(padj_treatment < 0.05) %>% arrange(padj_treatment)
cat("All disease-significant genes:", nrow(disease_all), "\n")
cat("All treatment-significant genes (no cross-filter):", nrow(treatment_all), "\n")

# ---- 9b. Additional comparison: nTg Veh vs nTg 915 ----
# Expected/confirmed result: 0 significant genes, consistent with the
# chromatin-gated model (KAT activation has minimal transcriptional effect
# in mice with intact chromatin acetylation homeostasis).
res_nTg <- results(dds, contrast = c("condition", "nTg_915", "nTg_Veh"))
df_nTg <- as.data.frame(res_nTg) %>%
  rownames_to_column("gene") %>%
  arrange(padj)
cat("nTg Veh vs nTg 915 significant genes:", sum(df_nTg$padj < 0.05, na.rm = TRUE), "\n")

# ---- 10. Export everything to Excel ----
wb <- createWorkbook()

addWorksheet(wb, "Rescued_UP_in_disease")
writeData(wb, "Rescued_UP_in_disease", rescued_up)

addWorksheet(wb, "Rescued_DOWN_in_disease")
writeData(wb, "Rescued_DOWN_in_disease", rescued_down)

addWorksheet(wb, "All_Disease_Sig")
writeData(wb, "All_Disease_Sig", disease_all)

addWorksheet(wb, "All_Treatment_Sig_NoFilter")
writeData(wb, "All_Treatment_Sig_NoFilter", treatment_all)

addWorksheet(wb, "Full_Merged_Results")
writeData(wb, "Full_Merged_Results", df_merged)

addWorksheet(wb, "nTg_Veh_vs_915_AllGenes")
writeData(wb, "nTg_Veh_vs_915_AllGenes", df_nTg)

addWorksheet(wb, "Sample_Design")
writeData(wb, "Sample_Design", data.frame(Sample = rownames(sample_info), sample_info))

saveWorkbook(wb, "RNAseq_FINAL_DESeq2_results.xlsx", overwrite = TRUE)

# ---- 11. Also save the raw collapsed count matrix (for GEO / source data) ----
counts_export <- data.frame(gene_symbol = rownames(counts(dds_collapsed)),
                              counts(dds_collapsed), check.names = FALSE)
write.csv(counts_export, "RNAseq_FINAL_collapsed_counts_16samples.csv", row.names = FALSE)

cat("\nDone! Files saved in", WORKDIR, ":\n",
    "- RNAseq_FINAL_DESeq2_results.xlsx\n",
    "- RNAseq_FINAL_collapsed_counts_16samples.csv\n",
    "- PCA_post_collapse.pdf\n",
    "\nRun 02_MAplots.R next to generate the MA plot figures.\n")
