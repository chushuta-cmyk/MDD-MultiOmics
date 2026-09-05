# ==============================================================================
# 03_normalization.R
# DESeq2 normalization: size factor estimation, VST transformation, 
# and export of normalized data for downstream (MOFA / visualization)
# ==============================================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cat("\n=== GSE80655 Normalization Started ===\n")

# ------------------------------------------------------------------------------
# 1. Load preprocessed data from 02
# ------------------------------------------------------------------------------
expr_counts <- readRDS("GSE80655_expr_counts_int.rds")
metadata <- readRDS("GSE80655_metadata_clean.rds")

cat(sprintf("Loaded count matrix: %d genes x %d samples\n", nrow(expr_counts), ncol(expr_counts)))
cat(sprintf("Loaded metadata: %d samples\n", nrow(metadata)))

# ------------------------------------------------------------------------------
# 2. Build DESeq2 object and estimate size factors (internal normalization)
# ------------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = expr_counts,
  colData = metadata,
  design = ~ Diagnosis  # Base design; batch will be added later in 04 if needed
)

# Pre-filter low-count genes (optional but recommended for DESeq2)
keep <- rowSums(counts(dds)) > 1  # Retain genes with at least 1 count across all samples
dds <- dds[keep, ]
cat(sprintf("Pre-filtered low-count genes. Remaining: %d\n", sum(keep)))

# Estimate size factors (normalization factors)
dds <- estimateSizeFactors(dds)
cat("Size factors estimated successfully.\n")
cat("Size factors (first 6 samples):\n")
print(head(sizeFactors(dds)))

# ------------------------------------------------------------------------------
# 3. Apply Variance Stabilizing Transformation (VST) for downstream use
#    VST is recommended for RNA-seq data for PCA, clustering, and MOFA input.
# ------------------------------------------------------------------------------
cat("\nApplying Variance Stabilizing Transformation (VST)...\n")
vsd <- vst(dds, blind = TRUE)  # blind = TRUE for unbiased QC
vsd_matrix <- assay(vsd)

cat(sprintf("VST matrix dimensions: %d genes x %d samples\n", nrow(vsd_matrix), ncol(vsd_matrix)))

# ------------------------------------------------------------------------------
# 4. QC Visualization (Post-Normalization) - PCA
# ------------------------------------------------------------------------------
save_sci_plot <- function(p, f, w = 8, h = 6) {
  ggsave(filename = f, plot = p, device = "tiff", 
         width = w, height = h, dpi = 300, compression = "lzw")
  ggsave(filename = gsub("\\.tiff", ".pdf", f), plot = p, 
         device = "pdf", width = w, height = h)
}

cat("\n--- Generating post-normalization PCA ---\n")
pca_data <- plotPCA(vsd, intgroup = c("Diagnosis", "Brain_Region"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"), 2)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Diagnosis, shape = Brain_Region)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "PCA (VST transformed) - GSE80655",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("Control" = "#2E86AB", "Major_Depression" = "#A23B72"))

save_sci_plot(p_pca, "QC_PCA_VST.tiff")
cat("✅ Saved: QC_PCA_VST.tiff / .pdf\n")

# ------------------------------------------------------------------------------
# 5. Export normalized objects for MOFA / Multi-omics integration
# ------------------------------------------------------------------------------
cat("\n--- Exporting normalized data for MOFA ---\n")

# Save the VST matrix (features x samples)
saveRDS(vsd_matrix, file = "GSE80655_expression_vst.rds")
cat("  - Saved: GSE80655_expression_vst.rds\n")

# Save metadata with proper formatting
saveRDS(metadata, file = "GSE80655_metadata_clean.rds")  # Overwrites with same, but ensures latest
cat("  - Saved: GSE80655_metadata_clean.rds\n")

# Feature annotation (gene symbols)
feature_annot <- data.frame(
  feature_id = rownames(vsd_matrix),
  gene_symbol = rownames(vsd_matrix)
)
saveRDS(feature_annot, file = "GSE80655_feature_annotation.rds")
cat("  - Saved: GSE80655_feature_annotation.rds\n")

# ------------------------------------------------------------------------------
# 6. Checklist
# ------------------------------------------------------------------------------
cat("\n=== CHECKLIST: 03_normalization.R ===\n")
cat("✅ DESeq2 object created with design ~ Diagnosis\n")
cat(sprintf("  - Genes after low-count filter: %d\n", nrow(dds)))
cat(sprintf("  - Samples: %d\n", ncol(dds)))
cat("✅ Size factors estimated (normalization complete)\n")
cat("✅ VST transformation applied\n")
cat(sprintf("  - VST matrix: %d x %d\n", nrow(vsd_matrix), ncol(vsd_matrix)))
cat("✅ Post-normalization PCA plot saved\n")
cat("✅ MOFA-ready objects exported:\n")
cat("  - GSE80655_expression_vst.rds\n")
cat("  - GSE80655_metadata_clean.rds\n")
cat("  - GSE80655_feature_annotation.rds\n")
cat("=== End of 03_normalization.R ===\n")

