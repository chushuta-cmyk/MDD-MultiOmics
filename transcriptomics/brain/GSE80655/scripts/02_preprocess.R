# ==============================================================================
# 02_preprocess.R
# Get phenotype from GEO, filter MDD/Control, align expression matrix,
# QC visualization (heatmap + PCA), batch effect detection, effect size (Cohen's d)
# ==============================================================================

library(GEOquery)
library(stringr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

cat("\n=== GSE80655 Preprocess Started ===\n")

# ------------------------------------------------------------------------------
# 1. Download phenotype data from GEO
# ------------------------------------------------------------------------------
gse_data <- getGEO("GSE80655", GSEMatrix = TRUE)
gse_eset <- gse_data[[1]]
phenotype_df <- pData(gse_eset)

cat(sprintf("Raw phenotype rows: %d\n", nrow(phenotype_df)))

# ------------------------------------------------------------------------------
# 2. Select relevant columns and filter MDD / Control
# ------------------------------------------------------------------------------
phenotype_selected <- phenotype_df[, c("geo_accession", "brain region:ch1", "clinical diagnosis:ch1")]

target_diagnoses <- c("Control", "Major Depression")

metadata_selected <- phenotype_selected[
  phenotype_selected$`clinical diagnosis:ch1` %in% target_diagnoses, 
]

cat(sprintf("Filtered phenotype: %d samples (%d MDD, %d Control)\n",
            nrow(metadata_selected),
            sum(metadata_selected$`clinical diagnosis:ch1` == "Major Depression"),
            sum(metadata_selected$`clinical diagnosis:ch1` == "Control")))

# ------------------------------------------------------------------------------
# 3. Clean column names and create final grouping factor
# ------------------------------------------------------------------------------
colnames(metadata_selected) <- c("GSM_ID", "Brain_Region_Raw", "Diagnosis_Raw")

metadata_selected$Brain_Region <- make.names(gsub(" ", "_", metadata_selected$Brain_Region_Raw))
metadata_selected$Diagnosis <- make.names(gsub(" ", "_", metadata_selected$Diagnosis_Raw))
metadata_selected$Final_Group <- factor(paste(metadata_selected$Brain_Region, 
                                              metadata_selected$Diagnosis, sep = "_"))

sample_metadata_final <- metadata_selected[, c("GSM_ID", "Final_Group", "Brain_Region", "Diagnosis")]

cat("\nFinal group counts:\n")
print(table(sample_metadata_final$Final_Group))

# ------------------------------------------------------------------------------
# 4. Batch effect detection (check if 'batch:ch1' exists in phenotype)
# ------------------------------------------------------------------------------
if ("batch:ch1" %in% colnames(phenotype_df)) {
  cat("\n⚠️ WARNING: 'batch:ch1' column detected in phenotype.\n")
  cat("   Consider including batch in DESeq2 design: '~ Batch + Diagnosis'\n")
  cat("   Please add the following code to 03/04 scripts:\n")
  cat("   sample_metadata_final$Batch <- phenotype_df[rownames(sample_metadata_final), 'batch:ch1']\n")
  cat("   dds <- DESeqDataSetFromMatrix(..., design = ~ Batch + Diagnosis)\n")
} else {
  cat("\n✅ No 'batch:ch1' column found. Assuming no known batch effects.\n")
}

# ------------------------------------------------------------------------------
# 5. Load deduplicated count matrix and align samples
# ------------------------------------------------------------------------------
final_matrix_unique <- readRDS("GSE80655_raw_counts_dedup.rds")

selected_gsm_ids <- sample_metadata_final$GSM_ID
matrix_cols <- colnames(final_matrix_unique)
matching_gsm_ids <- selected_gsm_ids[selected_gsm_ids %in% matrix_cols]

missing_gsm_ids <- selected_gsm_ids[!selected_gsm_ids %in% matrix_cols]
if (length(missing_gsm_ids) > 0) {
  cat(sprintf("⚠️ Warning: %d samples missing in count matrix. Ignored.\n", length(missing_gsm_ids)))
}

expr_data_final <- as.matrix(final_matrix_unique[, matching_gsm_ids])

sample_metadata_final_subset <- sample_metadata_final[
  match(colnames(expr_data_final), sample_metadata_final$GSM_ID), 
]

regions <- unique(sample_metadata_final_subset$Brain_Region)
regions_gsm_list_final <- list()
for (region in regions) {
  regions_gsm_list_final[[region]] <- sample_metadata_final_subset$GSM_ID[
    sample_metadata_final_subset$Brain_Region == region
  ]
}

cat(sprintf("Aligned expression matrix: %d genes x %d samples\n", 
            nrow(expr_data_final), ncol(expr_data_final)))

# ------------------------------------------------------------------------------
# 6. Convert to integer matrix (required for DESeq2)
# ------------------------------------------------------------------------------
expr_data_counts <- round(expr_data_final)
storage.mode(expr_data_counts) <- "integer"
cat("Converted to integer matrix for DESeq2.\n")

# ------------------------------------------------------------------------------
# 7. QC Visualization: Sample correlation heatmap + PCA
# ------------------------------------------------------------------------------
# Save SCI-style plots function (same as GSE54563)
save_sci_plot <- function(p, f, w = 8, h = 6) {
  ggsave(filename = f, plot = p, device = "tiff", 
         width = w, height = h, dpi = 300, compression = "lzw")
  ggsave(filename = gsub("\\.tiff", ".pdf", f), plot = p, 
         device = "pdf", width = w, height = h)
}

cat("\n--- Generating QC plots ---\n")

# 7a. Log2 transform for visualization
log_expr <- log2(expr_data_counts + 1)

# 7b. Sample correlation heatmap
cor_matrix <- cor(log_expr, use = "complete.obs")

annotation_col <- data.frame(
  Diagnosis = sample_metadata_final_subset$Diagnosis,
  Brain_Region = sample_metadata_final_subset$Brain_Region,
  row.names = colnames(log_expr)
)

p_heatmap <- pheatmap(
  cor_matrix,
  annotation_col = annotation_col,
  main = "Sample Correlation (GSE80655)",
  fontsize = 6,
  silent = TRUE
)

tiff("QC_Sample_Correlation_Heatmap.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")
grid::grid.draw(p_heatmap$gtable)
dev.off()

pdf("QC_Sample_Correlation_Heatmap.pdf", width = 10, height = 8)
grid::grid.draw(p_heatmap$gtable)
dev.off()

cat("✅ Saved: QC_Sample_Correlation_Heatmap.tiff / .pdf\n")


# ------------------------------------------------------------------------------
# 7c. PCA plot (with zero-variance filtering)
# ------------------------------------------------------------------------------
cat("\n--- Generating PCA plot ---\n")

# Remove genes with zero variance (constant across samples)
gene_var <- apply(log_expr, 1, var, na.rm = TRUE)
log_expr_var <- log_expr[gene_var > 0, ]
cat(sprintf("Removed %d genes with zero variance. Remaining: %d\n", 
            sum(gene_var == 0), nrow(log_expr_var)))

# PCA on filtered log2 expression matrix (transposed: samples as rows)
pca_res <- prcomp(t(log_expr_var), scale. = TRUE)

# Extract PC1 and PC2
pca_df <- as.data.frame(pca_res$x[, 1:2])
pca_df$Diagnosis <- sample_metadata_final_subset$Diagnosis
pca_df$Brain_Region <- sample_metadata_final_subset$Brain_Region

# Variance explained
var_explained <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Diagnosis, shape = Brain_Region)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of GSE80655 Samples",
       x = paste0("PC1 (", var_explained[1], "%)"),
       y = paste0("PC2 (", var_explained[2], "%)")) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values = c("Control" = "#2E86AB", "Major_Depression" = "#A23B72"))

save_sci_plot(p_pca, "QC_PCA.tiff", w = 8, h = 6)
cat("✅ Saved: QC_PCA.tiff / .pdf\n")

# ------------------------------------------------------------------------------
# 8. Effect Size (Cohen's d) calculation on log2-transformed data
# ------------------------------------------------------------------------------
cat("\n--- Calculating Cohen's d (per brain region) ---\n")

cohens_d_by_region <- list()
for (region in regions) {
  region_gsm <- regions_gsm_list_final[[region]]
  region_idx <- which(colnames(log_expr) %in% region_gsm)
  region_expr <- log_expr[, region_idx]
  
  region_meta <- sample_metadata_final_subset[
    sample_metadata_final_subset$GSM_ID %in% region_gsm, 
  ]
  region_meta <- region_meta[match(colnames(region_expr), region_meta$GSM_ID), ]
  
  ctrl_idx <- which(region_meta$Diagnosis == "Control")
  mdd_idx <- which(region_meta$Diagnosis == "Major_Depression")
  
  if (length(ctrl_idx) > 0 & length(mdd_idx) > 0) {
    ctrl_mean <- rowMeans(region_expr[, ctrl_idx, drop = FALSE], na.rm = TRUE)
    mdd_mean <- rowMeans(region_expr[, mdd_idx, drop = FALSE], na.rm = TRUE)
    ctrl_sd <- apply(region_expr[, ctrl_idx, drop = FALSE], 1, sd, na.rm = TRUE)
    mdd_sd <- apply(region_expr[, mdd_idx, drop = FALSE], 1, sd, na.rm = TRUE)
    
    cohens_d <- (mdd_mean - ctrl_mean) / sqrt((ctrl_sd^2 + mdd_sd^2) / 2)
    cohens_d_by_region[[region]] <- cohens_d
    
    cat(sprintf("  %s: mean |Cohen's d| = %.4f (median = %.4f)\n",
                region,
                mean(abs(cohens_d), na.rm = TRUE),
                median(abs(cohens_d), na.rm = TRUE)))
  } else {
    cat(sprintf("  %s: insufficient samples for Cohen's d\n", region))
  }
}

# ------------------------------------------------------------------------------
# 9. Save intermediate objects
# ------------------------------------------------------------------------------
saveRDS(expr_data_counts, file = "GSE80655_expr_counts_int.rds")
saveRDS(sample_metadata_final_subset, file = "GSE80655_metadata_clean.rds")
saveRDS(regions_gsm_list_final, file = "GSE80655_regions_list.rds")

cat("✅ Saved intermediate objects:\n")
cat("  - GSE80655_expr_counts_int.rds\n")
cat("  - GSE80655_metadata_clean.rds\n")
cat("  - GSE80655_regions_list.rds\n")

# ------------------------------------------------------------------------------
# 10. Checklist
# ------------------------------------------------------------------------------
cat("\n=== CHECKLIST: 02_preprocess.R ===\n")
cat("✅ Phenotype loaded from GEO (GSE80655)\n")
cat(sprintf("  - Total samples after filtering: %d\n", nrow(sample_metadata_final)))
cat(sprintf("  - MDD: %d, Control: %d\n", 
            sum(sample_metadata_final$Diagnosis == "Major_Depression"),
            sum(sample_metadata_final$Diagnosis == "Control")))
cat(sprintf("  - Brain regions: %s\n", paste(regions, collapse = ", ")))
cat("✅ Batch effect detection completed\n")
if ("batch:ch1" %in% colnames(phenotype_df)) {
  cat("  - ⚠️ Batch column found. See warning above.\n")
} else {
  cat("  - No batch column detected.\n")
}
cat("✅ Expression matrix aligned with metadata\n")
cat(sprintf("  - Final dimensions: %d genes x %d samples\n", 
            nrow(expr_data_counts), ncol(expr_data_counts)))
cat("✅ Converted to integer matrix\n")
cat("✅ QC plots generated and saved:\n")
cat("  - QC_Sample_Correlation_Heatmap.tiff / .pdf\n")
cat("  - QC_PCA.tiff / .pdf\n")
cat("✅ Cohen's d calculated per brain region\n")
for (region in regions) {
  if (!is.null(cohens_d_by_region[[region]])) {
    cat(sprintf("  - %s: mean |d| = %.4f\n", 
                region, mean(abs(cohens_d_by_region[[region]]), na.rm = TRUE)))
  }
}
cat("✅ Intermediate objects saved\n")
cat("=== End of 02_preprocess.R ===\n")

