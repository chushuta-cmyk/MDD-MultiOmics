# ==============================================================================
# GSE98793 Complete Analysis Pipeline (Local File Version)
# ==============================================================================
# Objective:
#   1. Read expression matrix and phenotype from local series_matrix.txt
#   2. QC filtering (remove NA, low expression)
#   3. Differential expression analysis (limma)
#   4. Probe ID to Gene Symbol annotation
#   5. Visualization: volcano, heatmap, PCA
#   6. Immune-related gene extraction and biomarker scoring
# ==============================================================================

cat("\n=== GSE98793 Analysis Pipeline (Local) Started ===\n")

# ------------------------------------------------------------------------------
# 1. Setup: load libraries and set working directory
# ------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(matrixStats)
library(limma)
library(AnnotationDbi)
library(hgu133plus2.db)   # platform-specific annotation

setwd("/Users/tan/Developer/projects/r_project/Study/GSE98793/GSE98793_analysis/")

# ------------------------------------------------------------------------------
# 2. Load expression matrix (skip 78 lines)
# ------------------------------------------------------------------------------
cat("\n--- Loading expression matrix ---\n")
expr_raw <- fread(
  "GSE98793_series_matrix.txt",
  skip = 78,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  fill = TRUE
)
cat(sprintf("Raw data: %d rows x %d columns\n", nrow(expr_raw), ncol(expr_raw)))

# Convert to numeric matrix with probe IDs as rownames
probe_ids <- expr_raw[, 1]
expr_matrix <- as.matrix(expr_raw[, -1])
rownames(expr_matrix) <- probe_ids
storage.mode(expr_matrix) <- "numeric"
cat(sprintf("Expression matrix: %d probes x %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))

# ------------------------------------------------------------------------------
# 3. Load phenotype information (lines 41-42: title and GEO accession)
# ------------------------------------------------------------------------------
cat("\n--- Loading phenotype ---\n")
pheno_raw <- fread(
  "GSE98793_series_matrix.txt",
  skip = 41,
  sep = "\t",
  header = FALSE,
  nrows = 2,
  data.table = FALSE,
  fill = TRUE
)

# Extract and clean sample titles and GEO IDs
sample_titles <- as.character(pheno_raw[1, -1])
sample_geo_ids <- as.character(pheno_raw[2, -1])
sample_titles <- gsub('"', '', sample_titles)
sample_geo_ids <- gsub('"', '', sample_geo_ids)

# Build phenotype data frame
pheno_data <- data.frame(
  Sample_ID   = sample_geo_ids,
  Sample_Title = sample_titles,
  stringsAsFactors = FALSE
)

# Assign disease group based on title pattern
pheno_data$Disease_Group <- ifelse(
  grepl("control", pheno_data$Sample_Title, ignore.case = TRUE),
  "Control", "MDD"
)

cat(sprintf("Phenotype: %d samples (Control: %d, MDD: %d)\n",
            nrow(pheno_data),
            sum(pheno_data$Disease_Group == "Control"),
            sum(pheno_data$Disease_Group == "MDD")))

# ------------------------------------------------------------------------------
# 4. Align expression columns to phenotype order
# ------------------------------------------------------------------------------
cat("\n--- Aligning expression and phenotype ---\n")
expr_matrix <- expr_matrix[, match(pheno_data$Sample_ID, colnames(expr_matrix))]
cat("Alignment complete.\n")

# ------------------------------------------------------------------------------
# 5. Quality control: remove NA rows and low-expression probes
# ------------------------------------------------------------------------------
cat("\n--- Quality control ---\n")
# Remove rows with any NA
na_per_row <- rowSums(is.na(expr_matrix))
keep_no_na <- na_per_row == 0
expr_clean <- expr_matrix[keep_no_na, ]
cat(sprintf("Rows with NA removed: %d -> %d\n", sum(!keep_no_na), nrow(expr_clean)))

# Keep probes expressed in at least 25% of samples (threshold: log2(8)=3)
keep_expr <- rowSums(expr_clean > log2(8)) >= 0.25 * ncol(expr_clean)
expr_filtered <- expr_clean[keep_expr, ]
cat(sprintf("Low-expression probes removed: %d -> %d\n",
            sum(!keep_expr), nrow(expr_filtered)))

# ------------------------------------------------------------------------------
# 6. Differential expression analysis (limma)
# ------------------------------------------------------------------------------
cat("\n--- Differential expression analysis ---\n")
group <- factor(pheno_data$Disease_Group, levels = c("Control", "MDD"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(expr_filtered, design)
contrast <- makeContrasts(MDD_vs_Control = MDD - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

deg_raw <- topTable(fit2, coef = "MDD_vs_Control", number = Inf, adjust.method = "BH")
deg_raw$probe_id <- rownames(deg_raw)
cat(sprintf("DEG analysis complete. %d probes tested.\n", nrow(deg_raw)))

# ------------------------------------------------------------------------------
# 7. Gene annotation (Probe ID → Gene Symbol)
# ------------------------------------------------------------------------------
cat("\n--- Gene annotation ---\n")
probes <- rownames(deg_raw)
symbols <- mapIds(hgu133plus2.db, keys = probes, column = "SYMBOL",
                  keytype = "PROBEID", multiVals = "first")
deg_raw$gene_symbol <- as.character(symbols)

# Remove genes without annotation
deg_annotated <- deg_raw[!is.na(deg_raw$gene_symbol), ]
# Deduplicate: keep the probe with highest absolute logFC per gene
deg_annotated <- deg_annotated[order(abs(deg_annotated$logFC), decreasing = TRUE), ]
deg_unique <- deg_annotated[!duplicated(deg_annotated$gene_symbol), ]

cat(sprintf("Annotated unique genes: %d / %d (%.1f%%)\n",
            nrow(deg_unique), nrow(deg_raw),
            100 * nrow(deg_unique) / nrow(deg_raw)))

# ------------------------------------------------------------------------------
# 8. Select significant DEGs (adj.P < 0.05, |logFC| > 0.5)
# ------------------------------------------------------------------------------
deg_sig <- deg_unique[deg_unique$adj.P.Val < 0.05 & abs(deg_unique$logFC) > 0.5, ]
deg_sig <- deg_sig[order(abs(deg_sig$logFC), decreasing = TRUE), ]
cat(sprintf("Significant DEGs: %d\n", nrow(deg_sig)))

# ------------------------------------------------------------------------------
# 9. Volcano plot
# ------------------------------------------------------------------------------
cat("\n--- Volcano plot ---\n")
volcano_data <- deg_unique %>%
  mutate(Color = case_when(
    adj.P.Val < 0.05 & logFC > 0.5 ~ "Up-regulated",
    adj.P.Val < 0.05 & logFC < -0.5 ~ "Down-regulated",
    TRUE ~ "Not significant"
  ))

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = Color)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Up-regulated" = "#FF6B6B",
                                "Down-regulated" = "#4ECDC4",
                                "Not significant" = "#CCCCCC")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(title = "GSE98793: MDD vs Control", x = "log2(Fold Change)", y = "-log10(adj.P)") +
  theme_minimal()
ggsave("GSE98793_volcano.pdf", p_volcano, width = 10, height = 8)

# ------------------------------------------------------------------------------
# 10. Heatmap of top 50 DEGs
# ------------------------------------------------------------------------------
cat("\n--- Heatmap (Top 50 DEGs) ---\n")
if (nrow(deg_sig) > 0) {
  top_genes <- rownames(deg_sig)[1:min(50, nrow(deg_sig))]
  top_expr <- expr_filtered[top_genes, ]
  rownames(top_expr) <- deg_unique[rownames(top_expr), "gene_symbol"]
  top_scaled <- t(scale(t(top_expr)))
  
  anno_col <- data.frame(
    Group = pheno_data$Disease_Group,
    row.names = colnames(top_scaled)
  )
  
  pheatmap(top_scaled,
           annotation_col = anno_col,
           show_colnames = FALSE,
           main = "Top 50 DEGs (MDD vs Control)",
           filename = "GSE98793_heatmap.pdf",
           width = 12, height = 14)
}

# ------------------------------------------------------------------------------
# 11. Immune-related gene extraction
# ------------------------------------------------------------------------------
cat("\n--- Immune-related genes ---\n")
immune_genes <- c(
  "CD4", "CD8A", "CD27", "CD28", "IL7R", "FOXP3",
  "CD19", "CD79A", "CD14", "CD16", "CSF1R",
  "IL1B", "IL6", "TNF", "IL10", "IL12A", "CXCL9", "CXCL10",
  "PDCD1", "LAG3", "CTLA4", "TIGIT", "CD274"
)
immune_deg <- deg_unique[deg_unique$gene_symbol %in% immune_genes, ]
immune_deg_sig <- immune_deg[immune_deg$adj.P.Val < 0.05, ]
cat(sprintf("Immune-related DEGs (adj.P<0.05): %d\n", nrow(immune_deg_sig)))

# ------------------------------------------------------------------------------
# 12. Biomarker candidate scoring
# ------------------------------------------------------------------------------
cat("\n--- Biomarker candidate scoring ---\n")
biomarkers <- deg_sig %>%
  mutate(
    is_immune = gene_symbol %in% immune_genes,
    score = -log10(adj.P.Val) * abs(logFC)
  ) %>%
  arrange(desc(score))

cat("Top 20 biomarker candidates:\n")
print(head(biomarkers[, c("gene_symbol", "logFC", "adj.P.Val", "score")], 20))

# ------------------------------------------------------------------------------
# 13. PCA plot using significant DEGs
# ------------------------------------------------------------------------------
cat("\n--- PCA plot ---\n")
if (nrow(deg_sig) > 1) {
  pca <- prcomp(t(expr_filtered[rownames(deg_sig), ]), scale. = TRUE)
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = pheno_data$Disease_Group
  )
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = "PCA plot (significant DEGs)") +
    theme_minimal()
  ggsave("GSE98793_pca.pdf", p_pca, width = 10, height = 8)
}

# ------------------------------------------------------------------------------
# 14. Export results
# ------------------------------------------------------------------------------
cat("\n--- Exporting results ---\n")
write.csv(deg_unique, "GSE98793_all_DEGs.csv", row.names = FALSE)
write.csv(deg_sig, "GSE98793_sig_DEGs.csv", row.names = FALSE)
write.csv(immune_deg_sig, "GSE98793_immune_DEGs.csv", row.names = FALSE)
write.csv(biomarkers, "GSE98793_biomarkers.csv", row.names = FALSE)
write.csv(pheno_data, "GSE98793_pheno.csv", row.names = FALSE)

# Save R objects for downstream use
save(expr_filtered, pheno_data, deg_unique, deg_sig, biomarkers,
     file = "GSE98793_analysis_objects.RData")

cat("\n=== Analysis complete. All outputs saved. ===\n")
