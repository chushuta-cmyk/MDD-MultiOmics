# ==============================================================================
# GSE54563 Visualization – Brain Tissue DEG Results (CORRECTED)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Setup & Reproducibility
# ------------------------------------------------------------------------------
set.seed(123)  # for reproducibility (VennDiagram / pheatmap)

library(ggplot2)
library(dplyr)
library(VennDiagram)
library(pheatmap)   # added for heatmap

out_dir <- "."  # change to "figures/" if needed
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Volcano Plot
# ------------------------------------------------------------------------------
volcano_data <- deg_54563 %>%
  mutate(
    threshold = case_when(
      P.Value < 0.01 & abs(cohens_d) > 0.5 ~ "Significant",
      P.Value < 0.05 ~ "Moderate",
      TRUE ~ "Not significant"
    )
  )

p_volcano <- ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value), color = threshold)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Significant"   = "#FF1493",
      "Moderate"      = "#FFB6C1",
      "Not significant" = "#CCCCCC"
    )
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "red") +
  labs(
    title = "Brain Tissue (GSE54563) – Volcano Plot",
    x = "log2(Fold Change)",
    y = "-log10(P‑value)",
    color = "Gene Category"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, "Volcano_Brain_Tissue.pdf"), p_volcano, width = 10, height = 8)

# ------------------------------------------------------------------------------
# 2. MA Plot – CORRECTED (uses actual mean expression, not log2(exp(logFC)))
# ------------------------------------------------------------------------------
if (exists("expr_filtered") && !is.null(expr_filtered)) {
  # Match expression rows to DEG table
  common_genes <- intersect(rownames(expr_filtered), rownames(deg_54563))
  avg_expr <- rowMeans(expr_filtered[common_genes, ], na.rm = TRUE)
  
  ma_data <- data.frame(
    logFC = deg_54563[common_genes, "logFC"],
    A = avg_expr
  )
  
  p_ma <- ggplot(ma_data, aes(x = logFC, y = A)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = "MA Plot: Brain Tissue Expression",
      x = "log2(Fold Change)",
      y = "log2(Mean Expression)"
    ) +
    theme_minimal()
  
  ggsave(file.path(out_dir, "MA_Plot_Brain.pdf"), p_ma, width = 10, height = 8)
} else {
  message("⚠️  expr_filtered not found. Skipping MA plot.")
}

# ------------------------------------------------------------------------------
# 3. Effect Size Distribution – FIXED labels using annotate()
# ------------------------------------------------------------------------------
effect_data <- data.frame(cohens_d = abs(cohens_d_54563))

p_effect <- ggplot(effect_data, aes(x = cohens_d)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "orange", size = 0.8) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "green", size = 0.8) +
  annotate("text", x = 0.25, y = Inf, label = "Small (0.2)", vjust = 2, color = "red", size = 4) +
  annotate("text", x = 0.55, y = Inf, label = "Medium (0.5)", vjust = 2, color = "orange", size = 4) +
  annotate("text", x = 0.85, y = Inf, label = "Large (0.8)", vjust = 2, color = "green", size = 4) +
  labs(
    title = "Distribution of Effect Sizes (Cohen's d)",
    x = "Absolute Cohen's d",
    y = "Gene Count"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "Effect_Size_Distribution_Brain.pdf"), p_effect, width = 10, height = 6)

# ------------------------------------------------------------------------------
# 4. Venn Diagram – DYNAMIC loading of blood genes + robust brain gene extraction
# ------------------------------------------------------------------------------
# Safely extract brain genes from DEG table
if (!is.null(deg_54563$gene_symbol)) {
  brain_genes <- deg_54563$gene_symbol[deg_54563$adj.P.Val < 0.05 & !is.na(deg_54563$gene_symbol)]
} else {
  # Fallback: use probe IDs if gene_symbol missing
  brain_genes <- rownames(deg_54563)[deg_54563$adj.P.Val < 0.05]
}

# Load actual blood results (if available)
if (file.exists("GSE98793_Blood_DEG_Results.csv")) {
  blood_res <- read.csv("GSE98793_Blood_DEG_Results.csv")
  blood_genes <- blood_res$gene_symbol[blood_res$adj.P.Val < 0.05 & !is.na(blood_res$gene_symbol)]
  blood_mean_logFC <- mean(blood_res$logFC, na.rm = TRUE)
  
  if (length(blood_genes) > 0 && length(brain_genes) > 0) {
    venn.diagram(
      list(Blood = blood_genes, Brain = brain_genes),
      filename = file.path(out_dir, "Venn_Blood_vs_Brain.pdf"),
      imagetype = "pdf",
      main = "Common Genes: Blood vs Brain"
    )
    cat("✓ Venn diagram saved.\n")
  } else {
    message("⚠️  Insufficient genes for Venn diagram (blood or brain list empty).")
  }
} else {
  message("⚠️  GSE98793_Blood_DEG_Results.csv not found. Skipping Venn diagram.")
}

# ------------------------------------------------------------------------------
# 5. Mean logFC Comparison (Blood vs Brain) – DYNAMIC
# ------------------------------------------------------------------------------
if (exists("blood_mean_logFC")) {
  comparison_df <- data.frame(
    logFC = c(blood_mean_logFC, mean(deg_54563$logFC[deg_54563$P.Value < 0.05], na.rm = TRUE)),
    Dataset = c("Blood (GSE98793)", "Brain (GSE54563)")
  )
  
  p_comparison <- ggplot(comparison_df, aes(x = Dataset, y = logFC, fill = Dataset)) +
    geom_bar(stat = "identity", alpha = 0.7) +
    labs(
      title = "Mean logFC Comparison: Blood vs Brain",
      y = "Mean logFC"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(out_dir, "Comparison_logFC_Blood_vs_Brain.pdf"), p_comparison, width = 8, height = 6)
} else {
  message("⚠️  Skipping blood vs brain bar plot (blood data not loaded).")
}

# ------------------------------------------------------------------------------
# 6. P-value Distribution
# ------------------------------------------------------------------------------
pval_data <- data.frame(minus_log10_p = -log10(deg_54563$P.Value))

p_pval <- ggplot(pval_data, aes(x = minus_log10_p)) +
  geom_histogram(bins = 50, fill = "darkblue", alpha = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
  geom_vline(xintercept = -log10(0.01), linetype = "dashed", color = "orange", size = 1) +
  labs(
    title = "Distribution of -log10(P‑values)",
    x = "-log10(P‑value)",
    y = "Gene Count"
  ) +
  theme_minimal()

ggsave(file.path(out_dir, "Pvalue_Distribution_Brain.pdf"), p_pval, width = 10, height = 6)

# ------------------------------------------------------------------------------
# 7. NEW: Heatmap of Top 50 DEGs (Reviewer expectation)
# ------------------------------------------------------------------------------
if (exists("expr_filtered") && nrow(expr_filtered) > 50) {
  # Get top 50 genes by P-value
  top50_genes <- rownames(deg_54563[order(deg_54563$P.Value), ])[1:50]
  top50_expr <- expr_filtered[rownames(expr_filtered) %in% top50_genes, ]
  
  # Match order
  top50_expr <- top50_expr[match(top50_genes, rownames(top50_expr)), ]
  
  # Scale rows for visualization
  heatmap_data <- t(scale(t(top50_expr)))
  
  # Prepare annotation
  if (exists("pheno_54563")) {
    annotation_col <- data.frame(
      Diagnosis = pheno_54563$disease,
      Pair = as.factor(pheno_54563$pair)
    )
    rownames(annotation_col) <- colnames(heatmap_data)
  } else {
    annotation_col <- NULL
  }
  
  pheatmap(
    heatmap_data,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = "Top 50 DEGs (GSE54563)",
    filename = file.path(out_dir, "Heatmap_Top50_DEGs.pdf"),
    width = 10,
    height = 8
  )
  cat("✓ Heatmap saved.\n")
} else {
  message("⚠️  expr_filtered missing or insufficient genes. Skipping heatmap.")
}

# ------------------------------------------------------------------------------
cat("\n✅ All corrected visualizations saved to:", out_dir, "\n")
# ------------------------------------------------------------------------------
