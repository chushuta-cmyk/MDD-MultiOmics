# ==============================================================================
# GSE98793 Blood MDD - Immune Infiltration Analysis (immunedeconv)
# ==============================================================================

cat("\n=== GSE98793 Immune Infiltration Pipeline Started ===\n")

# ------------------------------------------------------------------------------
# 1. Setup and data preparation
# ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(limma)
library(pheatmap)

# Assume expr_matrix_98793 and deg_98793 are pre-loaded
if (!exists("expr_matrix_98793")) stop("expr_matrix_98793 not found")
if (!exists("deg_98793")) stop("deg_98793 not found")

cat("GSE98793 expression matrix:", nrow(expr_matrix_98793), "genes x", 
    ncol(expr_matrix_98793), "samples\n")

# Significant genes from blood (for reference)
sig_genes_blood <- deg_98793[deg_98793$P.Value < 0.05 & abs(deg_98793$logFC) > 0.5, ]
cat("Significant blood genes (P<0.05, |logFC|>0.5):", nrow(sig_genes_blood), "\n")

# ------------------------------------------------------------------------------
# 2. Install/Load immunedeconv
# ------------------------------------------------------------------------------
if (!require("immunedeconv", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("omnideconv/immunedeconv")
}
library(immunedeconv)
cat("immunedeconv loaded.\n")

# ------------------------------------------------------------------------------
# 3. Prepare expression matrix (use all genes for accurate deconvolution)
# ------------------------------------------------------------------------------
expr_for_immune <- expr_matrix_98793
cat("Expression matrix for deconvolution:", nrow(expr_for_immune), "x", 
    ncol(expr_for_immune), "\n")

# ------------------------------------------------------------------------------
# 4. Run multiple deconvolution algorithms (quantiseq, timer, epic, xcell)
# ------------------------------------------------------------------------------
methods_to_run <- c("quantiseq", "timer", "epic", "xcell")
immune_results <- list()

for (m in methods_to_run) {
  cat("Running", m, "... ")
  result <- tryCatch(
    deconvolute(expr_for_immune, method = m),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    immune_results[[m]] <- result
    cat("OK\n")
  } else {
    cat("FAILED\n")
  }
}

if (length(immune_results) == 0) stop("No deconvolution method succeeded.")

# Use the first successful result for downstream analysis
immune_main <- immune_results[[1]]

# Convert to long format for plotting
immune_df <- immune_main %>%
  as.data.frame() %>%
  tibble::column_to_rownames("cell_type")

cell_types <- rownames(immune_df)
samples <- colnames(immune_df)

# Assume we have phenotype groups defined elsewhere; if not, derive from sample names or metadata.
# For simplicity, we assume a vector `mdd_samples` exists indicating MDD samples.
# If not, create a placeholder (user should supply real grouping).
if (!exists("mdd_samples")) {
  warning("mdd_samples not defined. Using first half as MDD for demo.")
  mdd_samples <- samples[1:floor(length(samples)/2)]
}

immune_long <- immune_df %>%
  rownames_to_column("Cell_type") %>%
  pivot_longer(-Cell_type, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Group = ifelse(Sample %in% mdd_samples, "MDD", "Control"))

cat("Immune data ready:", length(unique(immune_long$Cell_type)), "cell types,", 
    nrow(immune_long), "observations\n")

# ------------------------------------------------------------------------------
# 5. Statistical comparison: MDD vs Control per cell type
# ------------------------------------------------------------------------------
immune_comparison <- immune_long %>%
  group_by(Cell_type) %>%
  summarise(
    MDD_mean   = mean(Abundance[Group == "MDD"]),
    Control_mean = mean(Abundance[Group == "Control"]),
    P_value    = t.test(Abundance[Group == "MDD"], Abundance[Group == "Control"])$p.value,
    FC         = MDD_mean / Control_mean,
    .groups = "drop"
  ) %>%
  mutate(
    Sig = ifelse(P_value < 0.05, "***", "ns"),
    Direction = ifelse(FC > 1, "↑ in MDD", "↓ in MDD")
  ) %>%
  arrange(P_value)

cat("\nImmune cell types with significant differences (P<0.05):\n")
print(filter(immune_comparison, P_value < 0.05))

# ------------------------------------------------------------------------------
# 6. Save results and generate plots
# ------------------------------------------------------------------------------
dir.create("immune_results", showWarnings = FALSE)
write.csv(immune_comparison, "immune_results/Immune_Cell_Comparison.csv", row.names = FALSE)

# Boxplot
p_box <- ggplot(immune_long, aes(x = Cell_type, y = Abundance, fill = Group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  coord_flip() +
  labs(title = "Immune Cell Infiltration: MDD vs Control",
       x = "Cell Type", y = "Abundance Score") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("immune_results/Immune_Cell_Boxplot.pdf", p_box, width = 10, height = 8)

# Heatmap
sample_anno <- data.frame(
  Group = ifelse(colnames(immune_df) %in% mdd_samples, "MDD", "Control"),
  row.names = colnames(immune_df)
)
pheatmap(immune_df, annotation_col = sample_anno,
         main = "Immune Cell Infiltration Heatmap",
         scale = "row", filename = "immune_results/Immune_Cell_Heatmap.pdf",
         width = 12, height = 8)

# FC barplot
p_fc <- ggplot(immune_comparison, aes(x = reorder(Cell_type, FC), y = log2(FC), fill = Direction)) +
  geom_bar(stat = "identity", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = Sig), vjust = ifelse(log2(FC) > 0, -0.5, 1.5)) +
  coord_flip() +
  labs(title = "log2(Fold Change) in Immune Cells: MDD/Control",
       x = "Cell Type", y = "log2(Fold Change)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("immune_results/Immune_Cell_FC_Plot.pdf", p_fc, width = 10, height = 8)

# ------------------------------------------------------------------------------
# 7. Generate summary report
# ------------------------------------------------------------------------------
sig_cells <- filter(immune_comparison, P_value < 0.05)
summary_text <- sprintf("
GSE98793 Immune Infiltration Summary
=====================================
Expression matrix: %d genes x %d samples
MDD: %d, Control: %d
Deconvolution method: %s

Significant immune cell types (P<0.05):
%s

Key findings:
- %d cell types show differential infiltration between MDD and Control.
- Top altered: %s

Implications:
- Systemic immune dysregulation in MDD.
- Potential biomarkers and therapeutic targets.
",
nrow(expr_for_immune), ncol(expr_for_immune),
sum(immune_long$Group == "MDD") / length(unique(immune_long$Cell_type)),
sum(immune_long$Group == "Control") / length(unique(immune_long$Cell_type)),
names(immune_results)[1],
paste0(sig_cells$Cell_type, " (FC=", round(sig_cells$FC,2), ", P=", format(sig_cells$P_value, digits=3), ")", collapse="; "),
nrow(sig_cells),
paste(head(sig_cells$Cell_type, 3), collapse=", ")
)
writeLines(summary_text, "immune_results/Immune_Analysis_Summary.txt")

cat("\n=== Immune infiltration analysis complete ===\n")
cat("Output files saved in 'immune_results/' directory.\n")
