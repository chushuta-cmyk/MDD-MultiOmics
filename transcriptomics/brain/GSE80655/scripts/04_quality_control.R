# ==============================================================================
# 04_quality_control.R
# DESeq2 DEA per brain region, effect sizes, visualization, and functional enrichment
# ==============================================================================

library(DESeq2)
library(apeglm)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)

source("plot_generator.r")
cat("\n=== GSE80655 Differential Expression & Enrichment Started ===\n")

# ------------------------------------------------------------------------------
# 1. Load objects from previous steps
# ------------------------------------------------------------------------------
# Reload DESeq2 object (pre-normalized) or reconstruct it
expr_counts <- readRDS("GSE80655_expr_counts_int.rds")
metadata <- readRDS("GSE80655_metadata_clean.rds")
regions_list <- readRDS("GSE80655_regions_list.rds")
regions <- names(regions_list)

# VST matrix for effect size and plotting
vsd_matrix <- readRDS("GSE80655_expression_vst.rds")

cat(sprintf("Loaded VST matrix: %d x %d\n", nrow(vsd_matrix), ncol(vsd_matrix)))
cat(sprintf("Regions to analyze: %s\n", paste(regions, collapse = ", ")))

# Save function (already defined in 03, but we keep it here for safety if run standalone)
save_sci_plot <- function(p, f, w = 8, h = 6) {
  ggsave(filename = f, plot = p, device = "tiff", 
         width = w, height = h, dpi = 300, compression = "lzw")
  ggsave(filename = gsub("\\.tiff", ".pdf", f), plot = p, 
         device = "pdf", width = w, height = h)
}

# ------------------------------------------------------------------------------
# 2. Loop over brain regions for DESeq2 DEA
# ------------------------------------------------------------------------------
deg_results_by_region <- list()
cohens_d_by_region <- list()

cat("\n--- Starting DESeq2 DEA per region ---\n")

for (region in regions) {
  cat(paste("\n>>> Processing region:", region, "\n"))
  
  # Subset counts and metadata
  current_gsm <- regions_list[[region]]
  current_counts <- expr_counts[, current_gsm, drop = FALSE]
  current_metadata <- metadata[metadata$GSM_ID %in% current_gsm, ]
  current_metadata <- current_metadata[match(colnames(current_counts), current_metadata$GSM_ID), ]
  
  # Ensure Control is reference
  current_metadata$Diagnosis <- relevel(factor(current_metadata$Diagnosis), ref = "Control")
  
  # Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = current_counts,
    colData = current_metadata,
    design = ~ Diagnosis
# ------------------------------------------------------------------------------
# 2. Loop over brain regions for DESeq2 DEA
# ------------------------------------------------------------------------------
deg_results_by_region <- list()
cohens_d_by_region <- list()

cat("\n--- Starting DESeq2 DEA per region ---\n")

for (region in regions) {
  cat(paste("\n>>> Processing region:", region, "\n"))
  
  # Subset counts and metadata
  current_gsm <- regions_list[[region]]
  current_counts <- expr_counts[, current_gsm, drop = FALSE]
  current_metadata <- metadata[metadata$GSM_ID %in% current_gsm, ]
  current_metadata <- current_metadata[match(colnames(current_counts), current_metadata$GSM_ID), ]
  
  # Ensure Control is reference
  current_metadata$Diagnosis <- relevel(factor(current_metadata$Diagnosis), ref = "Control")
  
  # Build DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = current_counts,
    colData = current_metadata,
    design = ~ Diagnosis
  )
  
  # Filter low counts
  dds <- dds[rowSums(counts(dds)) > 0, ]
  
  # Run DESeq2
  dds <- DESeq(dds, quiet = TRUE)
  
  # Extract results with normal shrinkage (built-in, no apeglm dependency)
  res <- lfcShrink(
    dds, 
    contrast = c("Diagnosis", "Major_Depression", "Control"), 
    type = "normal"
  )
  
  # Format results
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("GeneSymbol") %>%
    dplyr::filter(pvalue < 0.01) %>%
    dplyr::rename(logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj) %>%
    dplyr::mutate(Region = region) %>%
    dplyr::select(GeneSymbol, logFC, P.Value, adj.P.Val, Region)
  
  cat(sprintf("  - DEGs (P < 0.01): %d\n", nrow(res_df)))
  deg_results_by_region[[region]] <- res_df
  
  # --- 2a. Compute Cohen's d for this region using VST data ---
  region_samples <- current_metadata$GSM_ID
  region_vst <- vsd_matrix[, colnames(vsd_matrix) %in% region_samples, drop = FALSE]
  region_vst <- region_vst[, match(region_samples, colnames(region_vst))]
  
  ctrl_idx <- current_metadata$Diagnosis == "Control"
  mdd_idx <- current_metadata$Diagnosis == "Major_Depression"
  
  if (sum(ctrl_idx) > 1 & sum(mdd_idx) > 1) {
    ctrl_expr <- region_vst[, ctrl_idx, drop = FALSE]
    mdd_expr <- region_vst[, mdd_idx, drop = FALSE]
    
    ctrl_mean <- rowMeans(ctrl_expr, na.rm = TRUE)
    mdd_mean <- rowMeans(mdd_expr, na.rm = TRUE)
    ctrl_sd <- apply(ctrl_expr, 1, sd, na.rm = TRUE)
    mdd_sd <- apply(mdd_expr, 1, sd, na.rm = TRUE)
    
    cohens_d <- (mdd_mean - ctrl_mean) / sqrt((ctrl_sd^2 + mdd_sd^2) / 2)
    cohens_d_by_region[[region]] <- cohens_d
    cat(sprintf("  - Mean |Cohen's d|: %.4f\n", mean(abs(cohens_d), na.rm = TRUE)))
  } else {
    cat("  - Skipping Cohen's d: insufficient samples.\n")
  }
}

# ------------------------------------------------------------------------------
# 3. Merge all DEG results
# ------------------------------------------------------------------------------
all_degs <- do.call(rbind, deg_results_by_region)
cat(sprintf("\nTotal merged DEGs (P < 0.01): %d\n", nrow(all_degs)))

# Save merged DEGs
write.csv(all_degs, "GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv", row.names = FALSE)
cat("✅ Saved: GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv\n")

# ==============================================================================
# 4. Visualization: Volcano Plots (Per Region + Merged)
# ==============================================================================
cat("\n--- Generating visualization plots ---\n")

# 4a. volcano map for each brain region
for (region in regions) {
  region_degs <- deg_results_by_region[[region]]
  
  # mark significant genes |logFC| > 0.5 且 P < 0.05
  region_degs$Significant <- ifelse(
    abs(region_degs$logFC) > 0.5 & region_degs$P.Value < 0.05,
    "Significant", "Not Significant"
  )
  
  # count significant genes number
  n_sig <- sum(region_degs$Significant == "Significant")
  
  p <- ggplot(region_degs, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c("Not Significant" = "grey70", "Significant" = "red")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
    labs(
      title = paste("GSE80655", region, "MDD vs Control"),
      subtitle = paste("Significant DEGs (|logFC| > 0.5, P < 0.05):", n_sig),
      x = "log2 Fold Change",
      y = "-log10(P-value)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  save_sci_plot(p, paste0("Volcano_", region, ".tiff"))
  cat(paste("✅ Saved: Volcano_", region, ".tiff / .pdf (", n_sig, " DEGs)\n", sep = ""))
}

# 4b. combine volcano maps show regions in different colors.
all_degs$Region <- factor(all_degs$Region, levels = regions)
all_degs$Significant <- ifelse(
  abs(all_degs$logFC) > 0.5 & all_degs$P.Value < 0.05,
  "Significant", "Not Significant"
)

p_merged <- ggplot(all_degs, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = Region, shape = Significant), alpha = 0.5, size = 1.0) +
  scale_shape_manual(values = c("Not Significant" = 1, "Significant" = 19)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.3) +
  labs(
    title = "GSE80655 MDD vs Control (All Brain Regions)",
    x = "log2 Fold Change",
    y = "-log10(P-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

save_sci_plot(p_merged, "Volcano_GSE80655_Merged.tiff")
cat("✅ Saved: Volcano_GSE80655_Merged.tiff / .pdf\n")

# MA Plot
avg_expr <- rowMeans(vsd_matrix, na.rm = TRUE)
all_degs$AveExpr <- avg_expr[match(all_degs$GeneSymbol, rownames(vsd_matrix))]

p_ma <- ggplot(all_degs, aes(x = AveExpr, y = logFC, color = abs(logFC) > 0.5)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("grey60", "red")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue", alpha = 0.5) +
  labs(title = "MA Plot (GSE80655)", x = "Average Expression (VST)", y = "log2 Fold Change") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
save_sci_plot(p_ma, "MA_GSE80655.tiff")
cat("✅ Saved: MA_GSE80655.tiff / .pdf\n")

# ==============================================================================
# 4d. DEG Heatmap (Per Region) — using plot_generator.r
# ==============================================================================
cat("\n--- Generating DEG heatmaps per region (using plot_generator) ---\n")

for (region in regions) {
  region_degs <- deg_results_by_region[[region]]
  
  if (nrow(region_degs) < 2) {
    cat(sprintf("  - Skipping %s: only %d DEGs\n", region, nrow(region_degs)))
    next
  }
  
  # 1. Add ProbeID column for matching with expression matrix row names
  region_degs$ProbeID <- region_degs$GeneSymbol
  
  # 2. Extract VST expression matrix for this region
  current_gsm <- regions_list[[region]]
  expr_sub <- vsd_matrix[, colnames(vsd_matrix) %in% current_gsm, drop = FALSE]
  expr_sub <- expr_sub[, match(current_gsm, colnames(expr_sub)), drop = FALSE]
  
  # 3. Prepare group factor (named vector)
  current_metadata <- metadata[metadata$GSM_ID %in% current_gsm, ]
  current_metadata <- current_metadata[match(colnames(expr_sub), current_metadata$GSM_ID), ]
  group_vec <- current_metadata$Diagnosis
  names(group_vec) <- current_metadata$GSM_ID
  
  # 4. Call heatmap function from plot_generator.r
  plot_deg_heatmap(
    deg_df = region_degs,
    expr_matrix = expr_sub,
    group_factor = group_vec,
    top_n = 40,
    plot_title = paste("GSE80655 Top 40 DEGs -", region),
    output_prefix = paste0("Heatmap_", region),
    output_dir = "."
  )
}

# ------------------------------------------------------------------------------
# 5. Functional Enrichment (GO, KEGG, Reactome)
# ------------------------------------------------------------------------------
cat("\n--- Running Functional Enrichment ---\n")

genes_for_enrich <- all_degs[all_degs$P.Value < 0.05, "GeneSymbol"]
cat(sprintf("Genes with P < 0.05 for enrichment: %d\n", length(genes_for_enrich)))

if (length(genes_for_enrich) >= 10) {
  entrez_ids <- bitr(genes_for_enrich, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat(sprintf("Mapped to ENTREZID: %d\n", nrow(entrez_ids)))
  
  if (nrow(entrez_ids) > 0) {
    # 5a. GO Biological Process
    cat("\n--- Running GO BP ---\n")
    go <- enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, 
                   ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
    if (!is.null(go) && nrow(go) > 0) {
      write.csv(as.data.frame(go), "GSE80655_GO_BP_Enrichment.csv", row.names = FALSE)
      save_sci_plot(dotplot(go, showCategory = 15), "GSE80655_GO_BP_Dotplot.tiff")
      cat("✅ GO results saved.\n")
    } else {
      cat("No significant GO terms.\n")
    }
    
    # 5b. KEGG
    cat("\n--- Running KEGG ---\n")
    options(timeout = 600)
    if (requireNamespace("KEGG.db", quietly = TRUE)) {
      kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa", 
                         pvalueCutoff = 0.1, qvalueCutoff = 0.2, use_internal_data = TRUE)
    } else {
      kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, organism = "hsa", 
                         pvalueCutoff = 0.1, qvalueCutoff = 0.2)
    }
    if (!is.null(kegg) && nrow(kegg) > 0) {
      kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      write.csv(as.data.frame(kegg_readable), "GSE80655_KEGG_Enrichment.csv", row.names = FALSE)
      save_sci_plot(barplot(kegg, showCategory = 15), "GSE80655_KEGG_Barplot.tiff")
      cat("✅ KEGG results saved.\n")
    } else {
      cat("No significant KEGG pathways.\n")
    }
    
    # 5c. Reactome
    cat("\n--- Running Reactome ---\n")
    react <- enrichPathway(gene = entrez_ids$ENTREZID, organism = "human", 
                           pvalueCutoff = 0.05, readable = TRUE)
    if (!is.null(react) && nrow(react) > 0) {
      write.csv(as.data.frame(react), "GSE80655_Reactome_Enrichment.csv", row.names = FALSE)
      save_sci_plot(dotplot(react, showCategory = 15), "GSE80655_Reactome_Dotplot.tiff")
      cat("✅ Reactome results saved.\n")
    } else {
      cat("No significant Reactome pathways.\n")
    }
  }
} else {
  cat("Skipping enrichment: too few genes (<10).\n")
}

# ------------------------------------------------------------------------------
# 6. Final Checklist
# ------------------------------------------------------------------------------
cat("\n=== CHECKLIST: 04_quality_control.R ===\n")
cat("✅ DESeq2 DEA completed for all brain regions\n")
for (r in regions) {
  cat(sprintf("  - %s: %d DEGs (P<0.01)\n", r, nrow(deg_results_by_region[[r]])))
}
cat("✅ Merged DEG table saved\n")
cat("✅ Cohen's d computed (per region, using VST data)\n")
cat("✅ Visualization plots saved: Volcano, MA\n")
cat("✅ Functional Enrichment (GO/KEGG/Reactome) attempted with P<0.05 genes\n")
cat("  - Check log files for significant terms.\n")
cat("=== End of 04_quality_control.R ===\n")