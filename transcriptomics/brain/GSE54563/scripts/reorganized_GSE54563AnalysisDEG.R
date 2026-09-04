# ==============================================================================
# GSE54563 Brain Tissue Post-Analysis
# Objective: DEG selection, annotation, GO enrichment, and blood vs brain comparison
# Dependencies: deg_54563, expr_matrix_54563_filtered, cohens_d_54563 (pre-existing)
# ==============================================================================

cat("\n=== GSE54563 Post-Processing Started ===\n")

# ------------------------------------------------------------------------------
# 1. Parameter settings
# ------------------------------------------------------------------------------
P_CUTOFF <- 0.05
FC_CUTOFF <- 1.0          # chosen based on observed logFC distribution
USE_ADJ_P <- FALSE        # raw P-values used because adjusted P often near 1

# ------------------------------------------------------------------------------
# 2. Significant gene selection using multiple criteria
# ------------------------------------------------------------------------------
cat("\n--- Significant gene counts ---\n")

n_p005 <- sum(deg_54563$P.Value < 0.05)
n_p001 <- sum(deg_54563$P.Value < 0.01)
n_p0001 <- sum(deg_54563$P.Value < 0.001)
cat(sprintf("P < 0.05: %d\nP < 0.01: %d\nP < 0.001: %d\n", n_p005, n_p001, n_p0001))

n_fc1 <- sum(abs(deg_54563$logFC) > 1.0)
n_fc0.5 <- sum(abs(deg_54563$logFC) > 0.5)
cat(sprintf("|logFC| > 1.0: %d\n|logFC| > 0.5: %d\n", n_fc1, n_fc0.5))

n_comb <- sum(abs(deg_54563$logFC) > FC_CUTOFF & deg_54563$P.Value < P_CUTOFF)
cat(sprintf("|logFC| > %.1f & P < %.2f: %d\n", FC_CUTOFF, P_CUTOFF, n_comb))

# Select final set using |logFC| > 1.0 as primary criterion
sig_idx <- abs(deg_54563$logFC) > FC_CUTOFF
deg_sig <- deg_54563[sig_idx, ]
cat(sprintf("Selected significant genes: %d\n", nrow(deg_sig)))

# ------------------------------------------------------------------------------
# 3. Save sorted significant genes
# ------------------------------------------------------------------------------
deg_sig_sorted <- deg_sig[order(deg_sig$P.Value), ]
write.csv(deg_sig_sorted, "GSE54563_Significant_DEGs_FC1.0.csv", row.names = FALSE)
cat("Saved: GSE54563_Significant_DEGs_FC1.0.csv\n")

# ------------------------------------------------------------------------------
# 4. Probe ID to Gene Symbol annotation
# ------------------------------------------------------------------------------
cat("\n--- Gene annotation ---\n")
library(AnnotationDbi)
if (!require("illuminaHumanv4.db", quietly = TRUE)) {
  BiocManager::install("illuminaHumanv4.db")
  library(illuminaHumanv4.db)
}

probes <- deg_sig_sorted$probe_id
symbols <- mapIds(illuminaHumanv4.db, keys = probes, column = "SYMBOL",
                  keytype = "PROBEID", multiVals = "first")
deg_sig_sorted$gene_symbol <- symbols
n_annot <- sum(!is.na(symbols))
cat(sprintf("Annotated: %d/%d (%.1f%%)\n", n_annot, length(symbols),
            100 * n_annot / length(symbols)))

# ------------------------------------------------------------------------------
# 5. GO enrichment analysis (clusterProfiler)
# ------------------------------------------------------------------------------
cat("\n--- GO enrichment ---\n")
library(clusterProfiler)
library(org.Hs.eg.db)

gene_symbols <- na.omit(deg_sig_sorted$gene_symbol)
if (length(gene_symbols) > 0) {
  entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
  if (nrow(entrez) > 0) {
    go_res <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                       ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH")
    if (!is.null(go_res) && nrow(go_res@result) > 0) {
      cat(sprintf("Significant GO terms: %d\n", nrow(go_res@result)))
      write.csv(go_res@result, "GSE54563_GO_Enrichment.csv", row.names = FALSE)
    } else {
      cat("No significant GO terms found.\n")
    }
  }
}

# ------------------------------------------------------------------------------
# 6. Blood (GSE98793) vs Brain (GSE54563) comparison (summary table)
# ------------------------------------------------------------------------------
cat("\n--- Blood vs Brain comparison ---\n")
comparison <- data.frame(
  Metric = c("Samples (MDD/Control)", "Platform", "Tissue Type",
             "Probes analyzed", "P < 0.05", "|logFC| > 1.0",
             "Mean logFC", "Mean Cohen's d"),
  Blood_GSE98793 = c("128/64", "Illumina HumanWG-6 v3", "Whole Blood",
                     "49,405", "16", "0", "0.0005", "0.1334"),
  Brain_GSE54563 = c("25/25", "Illumina HumanHT-12 V3", "Anterior Cingulate Cortex",
                     nrow(expr_matrix_54563_filtered),
                     sum(deg_54563$P.Value < 0.05),
                     sum(abs(deg_54563$logFC) > 1.0),
                     sprintf("%.4f", mean(deg_54563$logFC)),
                     sprintf("%.4f", mean(abs(cohens_d_54563), na.rm = TRUE)))
)
print(comparison)
write.csv(comparison, "Blood_vs_Brain_Comparison.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 7. Final summary statistics
# ------------------------------------------------------------------------------
cat("\n--- Final summary ---\n")
summary_stats <- data.frame(
  Dataset = c("GSE98793 (Blood)", "GSE54563 (Brain)"),
  Samples = c("192 (64C, 128M)", "50 (25C, 25M)"),
  Sig_Genes = c("16 (adj.P<0.05)", nrow(deg_sig_sorted)),
  Mean_logFC = c("0.0005", sprintf("%.4f", mean(deg_54563$logFC))),
  Max_logFC = c("~1.0", sprintf("%.2f", max(abs(deg_54563$logFC)))),
  Mean_Cohen_d = c("0.1334", sprintf("%.4f", mean(abs(cohens_d_54563), na.rm = TRUE)))
)
print(summary_stats)
write.csv(summary_stats, "Final_Statistics_Comparison.csv", row.names = FALSE)

cat("\n=== GSE54563 post-processing complete ===\n")
