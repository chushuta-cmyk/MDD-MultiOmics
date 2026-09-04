# ==============================================================================
# GSE54563 Brain Tissue Analysis – Complete Pipeline (Part 1)
# Objective: Load data, parse paired metadata, run paired limma DEA
# ==============================================================================

cat("\n=== GSE54563 Full Pipeline Started ===\n")

# ------------------------------------------------------------------------------
# 0. Setup and libraries
# ------------------------------------------------------------------------------
library(data.table)
library(dplyr)
library(limma)
library(matrixStats)
library(ggplot2)

save_sci_plot <- function(p, f, w=8, h=6) {
  ggsave(filename = f, plot = p, device = "tiff", width = w, height = h, dpi = 300, compression = "lzw")
  ggsave(filename = gsub("\\.tiff", ".pdf", f), plot = p, device = "pdf", width = w, height = h)
}
# Set working directory (modify as needed)
# setwd("/path/to/GSE54563/")

# ------------------------------------------------------------------------------
# 1. Parse phenotype from Lines 29-30 (Sample Title & GEO ID)
# ------------------------------------------------------------------------------
cat("\n--- 1. Parsing phenotype (Lines 29-30) ---\n")

pheno_raw <- fread(
  "GSE54563_series_matrix.txt",
  skip = 28,          # Line 29
  sep = "\t",
  header = FALSE,
  nrows = 2,          # Lines 29 and 30
  data.table = FALSE,
  fill = TRUE
)

# Extract and clean sample titles and GEO IDs
sample_titles <- gsub('"', '', as.character(pheno_raw[1, -1]))
sample_geo_ids <- gsub('"', '', as.character(pheno_raw[2, -1]))

# Assign disease group based on "MDD matched" vs "Control matched"
disease_group <- ifelse(
  grepl("^MDD matched", sample_titles, ignore.case = TRUE),
  "MDD", "Control"
)

# Extract pair numbers (e.g., from "MDD matched #1")
pair_number <- as.numeric(gsub(".*#(\\d+)$", "\\1", sample_titles))

# Build metadata
pheno_54563 <- data.frame(
  sample_id = sample_geo_ids,
  title = sample_titles,
  disease = disease_group,
  pair = pair_number,
  stringsAsFactors = FALSE
)

cat(sprintf("Metadata: %d samples (%d MDD, %d Control, %d pairs)\n",
            nrow(pheno_54563),
            sum(pheno_54563$disease == "MDD"),
            sum(pheno_54563$disease == "Control"),
            length(unique(pheno_54563$pair))))

# ------------------------------------------------------------------------------
# 2. Load expression matrix (from Line 67 onwards)
# ------------------------------------------------------------------------------
cat("\n--- 2. Loading expression matrix ---\n")

expr_raw <- fread(
  "GSE54563_series_matrix.txt",
  skip = 66,           # Line 67
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  fill = TRUE
)

# Convert to matrix with probe IDs as rownames
probe_ids <- expr_raw[, 1]
expr_matrix <- as.matrix(expr_raw[, -1])
rownames(expr_matrix) <- probe_ids

# Enforce numeric storage (efficient)
storage.mode(expr_matrix) <- "numeric"

# Align expression columns to phenotype order
expr_matrix <- expr_matrix[, match(pheno_54563$sample_id, colnames(expr_matrix))]
cat(sprintf("Expression matrix: %d probes x %d samples\n",
            nrow(expr_matrix), ncol(expr_matrix)))

# ------------------------------------------------------------------------------
# 3. Quality control: remove NAs and low-expressed probes
# ------------------------------------------------------------------------------
cat("\n--- 3. Quality control ---\n")

# Remove rows with any NA (or impute, but here we remove)
na_rows <- rowSums(is.na(expr_matrix)) > 0
expr_clean <- expr_matrix[!na_rows, ]
cat(sprintf("Rows with NA removed: %d, remaining: %d\n",
            sum(na_rows), nrow(expr_clean)))

max_val <- max(expr_clean, na.rm = TRUE)
if (max_val > 50) {
  cat(sprintf("Raw scale detected (max=%.2f), applying log2(x+1)...\n", max_val))
  expr_clean <- log2(expr_clean + 1)
} else {
  cat("Data already in log2 scale, skipped.\n")
}

# Low-expression filter: require raw expression > 8 (log2 > 3) in at least 10% of samples
keep <- rowSums(expr_clean > 8) >= 0.1 * ncol(expr_clean)
expr_filtered <- expr_clean[keep, ]
cat(sprintf("Probes kept after filtering: %d (%.1f%%)\n",
            sum(keep), 100 * sum(keep) / nrow(expr_clean)))

# ------------------------------------------------------------------------------
# 4. Sample correlation (QC check for outliers)
# ------------------------------------------------------------------------------
cat("\n--- 4. Sample correlation ---\n")

cor_mat <- cor(expr_filtered, use = "complete.obs")
mdd_idx <- which(pheno_54563$disease == "MDD")
ctrl_idx <- which(pheno_54563$disease == "Control")

if (length(mdd_idx) > 1) {
  mdd_cor <- mean(cor_mat[mdd_idx, mdd_idx][lower.tri(cor_mat[mdd_idx, mdd_idx])])
  cat(sprintf("MDD within-group correlation: %.4f\n", mdd_cor))
}
if (length(ctrl_idx) > 1) {
  ctrl_cor <- mean(cor_mat[ctrl_idx, ctrl_idx][lower.tri(cor_mat[ctrl_idx, ctrl_idx])])
  cat(sprintf("Control within-group correlation: %.4f\n", ctrl_cor))
}
if (length(mdd_idx) > 0 && length(ctrl_idx) > 0) {
  cat(sprintf("Between-group correlation: %.4f\n",
              mean(cor_mat[mdd_idx, ctrl_idx], na.rm = TRUE)))
}

# ------------------------------------------------------------------------------
# 5. Differential expression analysis (Paired design)
# ------------------------------------------------------------------------------
cat("\n--- 5. Paired Limma DEA ---\n")

# Design: group + paired factor to account for matching
group <- factor(pheno_54563$disease, levels = c("Control", "MDD"))
pair <- factor(pheno_54563$pair)

# Model matrix: group + pair (this absorbs pair-specific variation)
design <- model.matrix(~ group + pair)
# Rename columns for clarity
colnames(design)[2] <- "MDD_vs_Control"

cat(sprintf("Design matrix: %d samples, %d coefficients (including pairs)\n",
            nrow(design), ncol(design)))

# Fit linear model
fit <- lmFit(expr_filtered, design)

# Define contrast (MDD - Control) – note the coefficient name
contrast_matrix <- makeContrasts(
  MDD_vs_Control = `MDD_vs_Control`,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract results
deg <- topTable(fit2, coef = "MDD_vs_Control", number = Inf, adjust.method = "BH")
deg$probe_id <- rownames(deg)
# ---- Gene annotation & dedup (keep highest AvgExpr) ----
library(AnnotationDbi)
library(illuminaHumanv3.db)
probe_means <- rowMeans(expr_filtered, na.rm = TRUE)
symbols <- mapIds(illuminaHumanv3.db, keys = rownames(deg), column = "SYMBOL", 
                  keytype = "PROBEID", multiVals = "first")
deg$GeneSymbol <- symbols
deg <- deg[!is.na(deg$GeneSymbol) & deg$GeneSymbol != "", ]
deg$AvgExpr <- probe_means[rownames(deg)]
deg <- deg[order(deg$AvgExpr, decreasing = TRUE), ]
deg <- deg[!duplicated(deg$GeneSymbol), ]
cat(sprintf("DEA complete. %d probes tested.\n", nrow(deg)))

# ------------------------------------------------------------------------------
# 6. Compute Cohen's d (effect size)
# ------------------------------------------------------------------------------
cat("\n--- 6. Effect size (Cohen's d) ---\n")

ctrl_expr <- expr_filtered[, ctrl_idx, drop = FALSE]
mdd_expr  <- expr_filtered[, mdd_idx, drop = FALSE]

ctrl_mean <- rowMeans(ctrl_expr, na.rm = TRUE)
mdd_mean  <- rowMeans(mdd_expr, na.rm = TRUE)
ctrl_sd   <- rowSds(ctrl_expr, na.rm = TRUE)
mdd_sd    <- rowSds(mdd_expr, na.rm = TRUE)

cohens_d <- (mdd_mean - ctrl_mean) / sqrt((ctrl_sd^2 + mdd_sd^2) / 2)
deg$cohens_d <- cohens_d[rownames(deg)]

cat(sprintf("Mean |Cohen's d|: %.4f\n", mean(abs(cohens_d), na.rm = TRUE)))

# ------------------------------------------------------------------------------
# 7. Extract Top 16 genes and summary statistics
# ------------------------------------------------------------------------------
cat("\n--- 7. Results summary ---\n")

n_sig_adj <- sum(deg$adj.P.Val < 0.05, na.rm = TRUE)
n_sig_fc <- sum(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 0.5, na.rm = TRUE)

cat(sprintf("Significant (adj.P<0.05): %d\n", n_sig_adj))
cat(sprintf("Significant + |logFC|>0.5: %d\n", n_sig_fc))

# Top 16
top16 <- deg[order(deg$P.Value), ][1:16, c("probe_id", "logFC", "cohens_d", "P.Value", "adj.P.Val")]
print(top16)

# ------------------------------------------------------------------------------
# 8. Save outputs
# ------------------------------------------------------------------------------
cat("\n--- 8. Saving results ---\n")

write.csv(top16, "GSE54563_Top16_DEGs.csv", row.names = FALSE)
write.csv(deg[, c("probe_id", "logFC", "cohens_d", "P.Value", "adj.P.Val")],
          "GSE54563_All_DEG_Results.csv", row.names = FALSE)

# Summary table
summary_df <- data.frame(
  Metric = c("Total Samples", "MDD", "Control", "Pairs",
             "Probes Analyzed", "Sig DEG (adj.P<0.05)",
             "Sig DEG (|logFC|>0.5 & adj.P<0.05)", "Mean Cohen's d"),
  Value = c(nrow(pheno_54563),
            sum(pheno_54563$disease == "MDD"),
            sum(pheno_54563$disease == "Control"),
            length(unique(pheno_54563$pair)),
            nrow(expr_filtered),
            n_sig_adj, n_sig_fc,
            sprintf("%.4f", mean(abs(cohens_d), na.rm = TRUE)))
)
write.csv(summary_df, "GSE54563_Analysis_Summary.csv", row.names = FALSE)
# ============================================================================
# Enrichment Analysis: GO BP, KEGG, Reactome
# ============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

# ---- 1. Prepare gene list (P < 0.05) ----
genes <- deg[deg$P.Value < 0.05, "GeneSymbol"]
cat(sprintf("Genes with P<0.05: %d\n", length(genes)))

if (length(genes) >= 10) {
  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  cat(sprintf("Mapped ENTREZ IDs: %d\n", nrow(entrez)))
  
  if (nrow(entrez) > 0) {
    # ---- 2. GO Biological Process ----
    cat("\n--- Running GO BP enrichment ---\n")
    go <- enrichGO(
      gene = entrez$ENTREZID,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    if (!is.null(go) && nrow(go) > 0) {
      cat(sprintf("GO enriched terms: %d\n", nrow(go)))
      write.csv(as.data.frame(go), "GO_BP_Enrichment.csv", row.names = FALSE)
      save_sci_plot(dotplot(go, showCategory = 15), "GO_BP_Dotplot.tiff")
    } else {
      cat("No significant GO terms.\n")
    }
    
    # ---- 3. KEGG Pathway (use local KEGG.db if available) ----
    cat("\n--- Running KEGG enrichment ---\n")
    options(timeout = 600)  # in case of network download
    if (requireNamespace("KEGG.db", quietly = TRUE)) {
      kegg <- enrichKEGG(
        gene = entrez$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.2,
        use_internal_data = TRUE
      )
    } else {
      kegg <- enrichKEGG(
        gene = entrez$ENTREZID,
        organism = "hsa",
        pvalueCutoff = 0.1,
        qvalueCutoff = 0.2
      )
    }
    if (!is.null(kegg) && nrow(kegg) > 0) {
      cat(sprintf("KEGG enriched pathways: %d\n", nrow(kegg)))
      kegg_readable <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
      write.csv(as.data.frame(kegg_readable), "KEGG_Enrichment.csv", row.names = FALSE)
      save_sci_plot(barplot(kegg, showCategory = 15), "KEGG_Barplot.tiff")
    } else {
      cat("No significant KEGG pathways.\n")
    }
    
    # ---- 4. Reactome Pathway ----
    cat("\n--- Running Reactome enrichment ---\n")
    react <- enrichPathway(
      gene = entrez$ENTREZID,
      organism = "human",
      pvalueCutoff = 0.05,
      readable = TRUE
    )
    if (!is.null(react) && nrow(react) > 0) {
      cat(sprintf("Reactome enriched pathways: %d\n", nrow(react)))
      write.csv(as.data.frame(react), "Reactome_Enrichment.csv", row.names = FALSE)
      save_sci_plot(dotplot(react, showCategory = 15), "Reactome_Dotplot.tiff")
    } else {
      cat("No significant Reactome pathways.\n")
    }
    
  } else {
    cat("No ENTREZ IDs to map.\n")
  }
} else {
  cat("Too few genes (<10) for enrichment.\n")
}

cat("\n=== Enrichment analysis completed ===\n")