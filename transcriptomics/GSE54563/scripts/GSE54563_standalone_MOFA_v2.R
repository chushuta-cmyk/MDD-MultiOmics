# ==============================================================================
# GSE54563 Standalone Preprocessing → MOFA Input
# Platform: Illumina HumanHT-12 V3.0 (GPL6947)
# Tissue: Anterior Cingulate Cortex (ACC)
# Design: 25 MDD-Control matched pairs
# ==============================================================================
# This script is COMPLETELY STANDALONE - does NOT depend on previous RStudio objects
# It loads the raw series matrix, applies all bugfixes, maps gene symbols, and saves
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("    GSE54563 Standalone Preprocessing → MOFA Input\n")
cat("    Platform: Illumina HumanHT-12 V3.0 (GPL6947)\n")
cat("====================================================================\n\n")

# ==============================================================================
# Step 0: Setup
# ==============================================================================

alt_lib <- "/Users/tan/Developer/envs/R/lib/R/library"
if (dir.exists(alt_lib)) {
  .libPaths(c(alt_lib, .libPaths()))
  cat("✓ Added library path:", alt_lib, "\n")
}

library(data.table)
library(dplyr)
library(limma)
library(matrixStats)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("illuminaHumanv3.db", quietly = TRUE)) {
  BiocManager::install("illuminaHumanv3.db", ask = FALSE, update = FALSE)
}
library(illuminaHumanv3.db)
library(AnnotationDbi)

cat("✓ All packages loaded\n\n")

# ==============================================================================
# Step 1: Read phenotype data (lines 29-30)
# ==============================================================================

cat("【Step 1】Reading phenotype data...\n")

setwd("/Users/tan/Developer/projects/r_project/Study/GSE54563/")

pheno_raw <- fread(
  "GSE54563_series_matrix.txt",
  skip = 28,
  sep = "\t",
  header = FALSE,
  nrows = 2,
  data.table = FALSE,
  fill = TRUE
)

sample_title_raw <- as.character(pheno_raw[1, -1])
sample_geo_id_raw <- as.character(pheno_raw[2, -1])
sample_title <- gsub('"', '', sample_title_raw)
sample_geo_id <- gsub('"', '', sample_geo_id_raw)

# Extract groups
disease_group <- ifelse(
  grepl("^MDD matched", sample_title, ignore.case = TRUE),
  "MDD",
  "Control"
)

# Extract pair numbers
pair_number <- as.numeric(gsub(".*#([0-9]+)$", "\1", sample_title))

cat(sprintf("  Total samples: %d\n", length(sample_title)))
cat(sprintf("  MDD: %d, Control: %d\n", sum(disease_group == "MDD"), sum(disease_group == "Control")))
cat(sprintf("  Pairs: %d\n", length(unique(pair_number))))

# Create phenotype dataframe
pheno_54563 <- data.frame(
  sample_id = sample_geo_id,
  sample_title = sample_title,
  disease_group = disease_group,
  pair_id = pair_number,
  stringsAsFactors = FALSE
)

# ==============================================================================
# Step 2: Read expression matrix (from line 67)
# ==============================================================================

cat("\n【Step 2】Reading expression matrix...\n")

expr_raw <- fread(
  "GSE54563_series_matrix.txt",
  skip = 66,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  fill = TRUE
)

probe_ids <- as.character(expr_raw[, 1])
expr_matrix <- as.matrix(expr_raw[, -1])
expr_matrix <- apply(expr_matrix, 2, as.numeric)
rownames(expr_matrix) <- probe_ids

cat(sprintf("  Raw matrix: %d probes × %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))

# Align columns with phenotype
expr_sample_ids <- colnames(expr_raw)[-1]
if (!all(expr_sample_ids == sample_geo_id)) {
  expr_matrix <- expr_matrix[, match(sample_geo_id, expr_sample_ids)]
  cat("  ✓ Reordered samples to match phenotype\n")
}

# ==============================================================================
# Step 3: BUG FIX - Handle NA values
# ==============================================================================

cat("\n【Step 3】Handling NA values...\n")

na_rows <- rowSums(is.na(expr_matrix)) > 0
cat(sprintf("  Rows with NA: %d\n", sum(na_rows)))

expr_clean <- expr_matrix[!na_rows, ]
cat(sprintf("  After NA removal: %d probes × %d samples\n", nrow(expr_clean), ncol(expr_clean)))

# ==============================================================================
# Step 4: BUG FIX - Correct probe filtering
# ==============================================================================

cat("\n【Step 4】Filtering low-expression probes...\n")

# Log2 transform first (add 1 to avoid log(0))
expr_log2 <- log2(expr_clean + 1)

# Filter: keep probes with expression > 3 (log2(8)) in at least 10% of samples
keep_probes <- rowSums(expr_log2 > 3) >= ncol(expr_log2) * 0.1

cat(sprintf("  Original probes: %d\n", nrow(expr_clean)))
cat(sprintf("  Retained probes: %d\n", sum(keep_probes)))
cat(sprintf("  Removed probes: %d (%.1f%%)\n", sum(!keep_probes), 100*sum(!keep_probes)/length(keep_probes)))

expr_filtered <- expr_clean[keep_probes, ]
cat(sprintf("  Filtered matrix: %d probes × %d samples\n", nrow(expr_filtered), ncol(expr_filtered)))

# ==============================================================================
# Step 5: Map ILMN probe IDs to gene symbols
# ==============================================================================

cat("\n【Step 5】Mapping probe IDs to gene symbols...\n")

probe_ids <- rownames(expr_filtered)
cat(sprintf("  First 10 probes: %s\n", paste(head(probe_ids, 10), collapse = ", ")))

# Check overlap with illuminaHumanv3.db
db_probes <- keys(illuminaHumanv3.db, keytype = "PROBEID")
overlap <- sum(probe_ids %in% db_probes)
cat(sprintf("  %d/%d probes found in illuminaHumanv3.db\n", overlap, length(probe_ids)))

# If no overlap, try without ILMN_ prefix
if (overlap == 0 && grepl("^ILMN_", probe_ids[1])) {
  cat("  Trying without ILMN_ prefix...\n")
  probe_ids_noprefix <- gsub("^ILMN_", "", probe_ids)
  overlap2 <- sum(probe_ids_noprefix %in% db_probes)
  cat(sprintf("  Without prefix: %d matches\n", overlap2))
  if (overlap2 > 0) {
    probe_ids <- probe_ids_noprefix
  }
}

# Map to gene symbols
mapped <- AnnotationDbi::select(illuminaHumanv3.db,
                                  keys = probe_ids[probe_ids %in% db_probes],
                                  columns = c("SYMBOL", "ENTREZID"),
                                  keytype = "PROBEID")

cat(sprintf("  Retrieved %d mappings\n", nrow(mapped)))

# Remove duplicate probes (keep first)
mapped_dedup <- mapped[!duplicated(mapped$PROBEID), ]
gene_symbols <- mapped_dedup$SYMBOL
names(gene_symbols) <- mapped_dedup$PROBEID

cat(sprintf("  Unique gene symbols: %d\n", sum(!is.na(gene_symbols))))

# ==============================================================================
# Step 6: Create gene-level expression matrix
# ==============================================================================

cat("\n【Step 6】Creating gene-level matrix...\n")

# Build dataframe with gene symbols
expr_df <- data.frame(
  ProbeID = rownames(expr_filtered),
  GeneSymbol = gene_symbols[rownames(expr_filtered)],
  expr_filtered,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Remove unmapped probes
expr_mapped <- expr_df[!is.na(expr_df$GeneSymbol), ]
cat(sprintf("  Probes with symbols: %d\n", nrow(expr_mapped)))

# Deduplicate: keep highest mean expression per gene
expr_mapped$mean_expr <- rowMeans(expr_mapped[, -(1:2), drop = FALSE])
expr_mapped <- expr_mapped[order(expr_mapped$GeneSymbol, -expr_mapped$mean_expr), ]
expr_unique <- expr_mapped[!duplicated(expr_mapped$GeneSymbol), ]
expr_unique$mean_expr <- NULL

gene_vec <- expr_unique$GeneSymbol
expr_final <- as.matrix(expr_unique[, -(1:2), drop = FALSE])
rownames(expr_final) <- gene_vec

cat(sprintf("  Final: %d genes × %d samples\n", nrow(expr_final), ncol(expr_final)))

# ==============================================================================
# Step 7: Save to MOFA inputs
# ==============================================================================

cat("\n【Step 7】Saving to MOFA inputs...\n")

output_dir <- "/Users/tan/Desktop/MDD/next_steps/MOFA/inputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(expr_final, file.path(output_dir, "GSE54563_expression_cleaned.rds"))
saveRDS(pheno_54563, file.path(output_dir, "GSE54563_metadata.rds"))
write.csv(expr_final, file.path(output_dir, "GSE54563_expression_cleaned.csv"))

cat("  ✓ Saved: GSE54563_expression_cleaned.rds\n")
cat("  ✓ Saved: GSE54563_metadata.rds\n")
cat("  ✓ Saved: GSE54563_expression_cleaned.csv\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n====================================================================\n")
cat("    GSE54563 → MOFA Input COMPLETE\n")
cat("====================================================================\n")
cat(sprintf("Input:  %d probes × %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))
cat(sprintf("Output: %d genes × %d samples\n", nrow(expr_final), ncol(expr_final)))
cat(sprintf("Groups: %d MDD, %d Control (25 pairs)\n", 
            sum(disease_group == "MDD"), sum(disease_group == "Control")))
cat("Files saved to:", output_dir, "\n")
cat("====================================================================\n")
