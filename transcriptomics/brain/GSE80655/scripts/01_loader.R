# ==============================================================================
# 01_loader.R
# Load raw count matrix, map ENTREZID to Gene Symbol, deduplicate by max total count
# ==============================================================================

library(AnnotationDbi)
library(org.Hs.eg.db)

cat("\n=== GSE80655 Loader Started ===\n")

# ------------------------------------------------------------------------------
# 1. Load raw count matrix
# ------------------------------------------------------------------------------
file_path <- "/Users/tan/Developer/projects/r_project/Study/Psychology_study/datasets/GSE80655_raw_counts_GRCh38.p13_NCBI.tsv"

GSE80655_data <- read.delim(
  file = file_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

cat(sprintf("Raw matrix dimensions: %d rows x %d columns\n", nrow(GSE80655_data), ncol(GSE80655_data)))

# ------------------------------------------------------------------------------
# 2. ID mapping: ENTREZID -> Gene Symbol
# ------------------------------------------------------------------------------
gene_ids <- rownames(GSE80655_data)

symbols <- mapIds(
  x = org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

cat(sprintf("Mapped %d out of %d genes to SYMBOL\n", sum(!is.na(symbols)), length(symbols)))

# ------------------------------------------------------------------------------
# 3. Remove unmapped genes (NA)
# ------------------------------------------------------------------------------
mapped_idx <- !is.na(symbols)
valid_symbols <- symbols[mapped_idx]
valid_matrix <- GSE80655_data[mapped_idx, ]

cat(sprintf("Removed %d unmapped genes. Remaining: %d\n", 
            sum(!mapped_idx), nrow(valid_matrix)))

# ------------------------------------------------------------------------------
# 4. Deduplicate: keep the probe with highest total row sum per gene symbol
# ------------------------------------------------------------------------------
gene_sums <- rowSums(valid_matrix)

symbol_info <- data.frame(
  symbol = valid_symbols,
  sum_count = gene_sums,
  index = 1:length(valid_symbols)
)

symbol_info_sorted <- symbol_info[order(symbol_info$sum_count, decreasing = TRUE), ]
unique_info <- symbol_info_sorted[!duplicated(symbol_info_sorted$symbol), ]

final_matrix_unique <- valid_matrix[unique_info$index, ]
rownames(final_matrix_unique) <- unique_info$symbol

cat(sprintf("Deduplicated: %d unique gene symbols retained (removed %d duplicates)\n",
            nrow(final_matrix_unique), 
            nrow(valid_matrix) - nrow(final_matrix_unique)))

# ------------------------------------------------------------------------------
# 5. Save intermediate objects
# ------------------------------------------------------------------------------
saveRDS(final_matrix_unique, file = "GSE80655_raw_counts_dedup.rds")
cat("Saved deduplicated count matrix to: GSE80655_raw_counts_dedup.rds\n")

# ------------------------------------------------------------------------------
# 6. Checklist
# ------------------------------------------------------------------------------
cat("\n=== CHECKLIST: 01_loader.R ===\n")
cat("✅ Raw count matrix loaded successfully\n")
cat(sprintf("  - Total probes: %d\n", nrow(GSE80655_data)))
cat(sprintf("  - Total samples: %d\n", ncol(GSE80655_data)))
cat("✅ ENTREZID -> SYMBOL mapping completed\n")
cat(sprintf("  - Mapped: %d / %d (%.1f%%)\n", 
            sum(!is.na(symbols)), length(symbols), 
            100 * sum(!is.na(symbols)) / length(symbols)))
cat("✅ Unmapped probes removed\n")
cat("✅ Deduplication by max total count completed\n")
cat(sprintf("  - Unique genes retained: %d\n", nrow(final_matrix_unique)))
cat("✅ Intermediate object saved: GSE80655_raw_counts_dedup.rds\n")
cat("=== End of 01_loader.R ===\n")

