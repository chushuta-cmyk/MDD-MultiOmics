# ==============================================================================
# GSE80655 → MOFA Input (VST Normalized) - NO GEO DOWNLOAD
# Avoids GEOquery timeout by using local files or skipping GPL
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("    GSE80655 → MOFA Input (VST Normalized)\n")
cat("    NO GEO DOWNLOAD VERSION\n")
cat("====================================================================\n\n")

# Setup
alt_lib <- "/Users/tan/Developer/envs/R/lib/R/library"
if (dir.exists(alt_lib)) {
  .libPaths(c(alt_lib, .libPaths()))
  cat("✓ Added library path\n")
}

library(data.table)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

if (!require("DESeq2", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
}
library(DESeq2)

cat("✓ Packages loaded\n\n")

# ==============================================================================
# Step 1: Load raw counts
# ==============================================================================

cat("【Step 1】Loading GSE80655 raw counts...\n")

file_path <- "/Users/tan/Developer/projects/r_project/Study/Psychology_study/datasets/GSE80655_raw_counts_GRCh38.p13_NCBI.tsv"

if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}

GSE80655_data <- read.delim(
  file = file_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

cat(sprintf("  Raw counts: %d genes × %d samples\n", nrow(GSE80655_data), ncol(GSE80655_data)))

# ==============================================================================
# Step 2: ENTREZID → SYMBOL mapping
# ==============================================================================

cat("\n【Step 2】Mapping ENTREZID to gene symbols...\n")

gene_ids <- rownames(GSE80655_data)

symbols <- AnnotationDbi::mapIds(
  x = org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Remove NA mappings
mapped_indices <- !is.na(symbols)
valid_symbols <- symbols[mapped_indices]
valid_matrix <- GSE80655_data[mapped_indices, ]

cat(sprintf("  Mapped: %d/%d genes\n", length(valid_symbols), length(symbols)))

# Deduplicate: keep highest sum per symbol
gene_sums <- rowSums(valid_matrix)
symbol_info <- data.frame(
  symbol = valid_symbols,
  sum_count = gene_sums,
  index = 1:length(valid_symbols),
  stringsAsFactors = FALSE
)

symbol_info_sorted <- symbol_info[order(symbol_info$sum_count, decreasing = TRUE), ]
unique_info <- symbol_info_sorted[!duplicated(symbol_info_sorted$symbol), ]

final_matrix_unique <- valid_matrix[unique_info$index, ]
rownames(final_matrix_unique) <- unique_info$symbol

cat(sprintf("  After deduplication: %d unique genes\n", nrow(final_matrix_unique)))

# ==============================================================================
# Step 3: Load metadata from LOCAL series matrix (avoid GEO download)
# ==============================================================================

cat("\n【Step 3】Loading metadata from local series matrix...\n")

# Try to find local series matrix file
local_series <- "/Users/tan/Developer/projects/r_project/Study/Psychology_study/datasets/GSE80655_series_matrix.txt"
local_series_gz <- "/Users/tan/Developer/projects/r_project/Study/Psychology_study/datasets/GSE80655_series_matrix.txt.gz"

series_file <- NULL
if (file.exists(local_series)) {
  series_file <- local_series
  cat("  Found local series matrix (uncompressed)\n")
} else if (file.exists(local_series_gz)) {
  series_file <- local_series_gz
  cat("  Found local series matrix (gzipped)\n")
}

metadata_loaded <- FALSE

if (!is.null(series_file)) {
  tryCatch({
    # Read series matrix header to extract metadata
    lines <- readLines(series_file, n = 100)

    # Extract sample titles
    title_line <- lines[grep("^!Sample_title", lines)]
    if (length(title_line) > 0) {
      title_parts <- unlist(strsplit(title_line[1], "\t"))
      sample_titles <- gsub('"', '', title_parts[-1])
      cat(sprintf("  Found %d sample titles\n", length(sample_titles)))
    }

    # Extract characteristics
    char_lines <- lines[grep("^!Sample_characteristics_ch1", lines)]
    cat(sprintf("  Found %d characteristic lines\n", length(char_lines)))

    # Try to extract diagnosis and brain region from titles
    # GSE80655 titles typically contain brain region and diagnosis info
    diagnosis <- rep(NA, length(sample_titles))
    brain_region <- rep(NA, length(sample_titles))

    for (i in seq_along(sample_titles)) {
      title <- sample_titles[i]

      # Extract diagnosis
      if (grepl("Control", title, ignore.case = TRUE)) {
        diagnosis[i] <- "Control"
      } else if (grepl("MDD|Depression|Major", title, ignore.case = TRUE)) {
        diagnosis[i] <- "MDD"
      } else if (grepl("BD|Bipolar", title, ignore.case = TRUE)) {
        diagnosis[i] <- "BD"
      } else if (grepl("SZ|Schizo", title, ignore.case = TRUE)) {
        diagnosis[i] <- "SZ"
      }

      # Extract brain region
      if (grepl("AnCg|Anterior.cingulate|anterior_cingulate", title, ignore.case = TRUE)) {
        brain_region[i] <- "AnCg"
      } else if (grepl("DLPFC|dorsolateral|Dorsolateral", title, ignore.case = TRUE)) {
        brain_region[i] <- "DLPFC"
      } else if (grepl("nAcc|nucleus.accumbens|Nucleus", title, ignore.case = TRUE)) {
        brain_region[i] <- "nAcc"
      }
    }

    # Create metadata dataframe
    sample_metadata <- data.frame(
      GSM_ID = colnames(final_matrix_unique)[1:length(sample_titles)],
      Sample_Title = sample_titles,
      Diagnosis = diagnosis,
      Brain_Region = brain_region,
      stringsAsFactors = FALSE
    )

    # Filter to MDD and Control only
    sample_metadata <- sample_metadata[sample_metadata$Diagnosis %in% c("MDD", "Control"), ]
    sample_metadata <- sample_metadata[!is.na(sample_metadata$Diagnosis), ]

    cat(sprintf("\n  MDD/Control samples: %d\n", nrow(sample_metadata)))
    cat("  Diagnosis distribution:\n")
    print(table(sample_metadata$Diagnosis))
    cat("  Brain region distribution:\n")
    print(table(sample_metadata$Brain_Region))

    metadata_loaded <- TRUE

  }, error = function(e) {
    cat("  ✗ Failed to parse local series matrix:", e$message, "\n")
  })
}

# If local parsing failed, try GEOquery with getGPL=FALSE (faster, no GPL download)
if (!metadata_loaded) {
  cat("\n【Step 3b】Trying GEOquery with getGPL=FALSE...\n")

  if (!require("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
    library(GEOquery)
  }

  tryCatch({
    gse_data <- getGEO("GSE80655", GSEMatrix = TRUE, getGPL = FALSE, AnnotGPL = FALSE)
    gse_eset <- gse_data[[1]]
    phenotype_df <- pData(gse_eset)

    # Extract relevant columns
    phenotype_selected <- phenotype_df[, c("geo_accession", "brain region:ch1", "clinical diagnosis:ch1")]

    target_diagnoses <- c("Control", "Major Depression")
    metadata_selected <- phenotype_selected[
      phenotype_selected$`clinical diagnosis:ch1` %in% target_diagnoses,
    ]

    sample_metadata <- data.frame(
      GSM_ID = metadata_selected$geo_accession,
      Diagnosis = ifelse(metadata_selected$`clinical diagnosis:ch1` == "Control", "Control", "MDD"),
      Brain_Region = gsub(" ", "_", metadata_selected$`brain region:ch1`),
      stringsAsFactors = FALSE
    )

    cat(sprintf("  Loaded %d samples from GEO\n", nrow(sample_metadata)))
    metadata_loaded <- TRUE

  }, error = function(e) {
    cat("  ✗ GEOquery also failed:", e$message, "\n")
  })
}

# If all else fails, create minimal metadata from column names
if (!metadata_loaded) {
  cat("\n【Step 3c】Creating minimal metadata from column names...\n")

  # Use all columns as samples
  sample_ids <- colnames(final_matrix_unique)
  sample_metadata <- data.frame(
    GSM_ID = sample_ids,
    Diagnosis = "Unknown",
    Brain_Region = "Unknown",
    stringsAsFactors = FALSE
  )

  cat("  ⚠ WARNING: Using minimal metadata (no group info)\n")
  cat("  All samples will be treated as one group\n")
}

# ==============================================================================
# Step 4: Subset expression matrix
# ==============================================================================

cat("\n【Step 4】Subsetting expression matrix...\n")

selected_gsm_ids <- sample_metadata$GSM_ID
matrix_cols <- colnames(final_matrix_unique)
matching_gsm_ids <- selected_gsm_ids[selected_gsm_ids %in% matrix_cols]

missing_gsm_ids <- selected_gsm_ids[!selected_gsm_ids %in% matrix_cols]
if (length(missing_gsm_ids) > 0) {
  cat(sprintf("  ⚠ %d samples not found in expression matrix\n", length(missing_gsm_ids)))
}

expr_data_final <- as.matrix(final_matrix_unique[, matching_gsm_ids, drop = FALSE])

# Align metadata
sample_metadata_final <- sample_metadata[
  match(colnames(expr_data_final), sample_metadata$GSM_ID),
]

cat(sprintf("  Final matrix: %d genes × %d samples\n", nrow(expr_data_final), ncol(expr_data_final)))

# ==============================================================================
# Step 5: DESeq2 VST normalization
# ==============================================================================

cat("\n【Step 5】DESeq2 VST normalization...\n")

# Round to integers and convert
expr_data_counts <- round(expr_data_final)
storage.mode(expr_data_counts) <- "integer"

# Remove zero-count genes
expr_data_counts <- expr_data_counts[rowSums(expr_data_counts) > 0, ]

cat(sprintf("  After removing zero-count genes: %d genes\n", nrow(expr_data_counts)))

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = expr_data_counts,
  colData = sample_metadata_final,
  design = ~ 1  # No design, just normalization
)

# VST normalization (blind = TRUE for unsupervised MOFA)
cat("  Running VST...\n")
vsd <- vst(dds, blind = TRUE)
expr_vst <- assay(vsd)

cat(sprintf("  VST matrix: %d genes × %d samples\n", nrow(expr_vst), ncol(expr_vst)))

# ==============================================================================
# Step 6: Save for MOFA
# ==============================================================================

cat("\n【Step 6】Saving to MOFA inputs...\n")

output_dir <- "/Users/tan/Desktop/MDD/next_steps/MOFA/inputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(expr_vst, file.path(output_dir, "GSE80655_expression_cleaned.rds"))
saveRDS(sample_metadata_final, file.path(output_dir, "GSE80655_metadata.rds"))
write.csv(expr_vst, file.path(output_dir, "GSE80655_expression_cleaned.csv"))

cat("  ✓ Saved: GSE80655_expression_cleaned.rds\n")
cat("  ✓ Saved: GSE80655_metadata.rds\n")
cat("  ✓ Saved: GSE80655_expression_cleaned.csv\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n====================================================================\n")
cat("    GSE80655 → MOFA Input COMPLETE\n")
cat("====================================================================\n")
cat(sprintf("Input:  %d genes × %d samples (raw counts)\n", nrow(GSE80655_data), ncol(GSE80655_data)))
cat(sprintf("Output: %d genes × %d samples (VST normalized)\n", nrow(expr_vst), ncol(expr_vst)))
if ("Diagnosis" %in% colnames(sample_metadata_final)) {
  cat("Diagnosis distribution:\n")
  print(table(sample_metadata_final$Diagnosis))
}
if ("Brain_Region" %in% colnames(sample_metadata_final)) {
  cat("Brain region distribution:\n")
  print(table(sample_metadata_final$Brain_Region))
}
cat("Files saved to:", output_dir, "\n")
cat("====================================================================\n")
