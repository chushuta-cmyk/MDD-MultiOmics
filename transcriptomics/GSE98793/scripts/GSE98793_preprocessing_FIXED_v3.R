# ==============================================================================
# GSE98793 Preprocessing - FIXED VERSION v3
# Fixes all issues from previous conversation
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("    GSE98793 Preprocessing - FIXED v3\n")
cat("    Platform: Affymetrix HG-U133 Plus 2.0 (GPL570)\n")
cat("    Issue Fixes: library path, annotation package, probe mapping, select() conflict\n")
cat("====================================================================\n\n")

# ==============================================================================
# FIX 1: Set correct library path BEFORE loading any packages
# ==============================================================================

# Your packages are installed here:
alt_lib <- "/Users/tan/Developer/envs/R/lib/R/library"
if (dir.exists(alt_lib)) {
  .libPaths(c(alt_lib, .libPaths()))
  cat("✓ Added library path:", alt_lib, "\n")
} else {
  cat("⚠ Warning: Alternative library path not found, using default\n")
}

cat("Current .libPaths():\n")
for (p in .libPaths()) {
  cat("  ", p, "\n")
}
cat("\n")

# ==============================================================================
# FIX 2: Load required packages with correct versions
# ==============================================================================

cat("【Step 0】Loading packages...\n")

required_packages <- c(
  "GEOquery",        # GEO data access
  "limma",           # Differential expression
  "hgu133plus2.db",  # CORRECT annotation package for Affymetrix HG-U133 Plus 2.0
  "org.Hs.eg.db",    # Human gene annotation
  "AnnotationDbi",   # Annotation infrastructure
  "ggplot2",         # Plotting
  "dplyr",           # Data manipulation
  "pheatmap",        # Heatmaps
  "matrixStats"      # Matrix operations
)

missing_pkgs <- c()
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, pkg)
    cat("  ✗", pkg, "- NOT FOUND\n")
  } else {
    cat("  ✓", pkg, "- loaded\n")
  }
}

if (length(missing_pkgs) > 0) {
  cat("\n⚠ Missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  cat("Attempting to install missing packages...\n")

  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  for (pkg in missing_pkgs) {
    cat("  Installing", pkg, "...\n")
    tryCatch({
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
      library(pkg, character.only = TRUE)
      cat("    ✓ Success\n")
    }, error = function(e) {
      cat("    ✗ Failed:", e$message, "\n")
    })
  }
}

cat("\n")

# ==============================================================================
# FIX 3: Load GSE98793 data with proper error handling
# ==============================================================================

cat("【Step 1】Loading GSE98793 data...\n")

# Try local file first, then GEO download
local_file <- "/Users/tan/Developer/projects/r_project/Study/GSE98793/GSE98793_analysis/GSE98793_series_matrix.txt"

gse <- NULL

if (file.exists(local_file)) {
  cat("  Found local file:", local_file, "\n")
  tryCatch({
    gse_list <- getGEO(filename = local_file, getGPL = TRUE)
    if (is(gse_list, "list")) {
      gse <- gse_list[[1]]
    } else {
      gse <- gse_list
    }
    cat("  ✓ Loaded from local file\n")
  }, error = function(e) {
    cat("  ✗ Local file failed:", e$message, "\n")
  })
}

if (is.null(gse)) {
  cat("  Attempting GEO download...\n")
  tryCatch({
    gse_list <- getGEO("GSE98793", GSEMatrix = TRUE, getGPL = TRUE)
    gse <- gse_list[[1]]
    cat("  ✓ Downloaded from GEO\n")
  }, error = function(e) {
    cat("  ✗ GEO download failed:", e$message, "\n")
    stop("Cannot load GSE98793 data. Please check file path or internet connection.")
  })
}

# Extract components
expr_matrix <- exprs(gse)
pheno_data <- pData(gse)
feat_data <- fData(gse)

cat(sprintf("  Expression matrix: %d probes × %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))
cat(sprintf("  Phenotype columns: %s\n", paste(colnames(pheno_data)[1:min(5, ncol(pheno_data))], collapse = ", ")))
cat(sprintf("  Feature columns: %s\n", paste(colnames(feat_data)[1:min(5, ncol(feat_data))], collapse = ", ")))

# ==============================================================================
# FIX 4: Check probe IDs and annotation
# ==============================================================================

cat("\n【Step 2】Checking probe IDs...\n")

probe_ids <- rownames(expr_matrix)
cat("  First 10 probe IDs:", paste(head(probe_ids, 10), collapse = ", "), "\n")

# Check if IDs match hgu133plus2.db format (should be like 1007_s_at)
if (grepl("^ILMN_", probe_ids[1])) {
  cat("  ⚠ WARNING: IDs look like Illumina format (ILMN_*)\n")
  cat("    But GSE98793 should be Affymetrix HG-U133 Plus 2.0\n")
  cat("    Please verify the platform!\n")
} else if (grepl("_at$", probe_ids[1])) {
  cat("  ✓ IDs match Affymetrix format (*_at)\n")
} else {
  cat("  ? Unknown ID format. First ID:", probe_ids[1], "\n")
}

# ==============================================================================
# FIX 5: Map probe IDs to gene symbols using CORRECT annotation package
# CRITICAL FIX: Use AnnotationDbi::select to avoid dplyr conflict
# ==============================================================================

cat("\n【Step 3】Mapping probe IDs to gene symbols...\n")

# Use hgu133plus2.db (Affymetrix HG-U133 Plus 2.0)
# NOT illuminaHumanv3.db (which was wrong!)

# Check if hgu133plus2.db is available
if (!require("hgu133plus2.db", quietly = TRUE)) {
  cat("  Installing hgu133plus2.db...\n")
  BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)
  library(hgu133plus2.db)
}

# Get gene symbols - CRITICAL: Use AnnotationDbi::select explicitly
# to avoid conflict with dplyr::select
probe_ids_clean <- probe_ids[probe_ids %in% keys(hgu133plus2.db, keytype = "PROBEID")]
cat(sprintf("  %d/%d probes found in hgu133plus2.db\n", length(probe_ids_clean), length(probe_ids)))

if (length(probe_ids_clean) == 0) {
  cat("  ⚠ No probes mapped! Checking alternative approaches...\n")

  # Try using feature_data from GEO
  if ("Gene Symbol" %in% colnames(feat_data) || "GENE_SYMBOL" %in% colnames(feat_data)) {
    cat("  ✓ Using gene symbols from GEO featureData\n")
    gene_col <- ifelse("Gene Symbol" %in% colnames(feat_data), "Gene Symbol", "GENE_SYMBOL")
    gene_symbols <- feat_data[, gene_col]
    names(gene_symbols) <- rownames(feat_data)
  } else {
    cat("  ✗ No annotation available. Please check data source.\n")
    stop("Cannot map probes to genes.")
  }
} else {
  # Use hgu133plus2.db with AnnotationDbi::select (NOT dplyr::select!)
  cat("  Using AnnotationDbi::select for mapping...\n")
  gene_map <- AnnotationDbi::select(hgu133plus2.db, 
                     keys = probe_ids_clean,
                     columns = c("SYMBOL", "ENTREZID"),
                     keytype = "PROBEID")

  # Handle duplicates: keep first match
  gene_map <- gene_map[!duplicated(gene_map$PROBEID), ]

  # Create annotation vector
  gene_symbols <- gene_map$SYMBOL
  names(gene_symbols) <- gene_map$PROBEID

  cat(sprintf("  ✓ Mapped %d probes to gene symbols\n", sum(!is.na(gene_symbols))))
}

# ==============================================================================
# FIX 6: Create clean expression matrix with gene symbols
# ==============================================================================

cat("\n【Step 4】Creating gene-level expression matrix...\n")

# Add gene symbols to expression matrix
expr_with_genes <- data.frame(
  ProbeID = rownames(expr_matrix),
  GeneSymbol = gene_symbols[rownames(expr_matrix)],
  expr_matrix,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Remove unmapped probes
expr_mapped <- expr_with_genes[!is.na(expr_with_genes$GeneSymbol), ]
cat(sprintf("  Probes with gene symbols: %d\n", nrow(expr_mapped)))

# Handle duplicate gene symbols: keep probe with highest mean expression
expr_mapped$mean_expr <- rowMeans(expr_mapped[, -(1:2), drop = FALSE])
expr_mapped <- expr_mapped[order(expr_mapped$GeneSymbol, -expr_mapped$mean_expr), ]
expr_unique <- expr_mapped[!duplicated(expr_mapped$GeneSymbol), ]
expr_unique$mean_expr <- NULL

cat(sprintf("  Unique genes after deduplication: %d\n", nrow(expr_unique)))

# Final clean matrix
gene_symbols_vec <- expr_unique$GeneSymbol
expr_clean <- as.matrix(expr_unique[, -(1:2), drop = FALSE])
rownames(expr_clean) <- gene_symbols_vec

cat(sprintf("  Final matrix: %d genes × %d samples\n", nrow(expr_clean), ncol(expr_clean)))

# ==============================================================================
# FIX 7: Extract sample groups correctly
# ==============================================================================

cat("\n【Step 5】Extracting sample groups...\n")

# Check title column for group info
cat("  Sample titles (first 10):\n")
print(head(pheno_data$title, 10))

# Extract MDD vs Control from titles
sample_groups <- ifelse(
  grepl("MDD|depression|depressed|patient", pheno_data$title, ignore.case = TRUE),
  "MDD",
  "Control"
)

# Add to phenotype data
pheno_data$disease_group <- sample_groups

cat("\n  Group counts:\n")
print(table(pheno_data$disease_group))

# ==============================================================================
# FIX 8: Save preprocessed data
# ==============================================================================

cat("\n【Step 6】Saving preprocessed data...\n")

output_dir <- "/Users/tan/Desktop/MDD/next_steps/MOFA/inputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save expression matrix
saveRDS(expr_clean, file.path(output_dir, "GSE98793_expression_cleaned.rds"))
cat("  ✓ Saved: GSE98793_expression_cleaned.rds\n")

# Save metadata
saveRDS(pheno_data, file.path(output_dir, "GSE98793_metadata.rds"))
cat("  ✓ Saved: GSE98793_metadata.rds\n")

# Save as CSV for inspection
write.csv(expr_clean, file.path(output_dir, "GSE98793_expression_cleaned.csv"))
cat("  ✓ Saved: GSE98793_expression_cleaned.csv\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n====================================================================\n")
cat("    PREPROCESSING COMPLETE\n")
cat("====================================================================\n")
cat(sprintf("Input:  %d probes × %d samples\n", nrow(expr_matrix), ncol(expr_matrix)))
cat(sprintf("Output: %d genes × %d samples\n", nrow(expr_clean), ncol(expr_clean)))
cat(sprintf("Groups: %d MDD, %d Control\n", 
            sum(pheno_data$disease_group == "MDD"),
            sum(pheno_data$disease_group == "Control")))
cat("\nFiles saved to:", output_dir, "\n")
cat("====================================================================\n")
