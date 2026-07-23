#!/usr/bin/env Rscript
# =============================================================================
# Brain MOFA+ Pipeline - Step by Step (Following Official Bioconductor Tutorial)
# https://bioconductor.org/packages/release/bioc/vignettes/MOFA2/inst/doc/getting_started_R.html
# 
# Datasets: GSE80655 (RNA-seq) + GSE98793 (Microarray) + GSE88890 (Methylation)
# Strategy: Single group, 2 views (RNA, Methylation), NA for missing data
# Format: Long data.frame (recommended for complex data with missing views)
# =============================================================================

# ============================================================================
# STEP 0: SETUP - Configure Python explicitly
# ============================================================================

.libPaths(c("/Users/tan/Developer/envs/R/lib/R/library", .libPaths()))
cat("Step 0: Loading packages...\n")

suppressPackageStartupMessages({
  library(data.table)
  library(MOFA2)
  library(dplyr)
  library(reticulate)
})

cat("  ✓ MOFA2 version:", as.character(packageVersion("MOFA2")), "\n")

# EXPLICITLY set Python path to where mofapy2 is installed
# Based on pip output: /Library/Frameworks/Python.framework/Versions/3.13/bin/python3
python_path <- "/Library/Frameworks/Python.framework/Versions/3.13/bin/python3"

if (!file.exists(python_path)) {
  # Try alternative paths
  alt_paths <- c(
    "/Library/Frameworks/Python.framework/Versions/3.13/bin/python",
    "/usr/local/bin/python3.13",
    "/opt/homebrew/bin/python3.13",
    Sys.which("python3"),
    Sys.which("python")
  )
  for (p in alt_paths) {
    if (file.exists(p)) {
      python_path <- p
      break
    }
  }
}

cat("  Using Python:", python_path, "\n")
use_python(python_path, required = TRUE)

# Verify mofapy2 is available
tryCatch({
  py_run_string("import mofapy2")
  cat("  ✓ mofapy2 is available in Python\n\n")
}, error = function(e) {
  cat("\n⚠️  ERROR: mofapy2 not found in", python_path, "\n")
  cat("  Python error:", e$message, "\n")
  cat("  Please run in terminal:\n")
  cat("    ", python_path, "-m pip install mofapy2\n\n")
  stop("mofapy2 not available")
})

# Configuration
INPUT_DIR <- "/Users/tan/Desktop/MDD/next_steps/MOFA/inputs"
OUTPUT_DIR <- "/Users/tan/Desktop/MDD/next_steps/MOFA/brain_mofa"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# STEP 1: LOAD PREPROCESSED DATA
# ============================================================================

cat("Step 1: Loading preprocessed datasets...\n\n")

# Helper to load expression matrix and metadata
load_dataset <- function(name, expr_file, meta_file, cohort) {
  expr_path <- file.path(INPUT_DIR, expr_file)
  meta_path <- file.path(INPUT_DIR, meta_file)

  expr <- readRDS(expr_path)
  meta <- readRDS(meta_path)

  cat("  ✓", name, ":", nrow(expr), "features x", ncol(expr), "samples\n")
  return(list(expr = as.matrix(expr), meta = meta, cohort = cohort))
}

gse80655 <- load_dataset("GSE80655 (Brain RNA-seq)", 
                         "GSE80655_expression_cleaned.rds", 
                         "GSE80655_metadata.rds", 
                         "GSE80655")

gse98793 <- load_dataset("GSE98793 (Microarray)", 
                         "GSE98793_expression_cleaned.rds", 
                         "GSE98793_metadata.rds", 
                         "GSE98793")

gse88890 <- load_dataset("GSE88890 (Methylation)", 
                         "GSE88890_methylation_gene_level.rds", 
                         "GSE88890_metadata.rds", 
                         "GSE88890")

# ============================================================================
# STEP 2: PREPARE METADATA (with cohort-prefixed sample names)
# ============================================================================

cat("\nStep 2: Preparing metadata...\n")

# Standardize metadata - extract sample ID and condition
std_meta <- function(meta, cohort) {
  # Find sample ID column
  sample_col <- NULL
  for (col in c("sample", "GSM_ID", "geo_accession", "Sample_ID", "Sample_Name", "title")) {
    if (col %in% colnames(meta)) {
      sample_col <- col
      break
    }
  }
  if (is.null(sample_col)) sample_col <- colnames(meta)[1]

  # Find condition column
  cond_col <- NULL
  for (col in c("condition", "Diagnosis", "disease", "status", "Phenotype", 
                "disease_group", "subject group", "characteristics_ch1")) {
    if (col %in% colnames(meta)) {
      cond_col <- col
      break
    }
  }

  df <- data.frame(
    sample = paste0(cohort, "_", as.character(meta[[sample_col]])),
    cohort = cohort,
    stringsAsFactors = FALSE
  )

  if (!is.null(cond_col)) {
    cond <- toupper(as.character(meta[[cond_col]]))
    cond <- gsub("MAJOR DEPRESSIVE DISORDER|DEPRESSION", "MDD", cond)
    cond <- ifelse(grepl("MDD|CASE|PATIENT", cond), "MDD", "Control")
    df$condition <- cond
  } else {
    df$condition <- "unknown"
  }

  return(df)
}

meta_80655 <- std_meta(gse80655$meta, "GSE80655")
meta_98793 <- std_meta(gse98793$meta, "GSE98793")
meta_88890 <- std_meta(gse88890$meta, "GSE88890")

all_meta <- rbind(meta_80655, meta_98793, meta_88890)

cat("  Sample counts by cohort and condition:\n")
print(table(all_meta$cohort, all_meta$condition))
cat("  Total samples:", nrow(all_meta), "\n")

# ============================================================================
# STEP 3: BUILD LONG DATA.FRAME (Official Recommended Format)
# ============================================================================

cat("\nStep 3: Building long data.frame for MOFA...\n")
cat("  (This is the recommended format for complex data with missing views)\n")

# Function to convert matrix to long format
matrix_to_long <- function(mat, view_name, cohort_prefix) {
  # Rename columns with cohort prefix
  colnames(mat) <- paste0(cohort_prefix, "_", colnames(mat))

  # Convert to data.frame
  df <- as.data.frame(mat)
  df$feature <- rownames(df)

  # Pivot to long format
  long_df <- df %>%
    tidyr::pivot_longer(
      cols = -feature,
      names_to = "sample",
      values_to = "value"
    ) %>%
    mutate(view = view_name) %>%
    dplyr::select(sample, feature, view, value)

  return(long_df)
}

# Convert each dataset to long format
# GSE80655: RNA view
long_80655 <- matrix_to_long(gse80655$expr, "RNA", "GSE80655")
long_80655$feature <- paste0("GEX_", long_80655$feature)

# GSE98793: RNA view  
long_98793 <- matrix_to_long(gse98793$expr, "RNA", "GSE98793")
long_98793$feature <- paste0("GEX_", long_98793$feature)

# GSE88890: Methylation view
long_88890 <- matrix_to_long(gse88890$expr, "Methylation", "GSE88890")
long_88890$feature <- paste0("MET_", long_88890$feature)

# Combine all
dt <- rbind(long_80655, long_98793, long_88890)
dt <- as.data.table(dt)

cat("  ✓ Long data.frame created\n")
cat("    Rows:", nrow(dt), "\n")
cat("    Samples:", length(unique(dt$sample)), "\n")
cat("    Views:", paste(unique(dt$view), collapse = ", "), "\n")
cat("    Features (RNA):", length(grep("^GEX_", unique(dt$feature))), "\n")
cat("    Features (Meth):", length(grep("^MET_", unique(dt$feature))), "\n")

# Check data completeness
cat("\n  Data completeness check:\n")
completeness <- dt %>%
  group_by(view) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_features = n_distinct(feature),
    n_values = n(),
    .groups = "drop"
  )
print(completeness)

# ============================================================================
# STEP 4: CREATE MOFA OBJECT FROM DATA.FRAME
# ============================================================================

cat("\nStep 4: Creating MOFA object from data.frame...\n")

MOFAobject <- create_mofa(dt)

cat("  ✓ MOFA object created\n")
print(MOFAobject)

# ============================================================================
# STEP 5: ADD SAMPLE METADATA (match exactly with MOFA object samples)
# ============================================================================

cat("\nStep 5: Adding sample metadata...\n")

# Get sample names from MOFA object - handle list return
mofa_samples_raw <- samples_names(MOFAobject)

# Convert to vector if it's a list
if (is.list(mofa_samples_raw)) {
  if (length(mofa_samples_raw) == 1) {
    mofa_samples <- mofa_samples_raw[[1]]
  } else {
    mofa_samples <- unlist(mofa_samples_raw)
  }
} else {
  mofa_samples <- mofa_samples_raw
}

cat("  MOFA object samples:", length(mofa_samples), "\n")
cat("  First 5:", paste(head(mofa_samples, 5), collapse = ", "), "\n")

# Check which metadata samples match
meta_samples <- all_meta$sample
cat("  Metadata samples:", length(meta_samples), "\n")
cat("  First 5:", paste(head(meta_samples, 5), collapse = ", "), "\n")

# Find matching samples
matched <- intersect(mofa_samples, meta_samples)
cat("  Matched samples:", length(matched), "/", length(mofa_samples), "\n")

# Find non-matching samples
missing_in_meta <- setdiff(mofa_samples, meta_samples)
missing_in_mofa <- setdiff(meta_samples, mofa_samples)

if (length(missing_in_meta) > 0) {
  cat("  ⚠️  Samples in MOFA but not in metadata (first 5):", 
      paste(head(missing_in_meta, 5), collapse = ", "), "\n")
}
if (length(missing_in_mofa) > 0) {
  cat("  ⚠️  Samples in metadata but not in MOFA (first 5):", 
      paste(head(missing_in_mofa, 5), collapse = ", "), "\n")
}

# Reorder metadata to match MOFA sample order
# Only keep samples that exist in both
all_meta_matched <- all_meta[all_meta$sample %in% mofa_samples, ]

# Create a complete metadata frame with all MOFA samples (fill unknown for missing)
meta_for_mofa <- data.frame(
  sample = mofa_samples,
  stringsAsFactors = FALSE
)

# Merge with our metadata
meta_for_mofa <- merge(meta_for_mofa, all_meta_matched, by = "sample", all.x = TRUE)

# Fill NAs for missing metadata
meta_for_mofa$cohort[is.na(meta_for_mofa$cohort)] <- "unknown"
meta_for_mofa$condition[is.na(meta_for_mofa$condition)] <- "unknown"

# Set rownames to sample names
rownames(meta_for_mofa) <- meta_for_mofa$sample

# Add to MOFA object
samples_metadata(MOFAobject) <- meta_for_mofa

cat("  ✓ Metadata added (", nrow(meta_for_mofa), "samples )\n")

# ============================================================================
# STEP 6: PLOT DATA OVERVIEW (Before Training)
# ============================================================================

cat("\nStep 6: Plotting data overview...\n")

pdf(file.path(OUTPUT_DIR, "00_data_overview_before_training.pdf"), width = 10, height = 6)
plot_data_overview(MOFAobject)
dev.off()

cat("  ✓ Saved to:", file.path(OUTPUT_DIR, "00_data_overview_before_training.pdf"), "\n")

# ============================================================================
# STEP 7: DEFINE DATA OPTIONS
# ============================================================================

cat("\nStep 7: Defining data options...\n")

data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE  # Scale each view to unit variance
data_opts$scale_groups <- FALSE

cat("  scale_views:", data_opts$scale_views, "\n")
cat("  scale_groups:", data_opts$scale_groups, "\n")

# ============================================================================
# STEP 8: DEFINE MODEL OPTIONS
# ============================================================================

cat("\nStep 8: Defining model options...\n")

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
model_opts$likelihoods <- c(RNA = "gaussian", Methylation = "gaussian")

cat("  num_factors:", model_opts$num_factors, "\n")
cat("  likelihoods:", paste(model_opts$likelihoods, collapse = ", "), "\n")
cat("  spikeslab_factors:", model_opts$spikeslab_factors, "\n")
cat("  spikeslab_weights:", model_opts$spikeslab_weights, "\n")
cat("  ard_factors:", model_opts$ard_factors, "\n")
cat("  ard_weights:", model_opts$ard_weights, "\n")

# ============================================================================
# STEP 9: DEFINE TRAINING OPTIONS
# ============================================================================

cat("\nStep 9: Defining training options...\n")

train_opts <- get_default_training_options(MOFAobject)
train_opts$maxiter <- 1000
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42
train_opts$verbose <- TRUE
train_opts$drop_factor_threshold <- -1  # Don't auto-drop factors

cat("  maxiter:", train_opts$maxiter, "\n")
cat("  convergence_mode:", train_opts$convergence_mode, "\n")
cat("  seed:", train_opts$seed, "\n")

# ============================================================================
# STEP 10: PREPARE MOFA OBJECT
# ============================================================================

cat("\nStep 10: Preparing MOFA object...\n")

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

cat("  ✓ MOFA object prepared for training\n")

# ============================================================================
# STEP 11: TRAIN MOFA MODEL (using existing Python, NOT basilisk)
# ============================================================================

cat("\nStep 11: Training MOFA model...\n")
cat("  (This may take 10-30 minutes. Watch for ELBO convergence.)\n\n")

outfile <- file.path(OUTPUT_DIR, "brain_mofa_model.hdf5")

# Use use_basilisk = FALSE to avoid downloading Miniforge
# Instead use the Python environment we configured in Step 0
MOFAobject.trained <- run_mofa(MOFAobject, outfile = outfile, use_basilisk = FALSE)

cat("\n  ✓ Training complete!\n")
cat("  Model saved to:", outfile, "\n")

# ============================================================================
# STEP 12: VARIANCE EXPLAINED
# ============================================================================

cat("\nStep 12: Calculating variance explained...\n")

r2 <- get_variance_explained(MOFAobject.trained)

cat("\n  Variance explained per factor:\n")
for (i in 1:ncol(r2$R2Total)) {
  cat("    Factor", i, ":", round(r2$R2Total[1, i] * 100, 2), "% total\n")
}

# Save variance plot
pdf(file.path(OUTPUT_DIR, "01_variance_explained.pdf"), width = 10, height = 6)
plot_variance_explained(MOFAobject.trained, x = "view", y = "factor")
dev.off()

pdf(file.path(OUTPUT_DIR, "02_variance_heatmap.pdf"), width = 8, height = 6)
plot_variance_explained(MOFAobject.trained, x = "factor", y = "view", plot_total = TRUE)
dev.off()

# ============================================================================
# STEP 13: FACTOR PLOTS
# ============================================================================

cat("\nStep 13: Generating factor plots...\n")

pdf(file.path(OUTPUT_DIR, "03_factors_by_cohort.pdf"), width = 14, height = 10)
plot_factor(MOFAobject.trained, factors = 1:6, color_by = "cohort", 
            add_violin = TRUE, dodge = TRUE)
dev.off()

pdf(file.path(OUTPUT_DIR, "04_factors_by_condition.pdf"), width = 14, height = 10)
plot_factor(MOFAobject.trained, factors = 1:6, color_by = "condition", 
            add_violin = TRUE, dodge = TRUE)
dev.off()

pdf(file.path(OUTPUT_DIR, "05_factor_correlation.pdf"), width = 8, height = 7)
plot_factor_correlation(MOFAobject.trained)
dev.off()

# ============================================================================
# STEP 14: FEATURE WEIGHTS
# ============================================================================

cat("\nStep 14: Extracting top feature weights...\n")

for (view in views_names(MOFAobject.trained)) {
  for (factor in 1:min(5, model_opts$num_factors)) {
    w <- get_weights(MOFAobject.trained, view = view, factor = factor)

    top_pos <- head(sort(w, decreasing = TRUE), 30)
    top_neg <- head(sort(w, decreasing = FALSE), 30)

    df <- data.frame(
      Feature = c(names(top_pos), names(top_neg)),
      Weight = c(as.numeric(top_pos), as.numeric(top_neg)),
      Direction = c(rep("Positive", 30), rep("Negative", 30))
    )

    write.csv(df, file.path(OUTPUT_DIR, 
      sprintf("weights_%s_Factor%02d.csv", view, factor)), row.names = FALSE)
  }
}

# ============================================================================
# STEP 15: SAVE RESULTS
# ============================================================================

cat("\nStep 15: Saving results...\n")

# Factor values
factors_df <- get_factors(MOFAobject.trained, as.data.frame = TRUE)
factors_meta <- merge(factors_df, meta_for_mofa, by = "sample")
write.csv(factors_meta, file.path(OUTPUT_DIR, "factor_values_all_samples.csv"), row.names = FALSE)

# Summary
summary_df <- data.frame(
  Factor = 1:model_opts$num_factors,
  Total_R2 = as.numeric(r2$R2Total[1, ])
)
for (view in views_names(MOFAobject.trained)) {
  summary_df[[paste0("R2_", view)]] <- as.numeric(r2$R2PerFactor[[view]][1, ])
}
write.csv(summary_df, file.path(OUTPUT_DIR, "factor_variance_summary.csv"), row.names = FALSE)

# Save trained model
saveRDS(MOFAobject.trained, file.path(OUTPUT_DIR, "brain_mofa_results.rds"))
saveRDS(r2, file.path(OUTPUT_DIR, "variance_explained.rds"))

cat("  ✓ All results saved to:", OUTPUT_DIR, "\n")

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║              Brain MOFA Pipeline Complete!                 ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat("\nTop 3 factors by variance explained:\n")
for (i in 1:3) {
  cat("  Factor", i, ":", round(summary_df$Total_R2[i] * 100, 2), "%\n")
}
cat("\nOutput files in:", OUTPUT_DIR, "\n")

invisible(list(mofa = MOFAobject.trained, r2 = r2))
