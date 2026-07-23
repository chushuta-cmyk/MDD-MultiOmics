# ==============================================================================
# MASTER MOFA PIPELINE v3 - FINAL VERSION
# Combines all datasets into MOFA-ready input
# ==============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════╗\n")
cat("║     MDD MULTI-OMICS MOFA+ PREPROCESSING PIPELINE         ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n\n")

# Setup
alt_lib <- "/Users/tan/Developer/envs/R/lib/R/library"
if (dir.exists(alt_lib)) .libPaths(c(alt_lib, .libPaths()))

library(MOFA2)
library(data.table)
library(dplyr)

cat("✓ MOFA2 loaded, version:", as.character(packageVersion("MOFA2")), "\n\n")

input_dir <- "/Users/tan/Desktop/MDD/next_steps/MOFA/inputs"
output_dir <- "/Users/tan/Desktop/MDD/next_steps/MOFA/outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# Step 1: Load available datasets
# ==============================================================================

cat("【Step 1】Loading preprocessed datasets...\n\n")

datasets <- list()
dataset_info <- list()

# GSE80655 - Brain RNA-seq
f <- file.path(input_dir, "GSE80655_expression_cleaned.rds")
if (file.exists(f)) {
  datasets$GSE80655_RNA <- readRDS(f)
  dataset_info$GSE80655 <- list(genes = nrow(datasets$GSE80655_RNA), samples = ncol(datasets$GSE80655_RNA), tissue = "Brain", type = "RNA-seq")
  cat(sprintf("  ✅ GSE80655 (Brain RNA-seq): %d genes × %d samples\n", 
              nrow(datasets$GSE80655_RNA), ncol(datasets$GSE80655_RNA)))
} else {
  cat("  ❌ GSE80655 not found\n")
}

# GSE98793 - Blood Microarray
f <- file.path(input_dir, "GSE98793_expression_cleaned.rds")
if (file.exists(f)) {
  datasets$GSE98793_RNA <- readRDS(f)
  dataset_info$GSE98793 <- list(genes = nrow(datasets$GSE98793_RNA), samples = ncol(datasets$GSE98793_RNA), tissue = "Blood", type = "Microarray")
  cat(sprintf("  ✅ GSE98793 (Blood Microarray): %d genes × %d samples\n",
              nrow(datasets$GSE98793_RNA), ncol(datasets$GSE98793_RNA)))
} else {
  cat("  ❌ GSE98793 not found\n")
}

# GSE54563 - Brain Microarray (ACC)
f <- file.path(input_dir, "GSE54563_expression_cleaned.rds")
if (file.exists(f)) {
  datasets$GSE54563_RNA <- readRDS(f)
  dataset_info$GSE54563 <- list(genes = nrow(datasets$GSE54563_RNA), samples = ncol(datasets$GSE54563_RNA), tissue = "Brain_ACC", type = "Microarray")
  cat(sprintf("  ✅ GSE54563 (Brain ACC Microarray): %d genes × %d samples\n",
              nrow(datasets$GSE54563_RNA), ncol(datasets$GSE54563_RNA)))
} else {
  cat("  ❌ GSE54563 not found\n")
}

# GSE88890 - Brain Methylation
f <- file.path(input_dir, "GSE88890_methylation_gene_level.rds")
if (file.exists(f)) {
  datasets$GSE88890_Meth <- readRDS(f)
  dataset_info$GSE88890 <- list(genes = nrow(datasets$GSE88890_Meth), samples = ncol(datasets$GSE88890_Meth), tissue = "Brain", type = "Methylation")
  cat(sprintf("  ✅ GSE88890 (Brain Methylation): %d genes × %d samples\n",
              nrow(datasets$GSE88890_Meth), ncol(datasets$GSE88890_Meth)))
} else {
  cat("  ❌ GSE88890 not found\n")
}

if (length(datasets) == 0) {
  stop("No datasets found! Please run preprocessing scripts first.")
}

cat(sprintf("\nTotal datasets loaded: %d\n", length(datasets)))

# ==============================================================================
# Step 2: Strategy - Multi-Group MOFA
# ==============================================================================

cat("\n【Step 2】MOFA Strategy: Multi-Group (independent cohorts)\n\n")

cat("These datasets have NO shared samples (different cohorts).\n")
cat("Using MOFA2 multi-group design to extract cross-cohort factors.\n\n")

# Separate RNA and Methylation views
rna_views <- datasets[grepl("_RNA$", names(datasets))]
meth_views <- datasets[grepl("_Meth$", names(datasets))]

cat("RNA datasets:", length(rna_views), "\n")
for (name in names(rna_views)) {
  cat(sprintf("  - %s: %d genes\n", name, nrow(rna_views[[name]])))
}

cat("\nMethylation datasets:", length(meth_views), "\n")
for (name in names(meth_views)) {
  cat(sprintf("  - %s: %d genes\n", name, nrow(meth_views[[name]])))
}

# ==============================================================================
# Step 3: Build MOFA data - Multi-Group Format
# ==============================================================================

cat("\n【Step 3】Building MOFA multi-group data structure...\n\n")

# For multi-group MOFA, data format is:
# list(group1 = list(view1 = matrix, view2 = matrix), 
#      group2 = list(view1 = matrix, view2 = matrix), ...)

# But since our datasets don't share samples, we use a simplified approach:
# Create a single group with all samples, but views only have data for their respective samples

groups <- list()

# Process each RNA dataset as a separate group
for (name in names(rna_views)) {
  group_name <- gsub("_RNA$", "", name)
  mat <- rna_views[[name]]

  # Prefix feature names to avoid conflicts
  rownames(mat) <- paste0("GEX_", rownames(mat))

  groups[[group_name]] <- list(RNA = mat)
  cat(sprintf("  Group '%s': %d RNA features × %d samples\n", group_name, nrow(mat), ncol(mat)))
}

# Add methylation data to groups that have matching samples
# Since GSE88890 is independent, add as separate group
for (name in names(meth_views)) {
  group_name <- gsub("_Meth$", "", name)
  mat <- meth_views[[name]]

  # Prefix feature names
  rownames(mat) <- paste0("MET_", rownames(mat))

  if (group_name %in% names(groups)) {
    groups[[group_name]]$Methylation <- mat
    cat(sprintf("  Group '%s': Added Methylation %d features × %d samples\n", group_name, nrow(mat), ncol(mat)))
  } else {
    groups[[group_name]] <- list(Methylation = mat)
    cat(sprintf("  Group '%s': %d Methylation features × %d samples\n", group_name, nrow(mat), ncol(mat)))
  }
}

cat(sprintf("\nTotal groups: %d\n", length(groups)))
for (g in names(groups)) {
  views <- names(groups[[g]])
  cat(sprintf("  %s: %s\n", g, paste(views, collapse = ", ")))
}

# ==============================================================================
# Step 4: Create MOFA object
# ==============================================================================

cat("\n【Step 4】Creating MOFA object...\n")

# Use create_mofa (not create_mofa_object)
MOFAobject <- create_mofa(groups)

cat(sprintf("\n  ✅ MOFA object created\n"))
cat(sprintf("  Groups: %s\n", paste(names(groups), collapse = ", ")))
cat(sprintf("  Views: %s\n", paste(MOFAobject@dimensions$M, collapse = ", ")))
cat(sprintf("  Total samples: %d\n", MOFAobject@dimensions$N))

# ==============================================================================
# Step 5: Configure and train
# ==============================================================================

cat("\n【Step 5】Configuring MOFA...\n")

opts <- get_default_training_options(MOFAobject)
opts$convergence_mode <- "fast"
opts$maxiter <- 1000
opts$verbose <- TRUE

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

# Handle missing values (sparse data expected)
model_opts$spikeslab_factors <- FALSE
model_opts$spikeslab_weights <- FALSE

MOFAobject <- prepare_mofa(MOFAobject,
                           model_options = model_opts,
                           training_options = opts)

cat("  ✅ MOFA configured\n")
cat(sprintf("  Factors: %d\n", model_opts$num_factors))
cat(sprintf("  Convergence: %s\n", opts$convergence_mode))

# ==============================================================================
# Step 6: Train
# ==============================================================================

cat("\n【Step 6】Training MOFA...\n")
cat("  This may take 10-30 minutes...\n\n")

MOFAmodel <- run_mofa(MOFAobject, outfile = file.path(output_dir, "mofa_model.hdf5"))

cat("\n  ✅ MOFA training complete!\n")

# ==============================================================================
# Step 7: Extract results
# ==============================================================================

cat("\n【Step 7】Extracting and saving results...\n")

# Variance explained
r2 <- get_variance_explained(MOFAmodel)
saveRDS(r2, file.path(output_dir, "variance_explained.rds"))
write.csv(r2$r2_per_factor, file.path(output_dir, "variance_per_factor.csv"))
cat("  ✅ Variance explained\n")

# Factor scores
factors <- get_factors(MOFAmodel)
saveRDS(factors, file.path(output_dir, "factor_scores.rds"))
for (g in names(factors)) {
  write.csv(factors[[g]], file.path(output_dir, sprintf("factor_scores_%s.csv", g)))
}
cat("  ✅ Factor scores\n")

# Weights
weights <- get_weights(MOFAmodel)
for (view in names(weights)) {
  w <- weights[[view]]
  write.csv(w, file.path(output_dir, sprintf("weights_%s.csv", view)))
}
cat("  ✅ Feature weights\n")

# ==============================================================================
# Step 8: Plots
# ==============================================================================

cat("\n【Step 8】Generating plots...\n")

# Variance explained
pdf(file.path(output_dir, "variance_explained.pdf"), width = 10, height = 6)
plot_variance_explained(MOFAmodel)
dev.off()
cat("  ✅ variance_explained.pdf\n")

# Factor correlation
tryCatch({
  pdf(file.path(output_dir, "factor_correlation.pdf"), width = 8, height = 8)
  plot_factor_cor(MOFAmodel)
  dev.off()
  cat("  ✅ factor_correlation.pdf\n")
}, error = function(e) {
  cat("  ⚠ Factor correlation plot skipped\n")
})

# ==============================================================================
# Summary
# ==============================================================================

cat("\n╔════════════════════════════════════════════════════════════╗\n")
cat("║                  MOFA ANALYSIS COMPLETE                    ║\n")
cat("╚════════════════════════════════════════════════════════════╝\n")
cat(sprintf("\nModel: %s\n", file.path(output_dir, "mofa_model.hdf5")))
cat(sprintf("Factors: %d\n", model_opts$num_factors))
cat(sprintf("Groups: %s\n", paste(names(groups), collapse = ", ")))
cat(sprintf("Total samples: %d\n", MOFAobject@dimensions$N))
cat("\nDatasets used:\n")
for (name in names(dataset_info)) {
  info <- dataset_info[[name]]
  cat(sprintf("  %s: %s, %s, %d genes × %d samples\n", 
              name, info$tissue, info$type, info$genes, info$samples))
}
cat("\nAll outputs saved to:", output_dir, "\n")
cat("\nKey files:\n")
cat("  - mofa_model.hdf5\n")
cat("  - variance_explained.pdf / .csv\n")
cat("  - factor_scores_*.csv\n")
cat("  - weights_*.csv\n")
cat("═══════════════════════════════════════════════════════════════\n")
