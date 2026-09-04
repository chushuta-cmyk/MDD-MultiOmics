# ==============================================================================
# GSE53987: Postmortem Brain Transcriptomic Analysis in MDD
# ==============================================================================
# Platform: Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
# Regions: Hippocampus, Pre-frontal Cortex BA46, Associative Striatum
# Groups: MDD vs. Controls
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------------------
library(stringr)
library(GEOquery)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(hgu133plus2.db)

# --- Paths (modify as needed) ---
data_dir   <- "data"
output_dir <- "results"
script_dir <- "."

file_path <- file.path(data_dir, "GSE53987_series_matrix.txt")
helper_file <- file.path(script_dir, "plot_generator.R")

if (!file.exists(file_path)) stop("Data file not found: ", file_path)
if (!file.exists(helper_file)) stop("Helper script not found: ", helper_file)

source(helper_file)

# ------------------------------------------------------------------------------
# 2. Data Loading
# ------------------------------------------------------------------------------
raw_data <- read.delim(
  file = file_path,
  skip = 66,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  check.names = FALSE,
  as.is = TRUE
)

if (!"ID_REF" %in% colnames(raw_data)) stop("'ID_REF' column not found.")

rownames(raw_data) <- raw_data[, "ID_REF"]
expr_matrix <- as.matrix(raw_data[, -which(colnames(raw_data) == "ID_REF")])
storage.mode(expr_matrix) <- "numeric"

# ------------------------------------------------------------------------------
# 3. Log2 Transformation (auto-detect)
# ------------------------------------------------------------------------------
if (max(expr_matrix, na.rm = TRUE) > 50) {
  expr_matrix <- log2(expr_matrix + 1)
  cat("✓ Raw linear data detected. Log2 transformation applied.\n")
}

# ------------------------------------------------------------------------------
# 4. Quantile Normalization
# ------------------------------------------------------------------------------
cat("Running quantile normalization...\n")

pdf(file.path(output_dir, "QC_Boxplot_BeforeNorm.pdf"), width = 10, height = 5)
boxplot(expr_matrix, outline = FALSE, las = 2, main = "Before normalization")
dev.off()

expr_matrix <- normalizeBetweenArrays(expr_matrix, method = "quantile")

pdf(file.path(output_dir, "QC_Boxplot_AfterNorm.pdf"), width = 10, height = 5)
boxplot(expr_matrix, outline = FALSE, las = 2, main = "After quantile normalization")
dev.off()

pdf(file.path(output_dir, "QC_DensityPlot.pdf"), width = 8, height = 5)
plotDensities(expr_matrix, main = "Expression density")
dev.off()

cat("✓ Quantile normalization and QC plots complete.\n")

# ------------------------------------------------------------------------------
# 5. Expression Filtering
# ------------------------------------------------------------------------------
# Remove probes with any NA
expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0, ]

# Remove low-expression probes (rowMeans < median)
median_cutoff <- median(expr_matrix)
expr_matrix <- expr_matrix[rowMeans(expr_matrix) > median_cutoff, ]

# Remove low-variance probes (below 25th percentile)
probe_var <- apply(expr_matrix, 1, var)
var_cutoff <- quantile(probe_var, 0.25)
expr_matrix <- expr_matrix[probe_var > var_cutoff, ]

cat("Filtered matrix:", dim(expr_matrix), "\n")

# ------------------------------------------------------------------------------
# 6. Phenotype Parsing
# ------------------------------------------------------------------------------
geo_lines <- readLines(file_path)
pheno_line <- grep("^!Sample_source_name_ch1", geo_lines, value = TRUE)
if (length(pheno_line) == 0) stop("Metadata line not found.")

labels <- gsub('"', '', strsplit(pheno_line, "\t")[[1]][-1])
keep <- grepl("major depressive disorder|control", labels, ignore.case = TRUE)
selected_labels <- labels[keep]

split_labels <- str_split_fixed(selected_labels, ", ", n = 2)
brain_area <- split_labels[, 1]
disease <- split_labels[, 2]

simple <- disease == ""
brain_area[simple] <- "hippocampus"
disease[simple] <- split_labels[simple, 1]

# Clean brain area names
brain_area <- gsub(" ", "_", brain_area)
brain_area <- gsub("\\(|\\)", "", brain_area)
brain_area[grepl("Pre-frontal_cortex", brain_area)] <- "Pre-frontal_cortex_BA46"
brain_area[grepl("Associative_striatum", brain_area)] <- "Associative_striatum"
brain_area[grepl("hippocampus", brain_area)] <- "hippocampus"

group <- factor(paste(brain_area, disease, sep = "_"))
levels(group) <- gsub(" ", "_", levels(group))
levels(group) <- gsub("\\(|\\)", "", levels(group))

# Subset expression matrix
expr_data <- expr_matrix[, keep]
all_gsm <- colnames(expr_data)

# Build metadata
sample_metadata <- data.frame(
  GSM_ID = all_gsm,
  Group = group,
  Brain_Area = brain_area,
  Disease = disease,
  stringsAsFactors = FALSE
)

regions <- c("hippocampus", "Pre-frontal_cortex_BA46", "Associative_striatum")

# ------------------------------------------------------------------------------
# 7. Region-Specific Differential Expression (Limma)
# ------------------------------------------------------------------------------
deg_results <- list()

for (region in regions) {
  cat("\n--- Processing region:", region, "---\n")
  
  # Subset samples for this region
  gsm_ids <- sample_metadata$GSM_ID[sample_metadata$Brain_Area == region]
  region_expr <- expr_data[, gsm_ids, drop = FALSE]
  region_meta <- sample_metadata[sample_metadata$Brain_Area == region, ]
  
  # Align order
  region_expr <- region_expr[, region_meta$GSM_ID, drop = FALSE]
  
  # Build design matrix
  region_group <- factor(region_meta$Group, levels = unique(region_meta$Group))
  valid_levels <- make.names(levels(region_group))
  region_group_valid <- factor(region_group, labels = valid_levels)
  
  design <- model.matrix(~0 + region_group_valid)
  colnames(design) <- gsub("region_group_valid", "", colnames(design))
  
  # Fit model
  fit <- lmFit(region_expr, design)
  
  # Define contrast
  mdd_name <- make.names(paste0(region, "_major_depressive_disorder"))
  ctrl_name <- make.names(paste0(region, "_control"))
  
  if (!all(c(mdd_name, ctrl_name) %in% colnames(design))) {
    warning("Contrast groups not found for region:", region)
    next
  }
  
  contrast <- makeContrasts(contrasts = paste(mdd_name, ctrl_name, sep = " - "),
                            levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  
  # Extract results (P < 0.01)
  tt <- topTable(fit2, coef = 1, number = Inf, adjust.method = "fdr")
  tt <- tt[tt$P.Value < 0.01, ]
  
  if (nrow(tt) == 0) {
    deg_results[[region]] <- data.frame()
    next
  }
  
  # Probe to Gene Symbol annotation
  probe_ids <- rownames(tt)
  symbols <- suppressWarnings(
    mapIds(hgu133plus2.db, keys = probe_ids, column = "SYMBOL",
           keytype = "PROBEID", multiVals = "first")
  )
  
  # Fallback: uppercase conversion
  if (sum(!is.na(symbols)) < 10) {
    symbols_alt <- suppressWarnings(
      mapIds(hgu133plus2.db, keys = toupper(probe_ids), column = "SYMBOL",
             keytype = "PROBEID", multiVals = "first")
    )
    if (sum(!is.na(symbols_alt)) > sum(!is.na(symbols))) {
      symbols <- symbols_alt
    }
  }
  
  # Build result table
  res <- tt[, c("logFC", "P.Value", "adj.P.Val")]
  res$ProbeID <- probe_ids
  res$GeneSymbol <- symbols
  res$Region <- region
  res <- res[!is.na(res$GeneSymbol), ]
  
  # Deduplicate: keep probe with highest |logFC|
  res <- res[order(res$GeneSymbol, -abs(res$logFC)), ]
  res <- res[!duplicated(res$GeneSymbol), ]
  
  deg_results[[region]] <- res
  
  cat("  ✓ DEGs:", nrow(res), "\n")
}

# ------------------------------------------------------------------------------
# 8. Merge and Export Results
# ------------------------------------------------------------------------------
all_degs <- do.call(rbind, deg_results)

output_file <- file.path(output_dir, "GSE53987_MDD_vs_Control_BrainRegion_DEGs.csv")
write.csv(all_degs, file = output_file, row.names = FALSE)
cat("\n✓ Results saved to:", output_file, "\n")

# ------------------------------------------------------------------------------
# 9. Visualization (requires plot_generator.R)
# ------------------------------------------------------------------------------
if (exists("plot_volcano") && exists("plot_deg_heatmap")) {
  
  # Volcano plot per region
  for (region in regions) {
    region_degs <- all_degs[all_degs$Region == region, ]
    if (nrow(region_degs) > 0) {
      plot_volcano(
        deg_df = region_degs,
        p_cutoff = 0.05,
        logfc_cutoff = 0.05,
        plot_title = paste(region, "Transcriptome Dysregulation"),
        output_prefix = file.path(output_dir, paste0("volcano_", region))
      )
    }
  }
  
  # Volcano plot: all regions combined
  plot_volcano(
    deg_df = all_degs,
    p_cutoff = 0.05,
    logfc_cutoff = 0.05,
    plot_title = "GSE53987 Entire Dataset Profile",
    output_prefix = file.path(output_dir, "volcano_all_dataset")
  )
  
  # Heatmap per region (using each region's expression data)
  for (region in regions) {
    region_meta <- sample_metadata[sample_metadata$Brain_Area == region, ]
    region_expr <- expr_data[, region_meta$GSM_ID, drop = FALSE]
    region_group <- factor(region_meta$Group, levels = unique(region_meta$Group))
    region_degs <- all_degs[all_degs$Region == region, ]
    
    if (nrow(region_degs) > 0 && ncol(region_expr) > 2) {
      plot_deg_heatmap(
        deg_df = region_degs,
        expr_matrix = region_expr,
        group_factor = region_group,
        top_n = 40,
        plot_title = paste("Top 40 DEGs -", region),
        output_prefix = file.path(output_dir, paste0("heatmap_", region))
      )
    }
  }
  
  cat("✓ Visualization complete.\n")
  
} else {
  cat("⚠️ plot_generator.R not loaded or missing functions. Skipping visualization.\n")
}

# ------------------------------------------------------------------------------
# 10. Session Info (for reproducibility)
# ------------------------------------------------------------------------------
sink(file.path(output_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n✓ Analysis complete. Outputs saved to:", output_dir, "\n")