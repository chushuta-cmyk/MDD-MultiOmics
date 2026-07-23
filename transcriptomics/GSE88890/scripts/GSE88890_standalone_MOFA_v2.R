# ==============================================================================
# GSE88890 Methylation → MOFA Input (Gene-Level Beta)
# Platform: Illumina HumanMethylation450K
# Tissue: Brain (BA11, BA25)
# ==============================================================================
# Based on your GSE88890_Complete_Analysis_Final.R
# Modified to save gene-level beta matrix for MOFA integration
# ==============================================================================

cat("\n")
cat("====================================================================\n")
cat("    GSE88890 Methylation → MOFA Input (Gene-Level)\n")
cat("    Platform: Illumina 450K\n")
cat("====================================================================\n\n")

# ==============================================================================
# Step 0: Setup
# ==============================================================================

alt_lib <- "/Users/tan/Developer/envs/R/lib/R/library"
if (dir.exists(alt_lib)) {
  .libPaths(c(alt_lib, .libPaths()))
  cat("✓ Added library path\n")
}

library(data.table)
library(dplyr)
library(limma)
library(ggplot2)
library(pheatmap)

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19", ask = FALSE, update = FALSE)
}
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

cat("✓ Packages loaded\n\n")

# ==============================================================================
# Step 1: Extract metadata (from your working script)
# ==============================================================================

cat("【Step 1】Extracting metadata...\n")

setwd("/Users/tan/Developer/projects/r_project/Study/GSE88890/")

series_matrix_file <- "GSE88890_series_matrix.txt"

if (!file.exists(series_matrix_file)) {
  stop("File not found: ", series_matrix_file)
}

all_lines <- readLines(series_matrix_file)

# Extract sample titles
title_line <- all_lines[grep("^!Sample_title", all_lines)]
title_parts <- unlist(strsplit(title_line[1], "\t"))
sample_names <- gsub('"', '', title_parts[-1])

# Extract characteristics
all_char <- grep("^!Sample_characteristics_ch1", all_lines, value = TRUE)

# groupid
groupid_line <- all_char[grep("groupid", all_char)]
groupid_parts <- unlist(strsplit(groupid_line[1], "\t"))[-1]
groupid_values <- gsub('"', '', groupid_parts)
groupid_values <- gsub('.*groupid: ', '', groupid_values)

# region (Brodmann area)
region_line <- all_char[grep("tissue|Brodmann", all_char)]
region_parts <- unlist(strsplit(region_line[1], "\t"))[-1]
region_values <- gsub('"', '', region_parts)
brodmann_numbers <- gsub('.*Brodmann area ([0-9]+).*', '\\1', region_values)
region_final <- paste0("BA", brodmann_numbers)

# gender
gender_line <- all_char[grep("gender", all_char)]
gender_parts <- unlist(strsplit(gender_line[1], "\t"))[-1]
gender_values <- gsub('"', '', gender_parts)
gender_values <- gsub('.*gender: ', '', gender_values)

# Create targets
phenotype_values <- ifelse(grepl("MDD suicide", groupid_values), "MDD_suicide", "Control")

sample_id <- sprintf("Sample_%03d", 1:length(sample_names))

targets <- data.frame(
  Sample_ID = sample_id,
  Sample_Name = sample_names,
  Phenotype = phenotype_values,
  Region = region_final,
  Gender = gender_values,
  stringsAsFactors = FALSE
)
rownames(targets) <- targets$Sample_ID

cat(sprintf("  Total samples: %d\n", nrow(targets)))
cat("  Phenotype distribution:\n")
print(table(targets$Phenotype, targets$Region))

# ==============================================================================
# Step 2: Load Beta matrix (from series matrix, skip=74)
# ==============================================================================

cat("\n【Step 2】Loading Beta matrix...\n")

Beta_raw <- fread(series_matrix_file,
                   sep = "\t",
                   header = TRUE,
                   data.table = FALSE,
                   skip = 74)

rownames(Beta_raw) <- as.character(Beta_raw[[1]])
Beta_raw <- Beta_raw[, -1]
Beta <- as.matrix(Beta_raw)

cat(sprintf("  Beta matrix: %d CpGs × %d samples\n", nrow(Beta), ncol(Beta)))

# Basic QC
na_count <- sum(is.na(Beta))
cat(sprintf("  Missing values: %d (%.2f%%)\n", na_count, 100*na_count/length(Beta)))

# Filter: remove probes with NA, keep variable probes
na_per_probe <- rowSums(is.na(Beta))
Beta_qc <- Beta[na_per_probe == 0, ]

probe_vars <- apply(Beta_qc, 1, var, na.rm = TRUE)
Beta_qc <- Beta_qc[probe_vars > 0.001, ]

cat(sprintf("  After QC: %d CpGs\n", nrow(Beta_qc)))

# ==============================================================================
# Step 3: Map CpG to genes and aggregate
# ==============================================================================

cat("\n【Step 3】Mapping CpG probes to genes...\n")

# Get 450K annotation
cpg_annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Find overlapping probes
probes_in_data <- rownames(Beta_qc)
probes_in_anno <- rownames(cpg_annot)
overlap <- sum(probes_in_data %in% probes_in_anno)
cat(sprintf("  %d/%d probes found in 450K annotation\n", overlap, length(probes_in_data)))

# Subset to annotated probes
Beta_annotated <- Beta_qc[probes_in_data %in% probes_in_anno, ]
cpg_annot_sub <- cpg_annot[rownames(Beta_annotated), ]

# Get gene names (UCSC_RefGene_Name column)
gene_names <- cpg_annot_sub$UCSC_RefGene_Name

# Build probe-to-gene mapping (handle multiple genes per probe)
mapping_list <- list()
for (i in seq_along(gene_names)) {
  if (!is.na(gene_names[i]) && gene_names[i] != "") {
    genes <- unique(unlist(strsplit(as.character(gene_names[i]), ";")))
    genes <- genes[genes != ""]
    if (length(genes) > 0) {
      for (g in genes) {
        mapping_list[[length(mapping_list) + 1]] <- data.frame(
          ProbeID = rownames(Beta_annotated)[i],
          GeneSymbol = g,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

probe_gene_map <- do.call(rbind, mapping_list)
cat(sprintf("  %d probe-gene mappings created\n", nrow(probe_gene_map)))

# ==============================================================================
# Step 4: Aggregate to gene-level (mean beta per gene)
# ==============================================================================

cat("\n【Step 4】Aggregating to gene-level (mean beta per gene)...\n")

# Convert beta matrix to long format
beta_df <- as.data.frame(Beta_annotated)
beta_df$ProbeID <- rownames(Beta_annotated)

# Merge with gene mapping
beta_long <- merge(probe_gene_map, beta_df, by = "ProbeID")

# Calculate mean beta per gene per sample
# Exclude ProbeID and GeneSymbol columns for aggregation
sample_cols <- setdiff(colnames(beta_long), c("ProbeID", "GeneSymbol"))

gene_beta <- beta_long %>%
  group_by(GeneSymbol) %>%
  summarise(across(all_of(sample_cols), ~mean(., na.rm = TRUE)), .groups = "drop")

gene_beta_matrix <- as.matrix(gene_beta[, -1])
rownames(gene_beta_matrix) <- gene_beta$GeneSymbol

cat(sprintf("  Gene-level beta: %d genes × %d samples\n", 
            nrow(gene_beta_matrix), ncol(gene_beta_matrix)))

# ==============================================================================
# Step 5: Save for MOFA
# ==============================================================================

cat("\n【Step 5】Saving to MOFA inputs...\n")

output_dir <- "/Users/tan/Desktop/MDD/next_steps/MOFA/inputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

saveRDS(gene_beta_matrix, file.path(output_dir, "GSE88890_methylation_gene_level.rds"))
saveRDS(targets, file.path(output_dir, "GSE88890_metadata.rds"))
write.csv(gene_beta_matrix, file.path(output_dir, "GSE88890_methylation_gene_level.csv"))

cat("  ✓ Saved: GSE88890_methylation_gene_level.rds\n")
cat("  ✓ Saved: GSE88890_metadata.rds\n")
cat("  ✓ Saved: GSE88890_methylation_gene_level.csv\n")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n====================================================================\n")
cat("    GSE88890 → MOFA Input COMPLETE\n")
cat("====================================================================\n")
cat(sprintf("Input:  %d CpGs × %d samples\n", nrow(Beta), ncol(Beta)))
cat(sprintf("QC:     %d CpGs after filtering\n", nrow(Beta_qc)))
cat(sprintf("Output: %d genes × %d samples (mean beta per gene)\n", 
            nrow(gene_beta_matrix), ncol(gene_beta_matrix)))
cat(sprintf("Groups: %d MDD_suicide, %d Control\n", 
            sum(targets$Phenotype == "MDD_suicide"),
            sum(targets$Phenotype == "Control")))
cat("Files saved to:", output_dir, "\n")
cat("====================================================================\n")
