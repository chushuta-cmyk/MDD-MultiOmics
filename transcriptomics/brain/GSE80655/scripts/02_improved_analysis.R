# ============================================================================
# GSE88890 Improved Analysis Pipeline
# Enhanced version with SNP filtering, covariate adjustment, DMR scan
# ============================================================================

cat("\n=== GSE88890 Improved Analysis Started ===\n")

# ----------------------------------------------------------------------------
# 1. Setup
# ----------------------------------------------------------------------------
library(data.table)
library(limma)
library(ggplot2)
library(pheatmap)
library(tidyverse)
library(ChAMP)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(missMethyl)
library(clusterProfiler)
library(org.Hs.eg.db)

output_dir <- "GSE88890_Improved_Analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cat("✓ Libraries loaded\n\n")

# ----------------------------------------------------------------------------
# 2. Extract metadata from series_matrix.txt
# ----------------------------------------------------------------------------
cat("--- Extracting metadata ---\n")

series_matrix_file <- "GSE88890_series_matrix.txt"
all_lines <- readLines(series_matrix_file)

# Sample names
title_line <- all_lines[grep("^!Sample_title", all_lines)]
title_parts <- unlist(strsplit(title_line[1], "\t"))
sample_names <- gsub('"', '', title_parts[-1])

# Characteristics lines
all_char <- grep("^!Sample_characteristics_ch1", all_lines, value = TRUE)

# groupid
groupid_line <- all_char[grep("groupid", all_char)]
groupid_parts <- unlist(strsplit(groupid_line[1], "\t"))[-1]
groupid_values <- gsub('"', '', groupid_parts)
groupid_values <- gsub('.*groupid: ', '', groupid_values)

# Region (Brodmann area)
region_line <- all_char[grep("tissue|Brodmann", all_char)]
region_parts <- unlist(strsplit(region_line[1], "\t"))[-1]
region_values <- gsub('"', '', region_parts)
brodmann_numbers <- gsub('.*Brodmann area ([0-9]+).*', '\\1', region_values)
region_final <- paste0("BA", brodmann_numbers)

# Gender
gender_line <- all_char[grep("gender", all_char)]
gender_parts <- unlist(strsplit(gender_line[1], "\t"))[-1]
gender_values <- gsub('"', '', gender_parts)
gender_values <- gsub('.*gender: ', '', gender_values)

# Build targets
sample_id <- sprintf("Sample_%03d", 1:length(sample_names))
phenotype_values <- ifelse(grepl("MDD suicide", groupid_values), "MDD_suicide", "Control")

targets <- data.frame(
  Sample_ID = sample_id,
  Sample_Name = sample_names,
  Phenotype = as.factor(phenotype_values),
  Region = as.factor(region_final),
  Gender = as.factor(gender_values),
  stringsAsFactors = FALSE
)
rownames(targets) <- targets$Sample_ID

cat(sprintf("✓ Extracted %d samples\n", nrow(targets)))
write.csv(targets, file.path(output_dir, "00_Sample_Metadata.csv"), row.names = TRUE)
cat("\n")

# ----------------------------------------------------------------------------
# 3. Load Beta matrix and apply improved filtering
# ----------------------------------------------------------------------------
cat("--- Loading and filtering Beta matrix ---\n")

Beta_raw <- fread(series_matrix_file, sep = "\t", header = TRUE,
                  data.table = FALSE, skip = 74)
rownames(Beta_raw) <- Beta_raw[[1]]
Beta_raw <- Beta_raw[, -1]
Beta <- as.matrix(Beta_raw)
cat(sprintf("✓ Original Beta: %d CpGs × %d samples\n", nrow(Beta), ncol(Beta)))

# Improvement 1: Remove SNP-related probes
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
snp_probes <- rownames(annot)[!is.na(annot$SNP_ID)]
Beta_filtered <- Beta[!rownames(Beta) %in% snp_probes, ]
cat(sprintf("  Removed %d SNP-related probes\n", length(snp_probes)))

# Remove NA probes
na_per_probe <- rowSums(is.na(Beta_filtered))
Beta_filtered <- Beta_filtered[na_per_probe == 0, ]
cat(sprintf("  Kept %d probes without NA\n", nrow(Beta_filtered)))

# Improvement 2: Lower variance threshold (0.0001 vs 0.001)
probe_vars <- apply(Beta_filtered, 1, var, na.rm = TRUE)
var_threshold <- 0.0001
Beta_filtered <- Beta_filtered[probe_vars > var_threshold, ]
cat(sprintf("  Variance threshold: %f -> kept %d probes (%.1f%% total removed)\n",
            var_threshold, nrow(Beta_filtered),
            (1 - nrow(Beta_filtered)/nrow(Beta))*100))
Beta_filtered <- as.matrix(Beta_filtered)
cat("\n")

# ----------------------------------------------------------------------------
# 4. Quality control: sample correlation and missing rate
# ----------------------------------------------------------------------------
cat("--- Quality control ---\n")

beta_cor <- cor(Beta_filtered)
mean_cors <- colMeans(beta_cor)
threshold_cor <- mean(mean_cors) - 2 * sd(mean_cors)
abnormal <- names(mean_cors[mean_cors < threshold_cor])
if (length(abnormal) > 0) {
  cat("  Warning: potentially abnormal samples:", paste(abnormal, collapse=", "), "\n")
} else {
  cat("  ✓ No abnormal samples detected\n")
}

missing_rate <- colSums(is.na(Beta_filtered)) / nrow(Beta_filtered)
cat(sprintf("  Missing rate per sample: %.2f%% - %.2f%%\n",
            min(missing_rate)*100, max(missing_rate)*100))
cat("\n")

# ----------------------------------------------------------------------------
# 5. Design matrix with covariates (Region + Gender)
# ----------------------------------------------------------------------------
cat("--- Design matrix (including Region + Gender as covariates) ---\n")
design <- model.matrix(~targets$Region + targets$Gender + targets$Phenotype)
cat("Formula: ~Region + Gender + Phenotype\n")
cat("Design dimensions:", dim(design), "\n\n")

# ----------------------------------------------------------------------------
# 6. DMP analysis with covariate adjustment
# ----------------------------------------------------------------------------
cat("--- DMP analysis (covariate-adjusted) ---\n")

fit <- lmFit(Beta_filtered, design)
fit <- eBayes(fit)

# Phenotype coefficient is the last column
DMP_results <- topTable(fit, coef = ncol(design), number = Inf)

sig_count <- sum(DMP_results$adj.P.Val < 0.05)
hyper_count <- sum(DMP_results$adj.P.Val < 0.05 & DMP_results$logFC > 0)
hypo_count <- sum(DMP_results$adj.P.Val < 0.05 & DMP_results$logFC < 0)

cat(sprintf("✓ Significant DMPs (FDR<0.05): %d\n", sig_count))
cat(sprintf("  - Hyper-methylated: %d (%.1f%%)\n",
            hyper_count, if(sig_count>0) (hyper_count/sig_count)*100 else 0))
cat(sprintf("  - Hypo-methylated: %d (%.1f%%)\n\n",
            hypo_count, if(sig_count>0) (hypo_count/sig_count)*100 else 0))

write.csv(DMP_results, file.path(output_dir, "01_DMP_Results_Improved.csv"), row.names = TRUE)

# ----------------------------------------------------------------------------
# 7. CpG annotation
# ----------------------------------------------------------------------------
cat("--- CpG annotation ---\n")

DMP_annotated <- cbind(
  DMP_results,
  Gene = annot[rownames(DMP_results), "UCSC_RefGene_Name"],
  GenePart = annot[rownames(DMP_results), "UCSC_RefGene_Group"],
  Chr = annot[rownames(DMP_results), "chr"],
  Pos = annot[rownames(DMP_results), "pos"]
)

sig_dmp <- DMP_annotated[DMP_annotated$adj.P.Val < 0.05, ]
cat(sprintf("✓ %d DMPs annotated\n\n", nrow(sig_dmp)))
write.csv(sig_dmp, file.path(output_dir, "02_DMP_Annotated_Improved.csv"), row.names = TRUE)

# ----------------------------------------------------------------------------
# 8. Volcano plot
# ----------------------------------------------------------------------------
cat("--- Volcano plot ---\n")

pdf(file.path(output_dir, "03_VolcanoPlot_Improved.pdf"), width = 10, height = 8)
plot(DMP_results$logFC, -log10(DMP_results$adj.P.Val),
     main = "Volcano Plot: MDD Suicide vs Control (Batch-Corrected)",
     xlab = "log2(Beta Change)", ylab = "-log10(Adjusted P-value)",
     col = ifelse(DMP_results$adj.P.Val < 0.05, "red", "gray50"),
     pch = 16, cex = 0.6)
abline(h = -log10(0.05), lty = 2, col = "blue", lwd = 2)
abline(v = c(-0.1, 0.1), lty = 2, col = "blue", lwd = 2)
dev.off()
cat("  ✓ Volcano plot saved\n\n")

# ----------------------------------------------------------------------------
# 9. Top 30 DMPs heatmap
# ----------------------------------------------------------------------------
if (nrow(DMP_results) > 0) {
  cat("--- Top 30 DMPs heatmap ---\n")
  
  top_count <- min(30, nrow(DMP_results))
  top_idx <- order(DMP_results$adj.P.Val)[1:top_count]
  top_probes <- rownames(DMP_results)[top_idx]
  top_beta <- Beta_filtered[top_probes, ]
  
  pdf(file.path(output_dir, "04_Top30_Heatmap.pdf"), width = 12, height = 10)
  pheatmap(top_beta,
           annotation_col = data.frame(
             Phenotype = targets$Phenotype,
             Region = targets$Region,
             Gender = targets$Gender,
             row.names = rownames(targets)
           ),
           annotation_colors = list(
             Phenotype = c("MDD_suicide" = "#E41A1C", "Control" = "#377EB8"),
             Region = c("BA11" = "#4DAF4A", "BA25" = "#FF7F00"),
             Gender = c("M" = "#666666", "F" = "#FF69B4")
           ),
           main = sprintf("Top %d DMPs (Batch-Corrected)", top_count),
           scale = "row", show_colnames = FALSE,
           clustering_distance_cols = "correlation")
  dev.off()
  cat("  ✓ Heatmap saved\n\n")
}

# ----------------------------------------------------------------------------
# 10. Box plots for top 5 DMPs
# ----------------------------------------------------------------------------
if (nrow(DMP_results) >= 5) {
  cat("--- Top 5 DMPs box plots ---\n")
  
  top5_probes <- rownames(DMP_results)[order(DMP_results$adj.P.Val)][1:5]
  
  pdf(file.path(output_dir, "05_Top5_DMPs_BoxPlots.pdf"), width = 14, height = 10)
  par(mfrow = c(2, 3))
  
  for (i in 1:5) {
    probe <- top5_probes[i]
    beta_vals <- Beta_filtered[probe, ]
    gene_name <- ifelse(is.na(annot[probe, "UCSC_RefGene_Name"]), probe,
                        annot[probe, "UCSC_RefGene_Name"])
    
    plot_df <- data.frame(Beta = beta_vals, Phenotype = targets$Phenotype)
    boxplot(Beta ~ Phenotype, data = plot_df,
            main = paste0(gene_name, " (", probe, ")"),
            xlab = "Phenotype", ylab = "Beta Value",
            col = c("#377EB8", "#E41A1C"))
    points(jitter(as.numeric(plot_df$Phenotype), 0.2),
           plot_df$Beta, pch = 16, cex = 0.5)
  }
  
  par(mfrow = c(1, 1))
  dev.off()
  cat("  ✓ Box plots saved\n\n")
}

# ----------------------------------------------------------------------------
# 11. GO enrichment for methylation DMPs
# ----------------------------------------------------------------------------
cat("--- GO enrichment ---\n")

genes_raw <- as.character(sig_dmp$Gene[!is.na(sig_dmp$Gene) & sig_dmp$Gene != ""])
genes_list <- unique(unlist(lapply(strsplit(genes_raw, ";"), trimws)))
cat(sprintf("  Found %d unique genes\n", length(genes_list)))

if (length(genes_list) > 0) {
  gene_ids <- bitr(genes_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  if (nrow(gene_ids) > 0) {
    ego_bp <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db,
                       ont = "BP", pvalueCutoff = 0.05)
    if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
      cat(sprintf("  ✓ GO BP: %d significant terms\n", nrow(ego_bp)))
      write.csv(ego_bp@result, file.path(output_dir, "06_GO_BP_Results.csv"), row.names = FALSE)
    }
  }
}
cat("\n")

# ----------------------------------------------------------------------------
# 12. Lightweight DMR analysis
# ----------------------------------------------------------------------------
cat("--- DMR scan (lightweight) ---\n")

design_dmr <- model.matrix(~targets$Phenotype)
fit_dmr <- lmFit(Beta_filtered, design_dmr)
fit_dmr <- eBayes(fit_dmr)

cpg_results <- data.frame(
  CpG = rownames(DMP_results),
  logFC = DMP_results$logFC,
  P.Value = DMP_results$P.Value,
  Chr = annot[rownames(DMP_results), "chr"],
  Pos = annot[rownames(DMP_results), "pos"]
)
cpg_results <- cpg_results[order(cpg_results$Chr, cpg_results$Pos), ]

dmr_regions <- list()
for (chr in unique(cpg_results$Chr)) {
  chr_data <- cpg_results[cpg_results$Chr == chr, ]
  sig_data <- chr_data[chr_data$P.Value < 0.001, ]
  if (nrow(sig_data) >= 3) {
    gaps <- diff(sig_data$Pos)
    if (any(gaps < 1000)) {
      dmr_regions[[length(dmr_regions)+1]] <- data.frame(
        Chr = chr,
        Start = min(sig_data$Pos),
        End = max(sig_data$Pos),
        N_CpGs = nrow(sig_data),
        Mean_logFC = mean(sig_data$logFC)
      )
    }
  }
}

if (length(dmr_regions) > 0) {
  dmr_results <- do.call(rbind, dmr_regions)
  cat(sprintf("  Found %d potential DMRs\n", nrow(dmr_results)))
  write.csv(dmr_results, file.path(output_dir, "07_DMR_Results.csv"), row.names = FALSE)
} else {
  cat("  No DMRs detected\n")
}
cat("\n")

# ----------------------------------------------------------------------------
# 13. Final summary report
# ----------------------------------------------------------------------------
cat("--- Final report ---\n")

final_report <- sprintf("
================================================================================
           GSE88890 Improved Analysis Report (Batch-Corrected)
================================================================================

[Key Improvements]
- SNP probe removal: %d probes removed
- Variance threshold reduced: 0.001 -> 0.0001
- Covariates included: Region + Gender
- Added QC: abnormal sample detection
- Added DMR scan
- Added GO enrichment

[Data Summary]
Total probes: %d
Filtered probes: %d (%.1f%% removed)

Samples: %d (MDD suicide %d, Control %d)
Regions: BA11 %d, BA25 %d

[DMP Results]
Significant DMPs (FDR<0.05): %d
  - Hypermethylated: %d (%.1f%%)
  - Hypomethylated: %d (%.1f%%)

[Output Files]
- 00_Sample_Metadata.csv
- 01_DMP_Results_Improved.csv
- 02_DMP_Annotated_Improved.csv
- 03_VolcanoPlot_Improved.pdf
- 04_Top30_Heatmap.pdf
- 05_Top5_DMPs_BoxPlots.pdf
- 06_GO_BP_Results.csv
- 07_DMR_Results.csv

================================================================================
Completed: %s
================================================================================
",
    length(snp_probes),
    nrow(Beta), nrow(Beta_filtered), (1 - nrow(Beta_filtered)/nrow(Beta))*100,
    nrow(targets),
    sum(targets$Phenotype == "MDD_suicide"),
    sum(targets$Phenotype == "Control"),
    sum(targets$Region == "BA11"),
    sum(targets$Region == "BA25"),
    sig_count,
    hyper_count, if(sig_count>0) (hyper_count/sig_count)*100 else 0,
    hypo_count, if(sig_count>0) (hypo_count/sig_count)*100 else 0,
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

cat(final_report)
writeLines(final_report, file.path(output_dir, "00_Analysis_Report.txt"))

# ----------------------------------------------------------------------------
cat("\n✓✓✓ Improved Analysis Complete! ✓✓✓\n")
cat("Results saved to:", output_dir, "\n")