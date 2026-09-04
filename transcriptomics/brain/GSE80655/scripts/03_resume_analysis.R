# ============================================================================
# GSE88890 Resume / Patch-Run Analysis
# Skips completed steps, runs only missing analyses (DMR, GO/KEGG, etc.)
# ============================================================================

cat("\n=== GSE88890 Resume Analysis Started ===\n")

# ----------------------------------------------------------------------------
# 0. Configuration
# ----------------------------------------------------------------------------
input_data_dir  <- "./data"
input_plots_dir <- "./plots"
output_dir      <- "./results"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Input data directory: ", input_data_dir, "\n")
cat("Output directory: ", output_dir, "\n\n")

# ----------------------------------------------------------------------------
# 1. Setup: load libraries and check optional packages
# ----------------------------------------------------------------------------
library(data.table)
library(limma)
library(ggplot2)
library(pheatmap)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

has_missMethyl  <- requireNamespace("missMethyl", quietly = TRUE)
has_bumphunter  <- requireNamespace("bumphunter", quietly = TRUE)
has_DMRcate     <- requireNamespace("DMRcate", quietly = TRUE)

cat("✓ Core libraries loaded\n")
if (has_missMethyl)  cat("  - missMethyl available\n")
if (has_bumphunter)  cat("  - bumphunter available\n")
if (has_DMRcate)     cat("  - DMRcate available\n")
cat("\n")

# ----------------------------------------------------------------------------
# 2. Load Beta matrix (if available)
# ----------------------------------------------------------------------------
has_beta <- FALSE
Beta <- NULL
Beta_filtered <- NULL

beta_file <- file.path(input_data_dir, "GSE88890_normalisedBetas.csv")
if (file.exists(beta_file)) {
  cat("Loading Beta matrix from:", beta_file, "\n")
  Beta_raw <- fread(beta_file, sep = ",", header = TRUE, data.table = FALSE)
  rownames(Beta_raw) <- Beta_raw[[1]]
  Beta_raw <- Beta_raw[, -1]
  Beta <- as.matrix(Beta_raw)
  cat(sprintf("✓ Beta matrix loaded: %d CpGs × %d samples\n", nrow(Beta), ncol(Beta)))
  has_beta <- TRUE
} else {
  cat("✗ Beta matrix not found:", beta_file, "\n")
  cat("  Skipping steps that require raw Beta values.\n")
}

# ----------------------------------------------------------------------------
# 3. Load or create sample metadata
# ----------------------------------------------------------------------------
targets <- NULL

targets_file <- file.path(input_data_dir, "00_Sample_Metadata.csv")
if (file.exists(targets_file)) {
  cat("Loading sample metadata from:", targets_file, "\n")
  targets <- read.csv(targets_file, stringsAsFactors = FALSE)
  cat(sprintf("✓ Loaded metadata for %d samples\n", nrow(targets)))
} else if (has_beta) {
  cat("Metadata not found. Creating dummy targets from column names.\n")
  sample_names <- colnames(Beta)
  targets <- data.frame(Sample_Name = sample_names, stringsAsFactors = FALSE)
  targets$Region <- ifelse(
    grepl("_25$|25_|BA25", sample_names), "BA25",
    ifelse(grepl("_11$|11_|BA11", sample_names), "BA11", "Unknown")
  )
  targets$Phenotype <- "Unknown"
}

# Print summary if targets exist
if (!is.null(targets)) {
  cat(sprintf("✓ Total samples: %d\n", nrow(targets)))
  if ("Region" %in% colnames(targets)) {
    cat("  Regions:", paste(table(targets$Region), collapse = ", "), "\n")
  }
  if ("Phenotype" %in% colnames(targets)) {
    cat("  Groups:", paste(table(targets$Phenotype), collapse = ", "), "\n")
  }
  cat("\n")
}

# Ensure standard column names
if (!is.null(targets)) {
  colnames(targets) <- trimws(colnames(targets))
  
  # Sample column
  if (!"Sample_Name" %in% colnames(targets)) {
    cand <- c("Sample_Name", "sample", "Sample", "sample_id", "ID", "Array")
    hit <- cand[cand %in% colnames(targets)]
    if (length(hit) > 0) targets$Sample_Name <- targets[[hit[1]]]
  }
  
  # Phenotype column
  if (!"Phenotype" %in% colnames(targets)) {
    cand <- c("Phenotype", "Group", "group", "Diagnosis", "diagnosis", "Status", "status")
    hit <- cand[cand %in% colnames(targets)]
    if (length(hit) > 0) targets$Phenotype <- targets[[hit[1]]]
  }
  
  # Region column
  if (!"Region" %in% colnames(targets)) {
    cand <- c("Region", "region", "BrainRegion", "brain_region", "Area")
    hit <- cand[cand %in% colnames(targets)]
    if (length(hit) > 0) targets$Region <- targets[[hit[1]]]
  }
  
  # Clean phenotype
  targets$Phenotype <- trimws(as.character(targets$Phenotype))
  targets$Phenotype[targets$Phenotype == ""] <- NA
  targets$Phenotype[grepl("^control$", targets$Phenotype, ignore.case = TRUE)] <- "Control"
  targets$Phenotype[grepl("^mdd_suicide$", targets$Phenotype, ignore.case = TRUE)] <- "MDD_suicide"
  targets$Phenotype <- factor(targets$Phenotype, levels = c("Control", "MDD_suicide"))
  
  if (length(unique(na.omit(targets$Phenotype))) < 2) {
    stop("Phenotype has fewer than 2 valid groups after parsing.")
  }
}

# ----------------------------------------------------------------------------
# 4. Quality control and filtering (only if Beta available)
# ----------------------------------------------------------------------------
if (has_beta && !is.null(Beta)) {
  cat("\n--- Quality control and filtering ---\n")
  
  beta_min <- min(Beta, na.rm = TRUE)
  beta_max <- max(Beta, na.rm = TRUE)
  na_count <- sum(is.na(Beta))
  cat(sprintf("Beta range: [%.4f, %.4f]\n", beta_min, beta_max))
  cat(sprintf("Missing values: %d (%.2f%%)\n", na_count,
              (na_count / (nrow(Beta) * ncol(Beta))) * 100))
  
  na_per_probe <- rowSums(is.na(Beta))
  Beta_filtered <- Beta[na_per_probe == 0, ]
  probe_vars <- apply(Beta_filtered, 1, var, na.rm = TRUE)
  Beta_filtered <- Beta_filtered[probe_vars > 0.001, ]
  
  cat(sprintf("Original probes: %d\n", nrow(Beta)))
  cat(sprintf("After filtering: %d (%.1f%%)\n", nrow(Beta_filtered),
              (nrow(Beta_filtered) / nrow(Beta)) * 100))
  
  # QC density plot (only if not exists)
  qc_plot <- file.path(output_dir, "01_QC_DensityPlot.pdf")
  if (!file.exists(qc_plot)) {
    pdf(qc_plot, width = 10, height = 6)
    plot(density(Beta_filtered[, 1], na.rm = TRUE),
         main = "Beta Value Distribution",
         xlab = "Beta", ylab = "Density", col = "blue")
    for (i in 2:min(20, ncol(Beta_filtered))) {
      lines(density(Beta_filtered[, i], na.rm = TRUE), col = rgb(0, 0, 1, 0.2))
    }
    dev.off()
    cat("✓ QC density plot saved:", qc_plot, "\n")
  } else {
    cat("✓ QC density plot already exists, skipping.\n")
  }
} else {
  cat("\n--- Quality control skipped (no Beta matrix) ---\n")
}

# ----------------------------------------------------------------------------
# 5. DMP analysis: load existing or compute
# ----------------------------------------------------------------------------
cat("\n--- DMP analysis ---\n")

DMP_global <- NULL
dmp_computed <- FALSE

dmp_file <- file.path(input_data_dir, "DMP_Global_Results.csv")

if (has_beta && !is.null(Beta_filtered) && !is.null(targets)) {
  cat("Computing DMP analysis (limma)...\n")
  design_global <- model.matrix(~ targets$Phenotype)
  fit <- lmFit(Beta_filtered, design_global)
  fit <- eBayes(fit)
  DMP_global <- topTable(fit, coef = 2, number = Inf)
  
  sig_count <- sum(DMP_global$adj.P.Val < 0.05, na.rm = TRUE)
  hyper_count <- sum(DMP_global$adj.P.Val < 0.05 & DMP_global$logFC > 0, na.rm = TRUE)
  hypo_count <- sum(DMP_global$adj.P.Val < 0.05 & DMP_global$logFC < 0, na.rm = TRUE)
  
  cat(sprintf("✓ Significant DMPs (FDR<0.05): %d\n", sig_count))
  cat(sprintf("  - Hypermethylated: %d (%.1f%%)\n", hyper_count,
              ifelse(sig_count > 0, (hyper_count / sig_count) * 100, 0)))
  cat(sprintf("  - Hypomethylated: %d (%.1f%%)\n", hypo_count,
              ifelse(sig_count > 0, (hypo_count / sig_count) * 100, 0)))
  
  dmp_out <- file.path(output_dir, "01_DMP_Global_Results.csv")
  write.csv(DMP_global, dmp_out, row.names = TRUE)
  cat("✓ DMP results saved:", dmp_out, "\n")
  dmp_computed <- TRUE
  
} else if (file.exists(dmp_file)) {
  cat("Loading existing DMP results from:", dmp_file, "\n")
  DMP_global <- read.csv(dmp_file, row.names = 1)
  cat(sprintf("✓ Loaded DMP results: %d CpGs\n", nrow(DMP_global)))
  dmp_computed <- FALSE
} else {
  cat("✗ Cannot perform DMP analysis: missing Beta matrix or targets.\n")
}

# ----------------------------------------------------------------------------
# 6. DMP annotation: load existing or compute
# ----------------------------------------------------------------------------
cat("\n--- DMP annotation ---\n")

sig_dmp_annotated <- NULL
annot_computed <- FALSE

if (!is.null(DMP_global)) {
  cat("Annotating DMPs with CHAMP annotation...\n")
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  cpg_annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  common_cpg <- intersect(rownames(DMP_global), rownames(cpg_annot))
  DMP_global_sub <- DMP_global[common_cpg, , drop = FALSE]
  
  DMP_annotated <- cbind(
    DMP_global_sub,
    Gene = cpg_annot[common_cpg, "UCSC_RefGene_Name"],
    Chr  = cpg_annot[common_cpg, "chr"],
    Pos  = cpg_annot[common_cpg, "pos"]
  )
  
  sig_dmp_annotated <- DMP_annotated[DMP_annotated$adj.P.Val < 0.05, , drop = FALSE]
  
  annot_all_out <- file.path(output_dir, "02_DMP_Annotated_All.csv")
  annot_sig_out <- file.path(output_dir, "03_DMP_Annotated_Significant.csv")
  
  write.csv(DMP_annotated, annot_all_out, row.names = TRUE)
  write.csv(sig_dmp_annotated, annot_sig_out, row.names = TRUE)
  
  cat(sprintf("✓ %d significant DMPs annotated\n", nrow(sig_dmp_annotated)))
  cat("✓ All annotated DMPs saved:", annot_all_out, "\n")
  cat("✓ Significant annotated DMPs saved:", annot_sig_out, "\n")
  annot_computed <- TRUE
  
} else if (file.exists(file.path(input_data_dir, "03_DMP_Annotated_Significant.csv"))) {
  cat("Loading existing annotated DMPs...\n")
  sig_dmp_annotated <- read.csv(file.path(input_data_dir, "03_DMP_Annotated_Significant.csv"), row.names = 1)
  cat(sprintf("✓ Loaded %d annotated DMPs\n", nrow(sig_dmp_annotated)))
  annot_computed <- FALSE
} else {
  cat("✗ Cannot annotate DMPs: DMP results not available.\n")
}

# ----------------------------------------------------------------------------
# 7. Visualization: volcano and heatmap
# ----------------------------------------------------------------------------
cat("\n--- Visualization ---\n")

# Volcano plot
if (!is.null(DMP_global)) {
  volcano_out <- file.path(output_dir, "03_VolcanoPlot.pdf")
  pdf(volcano_out, width = 10, height = 8)
  plot(DMP_global$logFC, -log10(DMP_global$adj.P.Val),
       main = "Volcano Plot: MDD Suicide vs Control",
       xlab = "log2(Beta Change)", ylab = "-log10(Adjusted P-value)",
       col = ifelse(DMP_global$adj.P.Val < 0.05, "red", "gray50"),
       pch = 16, cex = 0.6)
  abline(h = -log10(0.05), lty = 2, col = "blue", lwd = 2)
  abline(v = c(-0.1, 0.1), lty = 2, col = "blue", lwd = 2)
  legend("topright", c("Significant (FDR<0.05)", "Not significant"),
         col = c("red", "gray50"), pch = 16)
  dev.off()
  cat("✓ Volcano plot saved:", volcano_out, "\n")
} else {
  cat("✗ Cannot create volcano plot: DMP results missing.\n")
}

# Heatmap (Top 30 DMPs)
if (!is.null(DMP_global) && has_beta && !is.null(Beta_filtered) && !is.null(targets)) {
  top_dmp_idx <- order(DMP_global$adj.P.Val)[1:min(30, nrow(DMP_global))]
  top_dmp_probes <- rownames(DMP_global)[top_dmp_idx]
  top_dmp_beta <- Beta_filtered[top_dmp_probes, , drop = FALSE]
  
  heatmap_out <- file.path(output_dir, "04_Top30_DMP_Heatmap.pdf")
  pdf(heatmap_out, width = 12, height = 10)
  pheatmap(top_dmp_beta,
           annotation_col = targets[, c("Region", "Phenotype"), drop = FALSE],
           main = "Heatmap: Top 30 DMPs",
           scale = "row", show_colnames = FALSE,
           clustering_distance_cols = "correlation")
  dev.off()
  cat("✓ Heatmap saved:", heatmap_out, "\n")
} else {
  cat("✗ Cannot create heatmap: missing Beta matrix or targets.\n")
}

# ----------------------------------------------------------------------------
# 8. DMR analysis (only if bumphunter available and not already done)
# ----------------------------------------------------------------------------
cat("\n--- DMR analysis ---\n")

dmr_file <- file.path(output_dir, "05_DMR_Results.csv")
if (file.exists(dmr_file)) {
  cat("✓ DMR results already exist:", dmr_file, "\n")
} else if (has_bumphunter && has_beta && !is.null(Beta_filtered) && !is.null(targets)) {
  cat("Running DMR analysis with bumphunter...\n")
  design <- model.matrix(~ targets$Phenotype)
  data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  annot <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  common_cpgs <- intersect(rownames(Beta_filtered), rownames(annot))
  Beta_sub <- Beta_filtered[common_cpgs, ]
  annot_sub <- annot[common_cpgs, ]
  pos <- data.frame(chr = annot_sub$chr, pos = annot_sub$pos)
  
  tryCatch({
    dmrs <- bumphunter::bumphunter(Beta_sub, design = design, pos = pos,
                                   cutoff = 0.2, B = 100, smooth = TRUE)
    dmr_table <- dmrs$table
    if (nrow(dmr_table) > 0) {
      write.csv(dmr_table, dmr_file, row.names = FALSE)
      cat(sprintf("✓ Found %d DMRs, saved to: %s\n", nrow(dmr_table), dmr_file))
    } else {
      cat("✓ No DMRs found.\n")
    }
  }, error = function(e) {
    cat("✗ DMR analysis failed:", e$message, "\n")
  })
} else {
  cat("✗ DMR analysis skipped: missing bumphunter, Beta matrix, or targets.\n")
}

# ----------------------------------------------------------------------------
# 9. GO/KEGG enrichment (only if missMethyl available)
# ----------------------------------------------------------------------------
cat("\n--- Pathway enrichment analysis ---\n")

go_file <- file.path(output_dir, "06_GO_Enrichment.csv")
kegg_file <- file.path(output_dir, "06_KEGG_Enrichment.csv")

if (file.exists(go_file) && file.exists(kegg_file)) {
  cat("✓ GO and KEGG enrichment results already exist.\n")
} else if (has_missMethyl && !is.null(sig_dmp_annotated) && !is.null(DMP_global)) {
  cat("Running GO/KEGG enrichment with missMethyl...\n")
  sig_cpgs <- gsub("^X", "", rownames(sig_dmp_annotated))
  all_cpgs <- gsub("^X", "", rownames(DMP_global))
  
  # GO enrichment
  if (!file.exists(go_file)) {
    tryCatch({
      go <- missMethyl::gometh(sig.cpg = sig_cpgs, all.cpg = all_cpgs, collection = "GO")
      write.csv(go, go_file, row.names = TRUE)
      cat(sprintf("✓ GO enrichment saved: %s (%d terms)\n", go_file,
                  sum(go$P.Value < 0.05, na.rm = TRUE)))
    }, error = function(e) {
      cat("✗ GO enrichment failed:", e$message, "\n")
    })
  }
  
  # KEGG enrichment
  if (!file.exists(kegg_file)) {
    tryCatch({
      kegg <- missMethyl::gometh(sig.cpg = sig_cpgs, all.cpg = all_cpgs, collection = "KEGG")
      write.csv(kegg, kegg_file, row.names = TRUE)
      cat(sprintf("✓ KEGG enrichment saved: %s (%d pathways)\n", kegg_file,
                  sum(kegg$P.Value < 0.05, na.rm = TRUE)))
    }, error = function(e) {
      cat("✗ KEGG enrichment failed:", e$message, "\n")
    })
  }
} else {
  cat("✗ Pathway enrichment skipped: missing missMethyl or significant DMPs.\n")
}

# ----------------------------------------------------------------------------
# 10. Final report
# ----------------------------------------------------------------------------
cat("\n--- Final report ---\n")

report_file <- file.path(output_dir, "FINAL_REPORT.txt")
if (file.exists(report_file)) {
  cat("✓ Final report already exists:", report_file, "\n")
} else {
  cat("Generating final report...\n")
  
  n_samples  <- ifelse(!is.null(targets), nrow(targets), 0)
  n_mdd     <- ifelse(!is.null(targets) && "Phenotype" %in% colnames(targets),
                      sum(targets$Phenotype == "MDD_suicide", na.rm = TRUE), 0)
  n_control <- ifelse(!is.null(targets) && "Phenotype" %in% colnames(targets),
                      sum(targets$Phenotype == "Control", na.rm = TRUE), 0)
  n_ba11    <- ifelse(!is.null(targets) && "Region" %in% colnames(targets),
                      sum(targets$Region == "BA11", na.rm = TRUE), 0)
  n_ba25    <- ifelse(!is.null(targets) && "Region" %in% colnames(targets),
                      sum(targets$Region == "BA25", na.rm = TRUE), 0)
  n_probes  <- ifelse(!is.null(Beta), nrow(Beta), 0)
  n_filtered <- ifelse(!is.null(Beta_filtered), nrow(Beta_filtered), 0)
  sig_count <- ifelse(!is.null(DMP_global),
                      sum(DMP_global$adj.P.Val < 0.05, na.rm = TRUE), 0)
  hyper_count <- ifelse(!is.null(DMP_global),
                        sum(DMP_global$adj.P.Val < 0.05 & DMP_global$logFC > 0, na.rm = TRUE), 0)
  hypo_count <- ifelse(!is.null(DMP_global),
                       sum(DMP_global$adj.P.Val < 0.05 & DMP_global$logFC < 0, na.rm = TRUE), 0)
  gene_count <- ifelse(!is.null(sig_dmp_annotated),
                       length(unique(na.omit(sig_dmp_annotated$Gene))), 0)
  
  final_report <- sprintf("
================================================================================
                  GSE88890 Analysis Report (Resume/Patch-Run)
         Brain Methylation and MDD Suicide Epigenetic Signatures
================================================================================

[Study Design]
Samples: %d (MDD suicide %d, Control %d)
Regions: BA11 (%d), BA25 (%d)
Platform: Illumina HumanMethylation450
Total probes: %d

[Quality Control]
Filtered probes: %d (%.1f%%)

[DMP Results]
Significant DMPs (FDR<0.05): %d
  - Hypermethylated: %d (%.1f%%)
  - Hypomethylated: %d (%.1f%%)

[Annotation]
Genes affected: %d

[Generated Files]
%s
%s
%s
%s
%s
%s

[Analysis Status]
- DMP: %s
- Annotation: %s
- DMR: %s
- Pathway enrichment: %s

[Note]
%s

================================================================================
Completed: %s
================================================================================
",
    n_samples, n_mdd, n_control,
    n_ba11, n_ba25,
    n_probes,
    n_filtered, ifelse(n_probes > 0, (n_filtered / n_probes) * 100, 0),
    sig_count,
    hyper_count, ifelse(sig_count > 0, (hyper_count / sig_count) * 100, 0),
    hypo_count, ifelse(sig_count > 0, (hypo_count / sig_count) * 100, 0),
    gene_count,
    ifelse(file.exists(file.path(output_dir, "01_QC_DensityPlot.pdf")), "✓ 01_QC_DensityPlot.pdf", "□ 01_QC_DensityPlot.pdf"),
    ifelse(file.exists(file.path(output_dir, "01_DMP_Global_Results.csv")), "✓ 01_DMP_Global_Results.csv", "□ 01_DMP_Global_Results.csv"),
    ifelse(file.exists(file.path(output_dir, "02_DMP_Annotated_Significant.csv")), "✓ 02_DMP_Annotated_Significant.csv", "□ 02_DMP_Annotated_Significant.csv"),
    ifelse(file.exists(file.path(output_dir, "03_VolcanoPlot.pdf")), "✓ 03_VolcanoPlot.pdf", "□ 03_VolcanoPlot.pdf"),
    ifelse(file.exists(file.path(output_dir, "04_Top30_DMP_Heatmap.pdf")), "✓ 04_Top30_DMP_Heatmap.pdf", "□ 04_Top30_DMP_Heatmap.pdf"),
    ifelse(file.exists(file.path(output_dir, "05_DMR_Results.csv")), "✓ 05_DMR_Results.csv", "□ 05_DMR_Results.csv"),
    ifelse(dmp_computed, "Newly computed", "Loaded from existing"),
    ifelse(annot_computed, "Newly computed", "Loaded from existing"),
    ifelse(file.exists(dmr_file), "Completed", "Not run / skipped"),
    ifelse(file.exists(go_file) || file.exists(kegg_file), "Completed", "Not run / skipped"),
    ifelse(!has_beta, "WARNING: Beta matrix not found - some analyses skipped.", "Normal"),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  
  writeLines(final_report, report_file)
  cat("✓ Final report saved:", report_file, "\n")
}

# ----------------------------------------------------------------------------
cat("\n✓✓✓ Resume Analysis Complete! ✓✓✓\n")
cat("Results saved to:", output_dir, "\n\n")