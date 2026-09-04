# ============================================================================
# FUTURE WORK: Integrated Follow-Up Analysis Pipeline
# Based on: Immune Exhaustion and Ciliary Dysfunction Hypothesis
# Assumes: expr_data, sample_metadata, hub_genes, etc. are loaded from file
# ============================================================================

cat("\n=== Follow-Up Analysis Pipeline Started ===\n")

# ----------------------------------------------------------------------------
# 0. Setup: load required packages and data
# ----------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(igraph)
library(STRINGdb)
library(immunedeconv)
library(clusterProfiler)
library(org.Hs.eg.db)
library(mediation)    # optional

# Load pre-computed objects (must exist)
data_file <- "GSE54563_analysis_objects.RData"
if (file.exists(data_file)) {
  load(data_file)
  cat("Loaded expression and metadata from", data_file, "\n")
} else {
  stop("Required RData file not found. Run GSE54563 pipeline first.")
}

# If hub genes file exists, load it
hub_file <- "WGCNA_Hub_Genes_Final_Complete.csv"
if (file.exists(hub_file)) {
  hub_genes <- read.csv(hub_file, stringsAsFactors = FALSE)
} else {
  cat("WARNING: Hub genes file not found. Skipping ciliary hub overlap.\n")
  hub_genes <- NULL
}

# ----------------------------------------------------------------------------
# PRIORITY 1: Validate B cell findings in blood (GSE39653)
# ----------------------------------------------------------------------------
cat("\n=== PRIORITY 1: Blood B Cell Validation ===\n")

# 1.1 Fetch or load GSE39653
if (!exists("blood_expr")) {
  tryCatch({
    gse <- getGEO("GSE39653", GSEMatrix = TRUE)[[1]]
    blood_expr <- exprs(gse)
    blood_meta <- pData(gse)
    save(blood_expr, blood_meta, file = "GSE39653_blood_data.RData")
  }, error = function(e) {
    cat("Falling back to local file GSE39653_blood_data.RData\n")
    load("GSE39653_blood_data.RData")
  })
}

# 1.2 Prepare expression matrix (log2 if needed)
blood_mat <- as.matrix(blood_expr)
if (min(blood_mat) > 0 & quantile(blood_mat, 0.99) > 100) {
  blood_mat <- log2(blood_mat + 1)
}

# 1.3 MCP-counter deconvolution
blood_immune <- deconvolute(blood_mat, method = "mcp_counter")
blood_scores <- as.matrix(blood_immune[, -1])
rownames(blood_scores) <- blood_immune$cell_type

# 1.4 Extract B cell row (flexible matching)
bcell_row <- grep("B.cell|B cell", rownames(blood_scores), ignore.case = TRUE)[1]
if (!is.na(bcell_row)) {
  bcell <- as.numeric(blood_scores[bcell_row, ])
  # Extract diagnosis from title (adjust pattern)
  diag <- ifelse(grepl("^MDD", blood_meta$title), "MDD", 
                 ifelse(grepl("^HC", blood_meta$title), "Control", NA))
  keep <- !is.na(diag)
  bcell <- bcell[keep]
  diag <- diag[keep]
  w <- wilcox.test(bcell[diag == "MDD"], bcell[diag == "Control"])
  fc <- mean(bcell[diag == "MDD"]) / mean(bcell[diag == "Control"])
  cat(sprintf("Blood B cell: FC=%.3f, p=%.4f\n", fc, w$p.value))
  # Save result
  write.csv(data.frame(FC=fc, Pvalue=w$p.value, N_MDD=sum(diag=="MDD"), N_Control=sum(diag=="Control")),
            "Blood_BCell_Result.csv", row.names = FALSE)
}

# ----------------------------------------------------------------------------
# PRIORITY 2: Ciliary gene functional analysis
# ----------------------------------------------------------------------------
cat("\n=== PRIORITY 2: Ciliary Gene Analysis ===\n")

ciliary_genes <- c("CFAP126","CFAP43","CFAP300","SPAG6","DNAI4","DYNLRB2",
                   "CFAP45","CFAP53","RSPH1","RSPH3","CCDC39","CCDC40",
                   "LRRC6","DRC1","GAS8","HEATR2","HYDIN","WDR16")

# 2.1 Check overlap with hub genes
if (!is.null(hub_genes)) {
  overlap <- hub_genes[hub_genes$GeneSymbol %in% ciliary_genes, ]
  write.csv(overlap, "Ciliary_Overlap_Hub.csv", row.names = FALSE)
  cat("Ciliary genes in hub set:", nrow(overlap), "\n")
}

# 2.2 STRING PPI (with fallback)
ppi <- NULL
tryCatch({
  string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400)
  mapped <- string_db$map(data.frame(genes=ciliary_genes), "genes", removeUnmappedRows=TRUE)
  ppi <- string_db$get_interactions(mapped$STRING_id)
  write.csv(ppi, "Ciliary_PPI.csv", row.names = FALSE)
  # Network metrics
  g <- graph_from_data_frame(ppi[, c("from","to")], directed=FALSE)
  metrics <- data.frame(Gene=names(degree(g)), Degree=degree(g), Betweenness=betweenness(g))
  write.csv(metrics, "Ciliary_Network_Metrics.csv", row.names = FALSE)
}, error = function(e) cat("STRING PPI skipped: ", e$message, "\n"))

# 2.3 GO enrichment
entrez <- bitr(ciliary_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
go <- enrichGO(gene=entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05)
if (nrow(go@result) > 0) {
  write.csv(go@result, "Ciliary_GO.csv", row.names = FALSE)
  pdf("Ciliary_GO_DotPlot.pdf")
  print(dotplot(go, showCategory=20))
  dev.off()
}

# ----------------------------------------------------------------------------
# PRIORITY 3: AQP4-Immune correlation
# ----------------------------------------------------------------------------
cat("\n=== PRIORITY 3: AQP4 - B Cell Correlation ===\n")

# AQP4 expression (assuming it's in expr_data)
if ("AQP4" %in% rownames(expr_data)) {
  aqp4 <- as.numeric(expr_data["AQP4", ])
  # Load B cell scores from brain (pre-computed, e.g., from MCP-counter on brain data)
  # Assume bcell_scores is available; otherwise skip
  if (exists("bcell_scores") && length(aqp4) == length(bcell_scores)) {
    corr <- cor.test(aqp4, bcell_scores, method="spearman")
    cat("AQP4-B cell correlation: Rho=", round(corr$estimate,3), " p=", round(corr$p.value,4), "\n")
    # Region-stratified
    if (exists("sample_metadata")) {
      regions <- unique(sample_metadata$Brain_Region)
      region_res <- data.frame()
      for (r in regions) {
        idx <- sample_metadata$Brain_Region == r
        if (sum(idx) > 3) {
          cr <- cor.test(aqp4[idx], bcell_scores[idx], method="spearman")
          region_res <- rbind(region_res, data.frame(Region=r, Rho=cr$estimate, P=cr$p.value))
        }
      }
      write.csv(region_res, "AQP4_BCell_by_Region.csv", row.names = FALSE)
    }
    # Mediation (optional)
    if (require(mediation, quietly=TRUE)) {
      med_df <- data.frame(AQP4=scale(aqp4), BCells=scale(bcell_scores),
                           Dx=ifelse(sample_metadata$Diagnosis=="MDD",1,0))
      med <- mediate(lm(BCells ~ AQP4, med_df), 
                     glm(Dx ~ AQP4 + BCells, family=binomial, med_df),
                     treat="AQP4", mediator="BCells", boot=TRUE, sims=500)
      write.csv(summary(med)$boot.ci, "Mediation_Results.csv")
    }
  }
}

# ----------------------------------------------------------------------------
# PRIORITY 4: Immune exhaustion markers
# ----------------------------------------------------------------------------
cat("\n=== PRIORITY 4: Exhaustion Markers ===\n")

exhaustion_markers <- c("PDCD1","LAG3","HAVCR2","TIGIT","CTLA4","TOX","TOX2","BATF","NR4A1","ENTPD1")
exhaustion_expr <- expr_data[rownames(expr_data) %in% exhaustion_markers, ]
if (nrow(exhaustion_expr) > 0) {
  exhaustion_expr <- exhaustion_expr[match(exhaustion_markers, rownames(exhaustion_expr)), ]
  exhaustion_score <- colMeans(exhaustion_expr, na.rm=TRUE)
  # Compare MDD vs Control (assuming Diagnosis in sample_metadata)
  if (exists("sample_metadata")) {
    diag <- sample_metadata$Diagnosis
    mdd <- exhaustion_score[diag == "MDD"]
    ctrl <- exhaustion_score[diag == "Control"]
    w <- wilcox.test(mdd, ctrl)
    fc <- mean(mdd)/mean(ctrl)
    cat("Exhaustion score: FC=", round(fc,3), " p=", round(w$p.value,4), "\n")
    # Per marker
    marker_res <- data.frame()
    for (g in exhaustion_markers) {
      if (g %in% rownames(expr_data)) {
        e <- as.numeric(expr_data[g, ])
        m <- e[diag=="MDD"]; c <- e[diag=="Control"]
        wt <- wilcox.test(m, c)
        marker_res <- rbind(marker_res, data.frame(Gene=g, FC=mean(m)/mean(c), P=wt$p.value))
      }
    }
    marker_res$FDR <- p.adjust(marker_res$P, method="BH")
    write.csv(marker_res, "Exhaustion_Markers.csv", row.names = FALSE)
  }
}

# ----------------------------------------------------------------------------
# PRIORITY 5: Single-cell deconvolution (framework only)
# ----------------------------------------------------------------------------
cat("\n=== PRIORITY 5: Single-Cell Deconvolution (Outline) ===\n")
cat("This requires a single-cell reference matrix. Implement when data available.\n")
cat("Suggested code:\n")
cat("  sc_res <- deconvolute(expr_data, method='quantiseq', reference=sc_ref)\n")
cat("  # Then extract B cell subtypes and compare as above.\n\n")

cat("\n=== All analyses completed ===\n")
