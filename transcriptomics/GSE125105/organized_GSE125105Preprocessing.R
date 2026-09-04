# =============================================================================
# GSE125105 DMP Analysis (Limma, no CHAMP)
# Chunk-based processing for large Beta matrix
# =============================================================================

cat("\n=== GSE125105 DMP Analysis Started ===\n")

library(data.table)
library(limma)

# -----------------------------------------------------------------------------
# 1. Load phenotype from series matrix (metadata lines)
# -----------------------------------------------------------------------------
meta_file <- "/data03/karama/projects/MDD_analysis/GSE125105/data/GSE125105_series_matrix.txt.gz"
meta_lines <- readLines(meta_file)

# Extract sample IDs and phenotype
sample_line <- grep("^!Series_sample_id", meta_lines, value = TRUE)
sample_ids <- gsub("^!Series_sample_id\\t\"|\"$", "", sample_line)
sample_ids <- unlist(strsplit(sample_ids, " "))

pheno_line <- grep("^!Sample_source_name_ch1", meta_lines, value = TRUE)
phenos <- gsub("^!Sample_source_name_ch1\\t|\"$", "", pheno_line)
phenos <- unlist(strsplit(phenos, "\"\t\""))
phenos <- gsub("_whole blood|\\s+", "", phenos)
phenos <- ifelse(grepl("^control", phenos, ignore.case = TRUE), "control", phenos)
phenos[sample_ids == "GSM3562834"] <- "control"   # manual fix for one sample
pheno <- factor(phenos, levels = c("control", "case"))

cat(sprintf("Phenotype: %d cases, %d controls\n", sum(pheno == "case"), sum(pheno == "control")))

# -----------------------------------------------------------------------------
# 2. Design matrix
# -----------------------------------------------------------------------------
design <- model.matrix(~pheno)

# -----------------------------------------------------------------------------
# 3. Process Beta matrix in chunks
# -----------------------------------------------------------------------------
# Note: the actual Beta file path should be corrected (the original had a duplication)
beta_file <- "/data03/karama/projects/MDD_analysis/GSE125105/data/GSE125105_matrix_normalized.txt.gz"
con <- gzfile(beta_file, "r")
header <- readLines(con, n = 1)

chunk_size <- 50000
results_list <- list()
chunk_idx <- 1

repeat {
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0) break
  
  # Write to temp file for fread (avoids parsing issues)
  tmp <- sprintf("/tmp/beta_chunk_%d.txt", chunk_idx)
  writeLines(c(header, lines), tmp)
  
  chunk <- fread(tmp, sep = "\t", header = TRUE, data.table = TRUE)
  cpg_ids <- chunk[[1]]
  chunk[, 1 := NULL]                    # remove first column (CpG names)
  
  # Convert to numeric matrix
  for (j in seq_along(chunk)) set(chunk, j = j, value = as.numeric(chunk[[j]]))
  beta_mat <- as.matrix(chunk)
  rownames(beta_mat) <- cpg_ids
  
  # Limma fit
  fit <- lmFit(beta_mat, design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = 2, number = Inf)
  results_list[[chunk_idx]] <- res
  
  # Clean up
  rm(chunk, beta_mat, fit, res)
  gc()
  file.remove(tmp)
  
  chunk_idx <- chunk_idx + 1
}

close(con)

cat(sprintf("Processed %d chunks.\n", chunk_idx - 1))

# -----------------------------------------------------------------------------
# 4. Combine and rank results
# -----------------------------------------------------------------------------
dmp_all <- do.call(rbind, results_list)
dmp_all <- dmp_all[order(dmp_all$adj.P.Val), ]

sig <- sum(dmp_all$adj.P.Val < 0.05)
hyper <- sum(dmp_all$adj.P.Val < 0.05 & dmp_all$logFC > 0)
hypo  <- sum(dmp_all$adj.P.Val < 0.05 & dmp_all$logFC < 0)

cat(sprintf("Significant DMPs (FDR<0.05): %d (hyper=%d, hypo=%d)\n", sig, hyper, hypo))

# -----------------------------------------------------------------------------
# 5. Save results
# -----------------------------------------------------------------------------
write.csv(dmp_all, "GSE125105_DMP_Results.csv", row.names = TRUE)
cat("Results saved to GSE125105_DMP_Results.csv\n\n")
