#!/usr/bin/env Rscript
# ============================================================================
# 01_load_and_qc_fast.R
# GSE125105 甲基化数据快速加载与QC
# 快速版本：不读取整个文件，仅进行基本检查
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("01_load_and_qc_fast.R: GSE125105 Methylation Data Fast Loading and QC\n")
cat("============================================================================\n")
cat("\n")

# 设置绝对路径
DATA_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/data"
RESULTS_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/results"
LOG_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/logs"

# 输出目录
QC_DIR <- file.path(RESULTS_DIR, "01_qc")
dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)

# 日志文件
log_file <- file.path(LOG_DIR, "01_load_and_qc.log")
sink(log_file, append = FALSE, split = TRUE)

cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Data Directory:", DATA_DIR, "\n")
cat("Results Directory:", RESULTS_DIR, "\n")
cat("QC Output Directory:", QC_DIR, "\n")
cat("\n")

# ============================================================================
# 1. 加载必要包
# ============================================================================
cat("1. Loading required packages...\n")
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}
if (!require("reshape2", quietly = TRUE)) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}
cat("✓ Packages loaded\n\n")

# ============================================================================
# 2. 加载样本信息
# ============================================================================
cat("2. Loading sample information...\n")
targets_file <- file.path(DATA_DIR, "GSE125105_CHAMP_targets.csv")
if (!file.exists(targets_file)) {
  stop("Targets file not found: ", targets_file)
}
targets <- fread(targets_file)
colnames(targets) <- c("Sample_ID", "Group")
targets$Group <- factor(targets$Group, levels = c("control", "case"))
cat("   Samples in targets:", nrow(targets), "\n")
cat("   Groups:", paste(levels(targets$Group), collapse = ", "), "\n")
cat("   Group counts:\n")
print(table(targets$Group))
cat("\n")

# ============================================================================
# 3. 从系列矩阵中提取样本映射
# ============================================================================
cat("3. Extracting sample mapping from series matrix...\n")
series_file <- file.path(DATA_DIR, "GSE125105_series_matrix.txt.gz")
if (!file.exists(series_file)) {
  stop("Series matrix file not found: ", series_file)
}

series_lines <- readLines(series_file)
title_line <- grep("^!Sample_title", series_lines, value = TRUE)
title_parts <- strsplit(title_line, "\t")[[1]]
titles <- gsub('^"|"$', '', title_parts[-1])

gsm_line <- grep("^!Sample_geo_accession", series_lines, value = TRUE)
gsm_parts <- strsplit(gsm_line, "\t")[[1]]
gsm_ids <- gsub('^"|"$', '', gsm_parts[-1])

sample_mapping <- data.frame(
  Sample_Name = titles,
  GSM_ID = gsm_ids,
  stringsAsFactors = FALSE
)
sample_mapping$Sample_Label <- gsub("genomic DNA from ", "", sample_mapping$Sample_Name)
cat("   Sample mapping created:", nrow(sample_mapping), "samples\n")
cat("   First few mappings:\n")
print(head(sample_mapping))
cat("\n")

# ============================================================================
# 4. 检查Beta矩阵文件
# ============================================================================
cat("4. Checking Beta matrix file...\n")
beta_file <- file.path(DATA_DIR, "GSE125105_matrix_normalized.txt.gz")
if (!file.exists(beta_file)) {
  stop("Beta matrix file not found: ", beta_file)
}

# 获取文件大小和行数
cat("   Getting file info...\n")
file_size <- file.info(beta_file)$size
cat("   File size:", round(file_size / 1024^3, 2), "GB\n")

# 使用系统命令获取行数（更快）
system_cmd <- paste("zcat", beta_file, "| wc -l")
line_count <- as.numeric(system(system_cmd, intern = TRUE))
cat("   Total lines (probes + header):", line_count, "\n")
total_probes <- line_count - 1
cat("   Estimated total probes:", total_probes, "\n")

# 读取标题行
cat("   Reading header...\n")
con <- gzfile(beta_file, "r")
header_line <- readLines(con, n = 1)
close(con)

header_cols <- strsplit(header_line, "\t")[[1]]
cat("   Total columns in file:", length(header_cols), "\n")

# 识别Beta值列（不包含"DetectionPval"的列）
beta_cols <- header_cols[!grepl("DetectionPval", header_cols)]
beta_cols <- beta_cols[beta_cols == "ID_REF" | !grepl("DetectionPval", beta_cols)]
sample_cols <- beta_cols[-1]
cat("   Beta columns (without DetectionPval):", length(beta_cols), "\n")
cat("   Sample columns for Beta values:", length(sample_cols), "\n")

# 创建从样本列到GSM ID的映射
col_to_gsm <- setNames(sample_mapping$GSM_ID, sample_mapping$Sample_Label)

# 重新排列targets以匹配矩阵列顺序
available_samples <- intersect(sample_cols, sample_mapping$Sample_Label)
cat("   Samples available in matrix:", length(available_samples), "\n")

targets_ordered <- data.frame(
  Sample_Label = available_samples,
  GSM_ID = col_to_gsm[available_samples],
  stringsAsFactors = FALSE
)
targets_ordered$Group <- targets$Group[match(targets_ordered$GSM_ID, targets$Sample_ID)]
targets_ordered <- targets_ordered[!is.na(targets_ordered$Group), ]
cat("   Samples with group information:", nrow(targets_ordered), "\n")
cat("   Final group counts:\n")
print(table(targets_ordered$Group))
cat("\n")

# ============================================================================
# 5. 快速QC：读取前1000行进行分析
# ============================================================================
cat("5. Performing quick QC on first 1000 probes...\n")
n_probes_qc <- 1000
cat("   Reading first", n_probes_qc, "probes...\n")

con <- gzfile(beta_file, "r")
header <- readLines(con, n = 1)
lines <- readLines(con, n = n_probes_qc)
close(con)

if (length(lines) > 0) {
  # 解析数据
  chunk_data <- do.call(rbind, strsplit(lines, "\t"))
  cpg_ids <- chunk_data[, 1]
  
  # 提取Beta值列
  beta_col_indices <- which(header_cols %in% beta_cols)
  beta_qc <- matrix(as.numeric(chunk_data[, beta_col_indices[-1]]), 
                    nrow = length(cpg_ids), 
                    ncol = length(beta_col_indices) - 1)
  colnames(beta_qc) <- sample_cols
  rownames(beta_qc) <- cpg_ids
  
  # 基本统计
  cat("   Beta value range in sample:", round(range(beta_qc, na.rm = TRUE), 4), "\n")
  na_count <- sum(is.na(beta_qc))
  cat("   NA count in sample:", na_count, "\n")
  cat("   NA percentage in sample:", round(na_count / length(beta_qc) * 100, 4), "%\n")
  
  # 5.1 Beta值分布密度图
  cat("   5.1 Beta value density plot...\n")
  pdf(file.path(QC_DIR, "beta_density_plot.pdf"), width = 10, height = 6)
  plot(density(as.vector(beta_qc), na.rm = TRUE),
       main = paste("Distribution of Beta Values (First", n_probes_qc, "probes)"),
       xlab = "Beta Value",
       ylab = "Density",
       col = "blue",
       lwd = 2)
  abline(v = c(0, 0.2, 0.8, 1), col = "gray", lty = 2)
  legend("topright", legend = c("Beta values"), col = "blue", lwd = 2)
  dev.off()
  cat("      Saved: beta_density_plot.pdf\n")
  
  # 5.2 按分组的Beta值分布箱线图
  cat("   5.2 Group-wise beta distribution...\n")
  beta_long <- melt(beta_qc)
  colnames(beta_long) <- c("Probe", "Sample_Label", "Beta")
  beta_long$GSM_ID <- col_to_gsm[as.character(beta_long$Sample_Label)]
  beta_long$Group <- targets$Group[match(beta_long$GSM_ID, targets$Sample_ID)]
  
  pdf(file.path(QC_DIR, "beta_boxplot_by_group.pdf"), width = 12, height = 6)
  p <- ggplot(beta_long, aes(x = Group, y = Beta, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    labs(title = paste("Beta Value Distribution by Group (First", n_probes_qc, "probes)"),
         x = "Group", y = "Beta Value") +
    theme_minimal() +
    theme(legend.position = "none")
  print(p)
  dev.off()
  cat("      Saved: beta_boxplot_by_group.pdf\n")
  
  # 5.3 MDS图
  cat("   5.3 MDS plot...\n")
  # 转置并移除有缺失值的列
  beta_t <- t(beta_qc)
  beta_t <- beta_t[, apply(beta_t, 2, function(x) !any(is.na(x)))]
  if (ncol(beta_t) > 0) {
    # 仅保留有分组信息的样本
    valid_samples <- intersect(rownames(beta_t), targets_ordered$Sample_Label)
    beta_t <- beta_t[valid_samples, ]
    
    if (nrow(beta_t) > 1) {
      dist_matrix <- dist(beta_t)
      mds_result <- cmdscale(dist_matrix, k = 2)
      
      mds_groups <- targets_ordered$Group[match(rownames(beta_t), targets_ordered$Sample_Label)]
      
      pdf(file.path(QC_DIR, "mds_plot.pdf"), width = 10, height = 8)
      plot(mds_result, 
           col = as.numeric(mds_groups),
           pch = 16,
           cex = 1.5,
           main = paste("MDS Plot of Samples (First", n_probes_qc, "probes)"),
           xlab = "MDS Dimension 1",
           ylab = "MDS Dimension 2")
      legend("topright", 
             legend = levels(mds_groups),
             col = 1:nlevels(mds_groups),
             pch = 16,
             title = "Group")
      dev.off()
      cat("      Saved: mds_plot.pdf\n")
    }
  }
  
} else {
  cat("   ⚠️  No lines read from Beta matrix\n")
}

cat("\n")

# ============================================================================
# 6. 保存QC摘要
# ============================================================================
cat("6. Saving QC summary...\n")

qc_summary <- data.frame(
  Metric = c("Total samples", "Total probes (estimated)", "Group: control", "Group: case",
             "Beta value min (sample)", "Beta value max (sample)", "NA percentage (sample)"),
  Value = c(nrow(targets_ordered),
            total_probes,
            sum(targets_ordered$Group == "control"),
            sum(targets_ordered$Group == "case"),
            if (exists("beta_qc")) round(min(beta_qc, na.rm = TRUE), 4) else "NA",
            if (exists("beta_qc")) round(max(beta_qc, na.rm = TRUE), 4) else "NA",
            if (exists("beta_qc")) paste0(round(na_count / length(beta_qc) * 100, 4), "%") else "NA")
)

write.csv(qc_summary, file.path(QC_DIR, "qc_summary.csv"), row.names = FALSE)
cat("      Saved: qc_summary.csv\n")

# 保存targets（已验证顺序）
write.csv(targets_ordered, file.path(QC_DIR, "targets_verified.csv"), row.names = FALSE)
cat("      Saved: targets_verified.csv\n")

# ============================================================================
# 7. 保存中间数据供后续步骤使用
# ============================================================================
cat("7. Saving intermediate data for downstream steps...\n")

# 保存targets为RDS
saveRDS(targets_ordered, file.path(QC_DIR, "targets.rds"))
cat("      Saved: targets.rds\n")

# 保存样本顺序
sample_order <- targets_ordered$Sample_Label
saveRDS(sample_order, file.path(QC_DIR, "sample_order.rds"))
cat("      Saved: sample_order.rds\n")

# 保存样本映射
saveRDS(col_to_gsm, file.path(QC_DIR, "sample_mapping.rds"))
cat("      Saved: sample_mapping.rds\n")

# 保存Beta矩阵文件路径和列信息
beta_info <- list(
  file_path = beta_file,
  total_probes = total_probes,
  sample_count = nrow(targets_ordered),
  beta_cols = beta_cols,
  sample_cols = sample_cols
)
saveRDS(beta_info, file.path(QC_DIR, "beta_info.rds"))
cat("      Saved: beta_info.rds\n")

# ============================================================================
# 8. 完成
# ============================================================================
cat("\n")
cat("QC Summary:\n")
print(qc_summary)
cat("\n")

cat("Output files saved to:", QC_DIR, "\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n✓ QC step completed successfully!\n")

sink()

# 返回成功状态
quit(status = 0)