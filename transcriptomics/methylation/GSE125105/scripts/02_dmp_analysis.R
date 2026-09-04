#!/usr/bin/env Rscript
# ============================================================================
# 02_dmp_analysis.R
# GSE125105 差异甲基化探针分析（分块limma） - 修复版
# 修复：处理DetectionPval列，使用正确的样本标签
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("02_dmp_analysis.R: GSE125105 Differential Methylation Analysis (Chunked Limma)\n")
cat("============================================================================\n")
cat("\n")

# 设置绝对路径
DATA_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/data"
RESULTS_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/results"
LOG_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/logs"

# 输出目录
DMP_DIR <- file.path(RESULTS_DIR, "02_dmp")
dir.create(DMP_DIR, showWarnings = FALSE, recursive = TRUE)

# 日志文件
log_file <- file.path(LOG_DIR, "02_dmp_analysis.log")
sink(log_file, append = FALSE, split = TRUE)

cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results Directory:", DMP_DIR, "\n")
cat("\n")

# ============================================================================
# 1. 加载必要包
# ============================================================================
cat("1. Loading required packages...\n")
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
if (!require("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("limma")
  library(limma)
}
cat("✓ Packages loaded\n\n")

# ============================================================================
# 2. 加载QC步骤的输出
# ============================================================================
cat("2. Loading QC step outputs...\n")
QC_DIR <- file.path(RESULTS_DIR, "01_qc")

targets <- readRDS(file.path(QC_DIR, "targets.rds"))
sample_order <- readRDS(file.path(QC_DIR, "sample_order.rds"))
beta_info <- readRDS(file.path(QC_DIR, "beta_info.rds"))
sample_mapping <- readRDS(file.path(QC_DIR, "sample_mapping.rds"))

beta_file <- beta_info$file_path
total_probes <- beta_info$total_probes
sample_count <- beta_info$sample_count
beta_cols <- beta_info$beta_cols
beta_col_indices <- beta_info$beta_col_indices

cat("   Beta file:", beta_file, "\n")
cat("   Total probes:", total_probes, "\n")
cat("   Sample count:", sample_count, "\n")
cat("   Beta columns:", length(beta_cols), "\n")
cat("   Group distribution:\n")
print(table(targets$Group))
cat("\n")

# ============================================================================
# 3. 准备设计矩阵
# ============================================================================
cat("3. Preparing design matrix...\n")
# 确保因子顺序：control为参考
targets$Group <- factor(targets$Group, levels = c("control", "case"))
design <- model.matrix(~ Group, data = targets)
rownames(design) <- targets$Sample_Label
colnames(design) <- c("Intercept", "Case_vs_Control")

cat("   Design matrix dimensions:", dim(design), "\n")
cat("   Coefficients:", colnames(design), "\n")
cat("   Design matrix head:\n")
print(head(design))
cat("\n")

# ============================================================================
# 4. 分块进行limma分析
# ============================================================================
cat("4. Performing chunked limma analysis...\n")

chunk_size <- 10000  # 每个块的探针数
all_results <- list()
chunk_num <- 1
processed_probes <- 0

# 打开文件连接
con <- gzfile(beta_file, "r")
header_line <- readLines(con, n = 1)  # 读取标题行
header <- strsplit(header_line, "\t")[[1]]

# 样本标签（来自beta列，排除ID_REF）
sample_labels <- beta_cols[-1]
cat("   Sample labels in matrix:", length(sample_labels), "\n")

# 验证样本顺序匹配
if (!identical(sample_labels, sample_order)) {
  cat("   WARNING: Sample order in file does not match saved order!\n")
  cat("   Reordering targets and design matrix to match file...\n")
  
  # 重新排序targets和设计矩阵以匹配文件顺序
  targets <- targets[match(sample_labels, targets$Sample_Label), ]
  design <- model.matrix(~ Group, data = targets)
  rownames(design) <- targets$Sample_Label
  colnames(design) <- c("Intercept", "Case_vs_Control")
}

cat("   Starting chunked processing...\n")
cat("   Chunk size:", chunk_size, "probes\n")

repeat {
  cat(sprintf("   Chunk %d: ", chunk_num))
  
  # 读取块
  lines <- readLines(con, n = chunk_size)
  if (length(lines) == 0) break
  
  # 解析为矩阵（仅Beta值列）
  chunk_data <- do.call(rbind, strsplit(lines, "\t"))
  
  # 检查列数是否匹配标题（可能有行号列）
  ncols <- ncol(chunk_data)
  header_len <- length(header)
  
  if (ncols == header_len) {
    # 没有行号列，第一列是CpG ID
    cpg_ids <- chunk_data[, 1]
    use_indices <- beta_col_indices
  } else if (ncols == header_len + 1) {
    # 有行号列，第一列是行号，第二列是CpG ID
    cpg_ids <- chunk_data[, 2]
    # 调整索引：所有列向右移动1位
    use_indices <- beta_col_indices + 1
    cat("Row number column detected, adjusting indices... ")
  } else {
    cat(sprintf("ERROR: chunk has %d columns, but header has %d columns. Cannot proceed.\n", ncols, header_len))
    stop("Column count mismatch")
  }
  
  # 提取Beta值列（使用调整后的索引）
  beta_chunk <- matrix(as.numeric(chunk_data[, use_indices[-1]]), 
                       nrow = length(cpg_ids), 
                       ncol = length(use_indices) - 1)
  rownames(beta_chunk) <- cpg_ids
  colnames(beta_chunk) <- sample_labels
  
  # 检查并处理缺失值
  na_count <- sum(is.na(beta_chunk))
  if (na_count > 0) {
    cat(sprintf("%d probes, %d NA values... ", length(cpg_ids), na_count))
    # 移除全部为NA的探针
    all_na <- apply(is.na(beta_chunk), 1, all)
    if (any(all_na)) {
      beta_chunk <- beta_chunk[!all_na, ]
      cpg_ids <- cpg_ids[!all_na]
      cat(sprintf("removed %d all-NA probes, ", sum(all_na)))
    }
  } else {
    cat(sprintf("%d probes... ", length(cpg_ids)))
  }
  
  # 如果块中还有探针，运行limma
  if (nrow(beta_chunk) > 0) {
    # 确保样本顺序一致
    if (!identical(colnames(beta_chunk), rownames(design))) {
      beta_chunk <- beta_chunk[, rownames(design)]
    }
    
    # 执行limma分析
    fit <- lmFit(beta_chunk, design)
    fit <- eBayes(fit)
    
    # 提取结果
    chunk_results <- topTable(fit, coef = "Case_vs_Control", number = Inf)
    
    # 确保有logFC列
    if (!"logFC" %in% colnames(chunk_results)) {
      # 尝试从列名中推断
      possible_fc_col <- grep("Groupcase|Case", colnames(chunk_results), value = TRUE)
      if (length(possible_fc_col) > 0) {
        chunk_results$logFC <- chunk_results[[possible_fc_col[1]]]
      }
    }
    
    # 添加CpG ID作为列
    chunk_results$CpG <- rownames(chunk_results)
    
    # 存储结果
    all_results[[chunk_num]] <- chunk_results
    
    processed_probes <- processed_probes + nrow(beta_chunk)
    cat(sprintf("analyzed %d probes (total: %d)\n", nrow(beta_chunk), processed_probes))
  } else {
    cat("no probes to analyze\n")
  }
  
  # 清理内存
  rm(chunk_data, beta_chunk, fit, chunk_results)
  gc()
  
  chunk_num <- chunk_num + 1
}

close(con)

cat("\n")
cat("   Total chunks processed:", chunk_num - 1, "\n")
cat("   Total probes analyzed:", processed_probes, "\n")
cat("   Percentage of total probes:", round(processed_probes / total_probes * 100, 2), "%\n")
cat("\n")

# ============================================================================
# 5. 合并所有结果
# ============================================================================
cat("5. Combining results from all chunks...\n")

if (length(all_results) == 0) {
  stop("No results generated from any chunk!")
}

# 合并所有数据框
DMP_results <- do.call(rbind, all_results)

# 确保列顺序一致
expected_cols <- c("CpG", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
available_cols <- colnames(DMP_results)

# 重新排序列
DMP_results <- DMP_results[, intersect(expected_cols, available_cols)]

# 按调整后p值排序
DMP_results <- DMP_results[order(DMP_results$adj.P.Val), ]

cat("   Combined results:", nrow(DMP_results), "probes\n")
cat("   Column names:", paste(colnames(DMP_results), collapse = ", "), "\n")
cat("\n")

# ============================================================================
# 6. 统计显著性
# ============================================================================
cat("6. Calculating significance statistics...\n")

# 定义显著性阈值
fdr_threshold <- 0.05
logfc_threshold <- 0.1  # 较小的logFC阈值，因为甲基化差异通常较小

# 统计
sig_dmp <- sum(DMP_results$adj.P.Val < fdr_threshold, na.rm = TRUE)
sig_dmp_logfc <- sum(DMP_results$adj.P.Val < fdr_threshold & 
                      abs(DMP_results$logFC) > logfc_threshold, na.rm = TRUE)
hyper_dmp <- sum(DMP_results$adj.P.Val < fdr_threshold & 
                  DMP_results$logFC > logfc_threshold, na.rm = TRUE)
hypo_dmp <- sum(DMP_results$adj.P.Val < fdr_threshold & 
                 DMP_results$logFC < -logfc_threshold, na.rm = TRUE)

cat("   Significant DMPs (FDR <", fdr_threshold, "):", sig_dmp, "\n")
cat("   Significant DMPs with |logFC| >", logfc_threshold, ":", sig_dmp_logfc, "\n")
cat("   Hypermethylated (logFC >", logfc_threshold, "):", hyper_dmp, "\n")
cat("   Hypomethylated (logFC < -", logfc_threshold, "):", hypo_dmp, "\n")
cat("\n")

# ============================================================================
# 7. 保存结果
# ============================================================================
cat("7. Saving results...\n")

# 7.1 完整DMP结果
dmp_file_all <- file.path(DMP_DIR, "GSE125105_DMP_all_results.csv")
write.csv(DMP_results, dmp_file_all, row.names = FALSE)
cat("   All DMP results:", dmp_file_all, "\n")

# 7.2 显著性DMP结果（FDR < 0.05）
sig_dmp_results <- DMP_results[DMP_results$adj.P.Val < fdr_threshold, ]
dmp_file_sig <- file.path(DMP_DIR, "GSE125105_DMP_significant.csv")
write.csv(sig_dmp_results, dmp_file_sig, row.names = FALSE)
cat("   Significant DMP results:", dmp_file_sig, "\n")

# 7.3 显著性且|logFC| > 0.1的DMP结果
sig_logfc_results <- DMP_results[DMP_results$adj.P.Val < fdr_threshold & 
                                  abs(DMP_results$logFC) > logfc_threshold, ]
dmp_file_sig_logfc <- file.path(DMP_DIR, "GSE125105_DMP_significant_logfc.csv")
write.csv(sig_logfc_results, dmp_file_sig_logfc, row.names = FALSE)
cat("   Significant DMP with |logFC| > 0.1:", dmp_file_sig_logfc, "\n")

# 7.4 Top 1000 DMP（按p值排序）
top_dmp <- head(DMP_results, 1000)
dmp_file_top <- file.path(DMP_DIR, "GSE125105_DMP_top1000.csv")
write.csv(top_dmp, dmp_file_top, row.names = FALSE)
cat("   Top 1000 DMP:", dmp_file_top, "\n")

# 7.5 统计摘要
stats_summary <- data.frame(
  Metric = c("Total probes analyzed", 
             "FDR threshold", 
             "logFC threshold",
             "Significant DMPs (FDR < 0.05)",
             "Significant DMPs (FDR < 0.05 & |logFC| > 0.1)",
             "Hypermethylated DMPs",
             "Hypomethylated DMPs",
             "Top significant CpG",
             "Top hypermethylated CpG",
             "Top hypomethylated CpG"),
  Value = c(
    processed_probes,
    fdr_threshold,
    logfc_threshold,
    sig_dmp,
    sig_dmp_logfc,
    hyper_dmp,
    hypo_dmp,
    ifelse(nrow(DMP_results) > 0, DMP_results$CpG[1], "NA"),
    ifelse(hyper_dmp > 0, DMP_results$CpG[which(DMP_results$adj.P.Val < fdr_threshold & 
                                                 DMP_results$logFC > logfc_threshold)[1]], "NA"),
    ifelse(hypo_dmp > 0, DMP_results$CpG[which(DMP_results$adj.P.Val < fdr_threshold & 
                                                DMP_results$logFC < -logfc_threshold)[1]], "NA")
  )
)

write.csv(stats_summary, file.path(DMP_DIR, "dmp_statistics_summary.csv"), row.names = FALSE)
cat("   Statistics summary:", file.path(DMP_DIR, "dmp_statistics_summary.csv"), "\n")

# 7.6 保存RDS格式供下游使用
saveRDS(DMP_results, file.path(DMP_DIR, "dmp_results.rds"))
cat("   RDS format:", file.path(DMP_DIR, "dmp_results.rds"), "\n")

# ============================================================================
# 8. 生成简单可视化
# ============================================================================
cat("8. Generating basic visualizations...\n")

if (nrow(DMP_results) > 0 && require("ggplot2", quietly = TRUE)) {
  # 8.1 Volcano plot（抽样，因为数据量太大）
  cat("   8.1 Volcano plot (sampled)...\n")
  
  # 随机抽样最多10000个点用于火山图
  set.seed(123)
  if (nrow(DMP_results) > 10000) {
    plot_data <- DMP_results[sample(1:nrow(DMP_results), 10000), ]
  } else {
    plot_data <- DMP_results
  }
  
  # 添加显著性分类
  plot_data$Significance <- "Not significant"
  plot_data$Significance[plot_data$adj.P.Val < fdr_threshold & 
                          abs(plot_data$logFC) > logfc_threshold] <- "Significant"
  plot_data$Significance[plot_data$adj.P.Val < fdr_threshold & 
                          plot_data$logFC > logfc_threshold] <- "Hypermethylated"
  plot_data$Significance[plot_data$adj.P.Val < fdr_threshold & 
                          plot_data$logFC < -logfc_threshold] <- "Hypomethylated"
  
  pdf(file.path(DMP_DIR, "volcano_plot.pdf"), width = 10, height = 8)
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Not significant" = "gray70",
                                  "Significant" = "orange",
                                  "Hypermethylated" = "red",
                                  "Hypomethylated" = "blue")) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = -log10(fdr_threshold), 
               linetype = "dashed", color = "gray50") +
    labs(title = "Volcano Plot of Differential Methylation",
         subtitle = paste("FDR <", fdr_threshold, ", |logFC| >", logfc_threshold),
         x = "log2 Fold Change (Case vs Control)",
         y = "-log10(p-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  print(p)
  dev.off()
  cat("      Saved: volcano_plot.pdf\n")
  
  # 8.2 P-value分布直方图
  cat("   8.2 P-value distribution histogram...\n")
  pdf(file.path(DMP_DIR, "pvalue_distribution.pdf"), width = 10, height = 6)
  hist(DMP_results$P.Value, breaks = 50, 
       main = "Distribution of P-values",
       xlab = "P-value", ylab = "Frequency",
       col = "skyblue", border = "white")
  abline(v = fdr_threshold, col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("FDR =", fdr_threshold), 
         col = "red", lwd = 2, lty = 2)
  dev.off()
  cat("      Saved: pvalue_distribution.pdf\n")
}

# ============================================================================
# 9. 完成
# ============================================================================
cat("\n")
cat("DMP Analysis Summary:\n")
print(stats_summary)
cat("\n")

cat("Output files saved to:", DMP_DIR, "\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n✓ DMP analysis completed successfully!\n")

sink()

# 返回成功状态
quit(status = 0)