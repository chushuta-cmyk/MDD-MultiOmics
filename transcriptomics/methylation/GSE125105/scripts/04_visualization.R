#!/usr/bin/env Rscript
# ============================================================================
# 04_visualization.R
# GSE125105 甲基化数据可视化
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("04_visualization.R: GSE125105 Methylation Data Visualization\n")
cat("============================================================================\n")
cat("\n")

# 设置绝对路径
DATA_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/data"
RESULTS_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/results"
LOG_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/logs"

# 输出目录
VIS_DIR <- file.path(RESULTS_DIR, "04_visualization")
dir.create(VIS_DIR, showWarnings = FALSE, recursive = TRUE)

# 日志文件
log_file <- file.path(LOG_DIR, "04_visualization.log")
sink(log_file, append = FALSE, split = TRUE)

cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Visualization Output Directory:", VIS_DIR, "\n")
cat("\n")

# ============================================================================
# 1. 加载必要包
# ============================================================================
cat("1. Loading required packages...\n")

# 检查并安装必要包
required_packages <- c("ggplot2", "reshape2", "RColorBrewer", "pheatmap", "ggrepel", "gridExtra")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 加载数据操作包
if (!require("data.table", quietly = TRUE)) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
cat("✓ Packages loaded\n\n")

# ============================================================================
# 2. 加载数据
# ============================================================================
cat("2. Loading data from previous steps...\n")

# 2.1 加载样本信息
QC_DIR <- file.path(RESULTS_DIR, "01_qc")
targets <- readRDS(file.path(QC_DIR, "targets.rds"))
beta_info <- readRDS(file.path(QC_DIR, "beta_info.rds"))
beta_file <- beta_info$file_path

cat("   Sample information: ", nrow(targets), " samples\n")
cat("   Beta file: ", beta_file, "\n")

# 2.2 加载注释的DMP结果
ANNOT_DIR <- file.path(RESULTS_DIR, "03_annotation")
annotated_file <- file.path(ANNOT_DIR, "GSE125105_DMP_annotated_all.csv")

if (file.exists(annotated_file)) {
  dmp_annotated <- fread(annotated_file)
  cat("   Annotated DMP results: ", nrow(dmp_annotated), " probes\n")
} else {
  # 如果注释文件不存在，加载原始DMP结果
  DMP_DIR <- file.path(RESULTS_DIR, "02_dmp")
  dmp_file <- file.path(DMP_DIR, "GSE125105_DMP_all_results.csv")
  if (file.exists(dmp_file)) {
    dmp_annotated <- fread(dmp_file)
    cat("   DMP results (not annotated): ", nrow(dmp_annotated), " probes\n")
  } else {
    stop("No DMP results found for visualization")
  }
}

# 2.3 加载显著性DMP
sig_dmp_file <- file.path(ANNOT_DIR, "GSE125105_DMP_annotated_significant.csv")
if (file.exists(sig_dmp_file)) {
  sig_dmp <- fread(sig_dmp_file)
  cat("   Significant DMPs: ", nrow(sig_dmp), " probes\n")
} else {
  # 从完整结果中提取显著性DMP
  sig_dmp <- dmp_annotated[dmp_annotated$adj.P.Val < 0.05, ]
  cat("   Significant DMPs (extracted): ", nrow(sig_dmp), " probes\n")
}

cat("\n")

# ============================================================================
# 3. 高级火山图
# ============================================================================
cat("3. Generating enhanced volcano plot...\n")

# 准备火山图数据
volcano_data <- dmp_annotated
if (!"logFC" %in% colnames(volcano_data)) {
  stop("logFC column not found in DMP results")
}

# 添加显著性标签
volcano_data$Significance <- "Not significant"
volcano_data$Significance[volcano_data$adj.P.Val < 0.05 & abs(volcano_data$logFC) > 0.1] <- "Significant"
volcano_data$Significance[volcano_data$adj.P.Val < 0.05 & volcano_data$logFC > 0.1] <- "Hypermethylated"
volcano_data$Significance[volcano_data$adj.P.Val < 0.05 & volcano_data$logFC < -0.1] <- "Hypomethylated"

# 限制数据量用于绘图（如果数据太大）
if (nrow(volcano_data) > 20000) {
  set.seed(123)
  plot_indices <- sample(1:nrow(volcano_data), 20000)
  plot_data <- volcano_data[plot_indices, ]
  cat("   Sampled ", nrow(plot_data), " probes for volcano plot\n")
} else {
  plot_data <- volcano_data
}

# 颜色方案
color_scheme <- c("Not significant" = "gray70",
                  "Significant" = "orange",
                  "Hypermethylated" = "red",
                  "Hypomethylated" = "blue")

# 3.1 基础火山图
pdf(file.path(VIS_DIR, "volcano_plot_enhanced.pdf"), width = 12, height = 10)
p <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = color_scheme) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.7) +
  labs(title = "Enhanced Volcano Plot of Differential Methylation",
       subtitle = paste("GSE125105: MDD vs Control (", nrow(plot_data), "probes shown)"),
       x = expression(log[2]~Fold~Change~"(Case vs Control)"),
       y = expression(-log[10]~"(p-value)")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(p)
dev.off()
cat("   Saved: volcano_plot_enhanced.pdf\n")

# 3.2 标注top DMP的火山图
if (nrow(sig_dmp) > 0) {
  # 选择top 20个最显著的DMP进行标注
  top_n <- min(20, nrow(sig_dmp))
  top_dmp <- sig_dmp[order(sig_dmp$adj.P.Val)[1:top_n], ]
  
  # 添加基因名称（如果有）
  if ("UCSC_RefGene_Name" %in% colnames(top_dmp)) {
    top_dmp$Label <- top_dmp$UCSC_RefGene_Name
    # 取第一个基因（如果多个用分号分隔）
    top_dmp$Label <- sapply(strsplit(as.character(top_dmp$Label), ";"), function(x) x[1])
  } else {
    top_dmp$Label <- top_dmp$CpG
  }
  
  # 创建标注图
  pdf(file.path(VIS_DIR, "volcano_plot_labeled.pdf"), width = 14, height = 10)
  p <- ggplot(plot_data, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = color_scheme) +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50", alpha = 0.7) +
    # 标注top DMP
    geom_point(data = top_dmp, aes(x = logFC, y = -log10(P.Value)), 
               color = "black", size = 3, shape = 1) +
    geom_text_repel(data = top_dmp, 
                    aes(x = logFC, y = -log10(P.Value), label = Label),
                    size = 3, box.padding = 0.5, max.overlaps = Inf) +
    labs(title = paste("Volcano Plot with Top", top_n, "Significant DMPs Labeled"),
         subtitle = paste("GSE125105: MDD vs Control"),
         x = expression(log[2]~Fold~Change~"(Case vs Control)"),
         y = expression(-log[10]~"(p-value)")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))
  print(p)
  dev.off()
  cat("   Saved: volcano_plot_labeled.pdf\n")
}

cat("\n")

# ============================================================================
# 4. 热图：Top DMP的甲基化模式
# ============================================================================
cat("4. Generating heatmap of top DMPs...\n")

# 选择top DMP（最多50个）
top_n_heatmap <- min(50, nrow(sig_dmp))
if (top_n_heatmap > 0) {
  top_dmp_heatmap <- sig_dmp[order(sig_dmp$adj.P.Val)[1:top_n_heatmap], ]
  top_cpgs <- top_dmp_heatmap$CpG
  
  cat("   Selected ", length(top_cpgs), " top DMPs for heatmap\n")
  
  # 从Beta矩阵中提取这些CpG的数据
  cat("   Extracting Beta values for selected CpGs...\n")
  
  # 读取Beta文件并提取特定CpG
  con <- gzfile(beta_file, "r")
  header_line <- readLines(con, n = 1)
  header <- strsplit(header_line, "\t")[[1]]
  sample_ids <- header[-1]
  
  # 创建结果矩阵
  heatmap_data <- matrix(NA, nrow = length(top_cpgs), ncol = length(sample_ids))
  rownames(heatmap_data) <- top_cpgs
  colnames(heatmap_data) <- sample_ids
  
  # 查找并提取目标CpG
  found_cpgs <- character(0)
  chunk_size <- 10000
  chunk_num <- 1
  
  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) break
    
    chunk_data <- do.call(rbind, strsplit(lines, "\t"))
    chunk_cpgs <- chunk_data[, 1]
    
    # 检查是否有目标CpG在此块中
    target_indices <- which(chunk_cpgs %in% top_cpgs)
    if (length(target_indices) > 0) {
      for (idx in target_indices) {
        cpg <- chunk_cpgs[idx]
        beta_values <- as.numeric(chunk_data[idx, -1])
        heatmap_data[cpg, ] <- beta_values
        found_cpgs <- c(found_cpgs, cpg)
      }
    }
    
    # 如果已找到所有CpG，则停止
    if (length(found_cpgs) == length(top_cpgs)) {
      cat("      Found all ", length(top_cpgs), " CpGs in chunk ", chunk_num, "\n")
      break
    }
    
    chunk_num <- chunk_num + 1
    if (chunk_num > 100) {  # 安全限制
      cat("      Warning: Exceeded chunk limit, found ", length(found_cpgs), " of ", length(top_cpgs), " CpGs\n")
      break
    }
  }
  close(con)
  
  # 移除未找到的CpG
  heatmap_data <- heatmap_data[found_cpgs, ]
  
  if (length(found_cpgs) > 0) {
    cat("   Successfully extracted ", nrow(heatmap_data), " CpGs\n")
    
    # 准备注释信息
    annotation_col <- data.frame(
      Group = targets$Group[match(colnames(heatmap_data), targets$Sample_ID)],
      row.names = colnames(heatmap_data)
    )
    
    # 为行添加基因注释（如果有）
    if ("UCSC_RefGene_Name" %in% colnames(top_dmp_heatmap)) {
      gene_names <- top_dmp_heatmap$UCSC_RefGene_Name[match(rownames(heatmap_data), top_dmp_heatmap$CpG)]
      gene_names <- sapply(strsplit(as.character(gene_names), ";"), function(x) ifelse(is.na(x[1]), "", x[1]))
      row_labels <- paste0(rownames(heatmap_data), "\n", gene_names)
    } else {
      row_labels <- rownames(heatmap_data)
    }
    
    # 生成热图
    cat("   Generating heatmap...\n")
    pdf(file.path(VIS_DIR, "top_dmp_heatmap.pdf"), width = 14, height = 12)
    
    # 设置颜色
    heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
    
    pheatmap(heatmap_data,
             main = paste("Heatmap of Top", nrow(heatmap_data), "Differentially Methylated Probes"),
             annotation_col = annotation_col,
             annotation_colors = list(Group = c(control = "blue", case = "red")),
             color = heatmap_colors,
             scale = "row",  # 按行标准化
             clustering_method = "complete",
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             show_rownames = TRUE,
             show_colnames = FALSE,
             fontsize_row = 8,
             fontsize_col = 6,
             labels_row = row_labels)
    
    dev.off()
    cat("   Saved: top_dmp_heatmap.pdf\n")
    
    # 保存热图数据
    write.csv(heatmap_data, file.path(VIS_DIR, "top_dmp_heatmap_data.csv"))
    cat("   Saved: top_dmp_heatmap_data.csv\n")
    
  } else {
    cat("   ⚠️  No CpGs found in beta matrix, skipping heatmap\n")
  }
} else {
  cat("   ⚠️  No significant DMPs for heatmap\n")
}

cat("\n")

# ============================================================================
# 5. 曼哈顿图（如果染色体信息可用）
# ============================================================================
cat("5. Generating Manhattan plot...\n")

if (all(c("chr", "pos", "P.Value") %in% colnames(dmp_annotated))) {
  manhattan_data <- dmp_annotated[, c("CpG", "chr", "pos", "P.Value", "logFC", "adj.P.Val")]
  
  # 清理染色体名称
  manhattan_data$chr <- gsub("chr", "", manhattan_data$chr)
  manhattan_data$chr[manhattan_data$chr == "X"] <- "23"
  manhattan_data$chr[manhattan_data$chr == "Y"] <- "24"
  manhattan_data$chr[manhattan_data$chr == "M" | manhattan_data$chr == "MT"] <- "25"
  manhattan_data$chr <- as.numeric(manhattan_data$chr)
  manhattan_data <- manhattan_data[!is.na(manhattan_data$chr), ]
  
  # 计算染色体累计位置
  manhattan_data <- manhattan_data[order(manhattan_data$chr, manhattan_data$pos), ]
  
  # 添加累计位置
  data_cum <- manhattan_data %>% 
    group_by(chr) %>% 
    summarize(max_pos = max(pos)) %>% 
    mutate(add_pos = lag(cumsum(max_pos), default = 0)) %>% 
    select(chr, add_pos)
  
  manhattan_data <- manhattan_data %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(pos_cum = pos + add_pos)
  
  # 准备轴标签
  axis_df <- manhattan_data %>% 
    group_by(chr) %>% 
    summarize(center = mean(pos_cum))
  
  # 染色体颜色
  manhattan_data$color <- ifelse(manhattan_data$chr %% 2 == 1, "black", "gray50")
  
  # 显著性阈值线
  sig_threshold <- -log10(0.05 / nrow(manhattan_data))  # Bonferroni校正
  sug_threshold <- -log10(1e-5)  # 建议阈值
  
  # 创建曼哈顿图
  pdf(file.path(VIS_DIR, "manhattan_plot.pdf"), width = 16, height = 6)
  p <- ggplot(manhattan_data, aes(x = pos_cum, y = -log10(P.Value), color = color)) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("black" = "black", "gray50" = "gray50")) +
    scale_x_continuous(label = c(1:22, "X", "Y", "M"), breaks = axis_df$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(manhattan_data$P.Value)) * 1.1)) +
    geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = sug_threshold, linetype = "dashed", color = "blue", alpha = 0.5) +
    labs(title = "Manhattan Plot of Methylation Association",
         x = "Chromosome",
         y = expression(-log[10]~"(p-value)")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  dev.off()
  cat("   Saved: manhattan_plot.pdf\n")
  
  # 按染色体着色的曼哈顿图
  pdf(file.path(VIS_DIR, "manhattan_plot_colored.pdf"), width = 16, height = 6)
  p <- ggplot(manhattan_data, aes(x = pos_cum, y = -log10(P.Value), color = as.factor(chr %% 2))) +
    geom_point(alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("0" = "deepskyblue", "1" = "orange")) +
    scale_x_continuous(label = c(1:22, "X", "Y", "M"), breaks = axis_df$center) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(-log10(manhattan_data$P.Value)) * 1.1)) +
    geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "red", alpha = 0.5) +
    geom_hline(yintercept = sug_threshold, linetype = "dashed", color = "blue", alpha = 0.5) +
    labs(title = "Manhattan Plot of Methylation Association (Colored by Chromosome)",
         x = "Chromosome",
         y = expression(-log[10]~"(p-value)")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  print(p)
  dev.off()
  cat("   Saved: manhattan_plot_colored.pdf\n")
  
} else {
  cat("   ⚠️  Chromosome or position information not available, skipping Manhattan plot\n")
}

cat("\n")

# ============================================================================
# 6. 其他可视化
# ============================================================================
cat("6. Generating additional visualizations...\n")

# 6.1 logFC分布图
pdf(file.path(VIS_DIR, "logfc_distribution.pdf"), width = 10, height = 6)
hist(dmp_annotated$logFC, breaks = 100, 
     main = "Distribution of log2 Fold Changes",
     xlab = "log2 Fold Change (Case vs Control)",
     ylab = "Frequency",
     col = "lightblue", border = "white")
abline(v = c(-0.1, 0, 0.1), col = c("blue", "black", "red"), lwd = 2, lty = c(2, 1, 2))
legend("topright", legend = c("Hypomethylated threshold", "No change", "Hypermethylated threshold"),
       col = c("blue", "black", "red"), lwd = 2, lty = c(2, 1, 2))
dev.off()
cat("   Saved: logfc_distribution.pdf\n")

# 6.2 -log10(p-value) vs logFC散点图（抽样）
if (nrow(dmp_annotated) > 10000) {
  set.seed(123)
  scatter_data <- dmp_annotated[sample(1:nrow(dmp_annotated), 10000), ]
} else {
  scatter_data <- dmp_annotated
}

pdf(file.path(VIS_DIR, "pvalue_vs_logfc_scatter.pdf"), width = 10, height = 8)
ggplot(scatter_data, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(alpha = 0.3, color = "darkgray") +
  geom_point(data = scatter_data[scatter_data$adj.P.Val < 0.05 & scatter_data$logFC > 0.1, ],
             aes(x = logFC, y = -log10(P.Value)), color = "red", alpha = 0.5) +
  geom_point(data = scatter_data[scatter_data$adj.P.Val < 0.05 & scatter_data$logFC < -0.1, ],
             aes(x = logFC, y = -log10(P.Value)), color = "blue", alpha = 0.5) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  labs(title = "Scatter Plot: -log10(p-value) vs logFC",
       x = "log2 Fold Change (Case vs Control)",
       y = "-log10(p-value)") +
  theme_minimal()
dev.off()
cat("   Saved: pvalue_vs_logfc_scatter.pdf\n")

# 6.3 显著性DMP的logFC分布（小提琴图）
if (nrow(sig_dmp) > 0) {
  sig_dmp$Direction <- ifelse(sig_dmp$logFC > 0, "Hypermethylated", "Hypomethylated")
  
  pdf(file.path(VIS_DIR, "sig_dmp_logfc_violin.pdf"), width = 8, height = 6)
  ggplot(sig_dmp, aes(x = Direction, y = logFC, fill = Direction)) +
    geom_violin(alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = c("Hypermethylated" = "red", "Hypomethylated" = "blue")) +
    labs(title = "Distribution of logFC in Significant DMPs",
         subtitle = paste("Total", nrow(sig_dmp), "significant DMPs (FDR < 0.05)"),
         x = "Methylation Direction",
         y = "log2 Fold Change") +
    theme_minimal() +
    theme(legend.position = "none")
  dev.off()
  cat("   Saved: sig_dmp_logfc_violin.pdf\n")
}

# 6.4 按染色体/基因区域统计的条形图（如果有注释）
if (all(c("chr", "adj.P.Val") %in% colnames(dmp_annotated))) {
  # 按染色体统计显著性DMP
  sig_by_chr <- dmp_annotated %>%
    filter(adj.P.Val < 0.05) %>%
    group_by(chr) %>%
    summarize(Count = n()) %>%
    arrange(desc(Count))
  
  if (nrow(sig_by_chr) > 0) {
    pdf(file.path(VIS_DIR, "sig_dmp_by_chromosome.pdf"), width = 12, height = 6)
    ggplot(sig_by_chr, aes(x = reorder(chr, Count), y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(title = "Significant DMPs by Chromosome",
           subtitle = "FDR < 0.05",
           x = "Chromosome",
           y = "Number of Significant DMPs") +
      theme_minimal()
    dev.off()
    cat("   Saved: sig_dmp_by_chromosome.pdf\n")
  }
}

cat("\n")

# ============================================================================
# 7. 创建可视化摘要
# ============================================================================
cat("7. Creating visualization summary...\n")

# 生成HTML报告（如果可用）
if (require("rmarkdown", quietly = TRUE)) {
  # 创建简单的R Markdown报告
  report_text <- paste0(
    "---\n",
    "title: \"GSE125105 Methylation Visualization Summary\"\n",
    "output: html_document\n",
    "---\n\n",
    "## Visualization Summary\n\n",
    "Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
    "### Data Overview\n",
    "- Total probes analyzed: ", nrow(dmp_annotated), "\n",
    "- Significant DMPs (FDR < 0.05): ", sum(dmp_annotated$adj.P.Val < 0.05, na.rm = TRUE), "\n",
    "- Hypermethylated DMPs: ", sum(dmp_annotated$adj.P.Val < 0.05 & dmp_annotated$logFC > 0.1, na.rm = TRUE), "\n",
    "- Hypomethylated DMPs: ", sum(dmp_annotated$adj.P.Val < 0.05 & dmp_annotated$logFC < -0.1, na.rm = TRUE), "\n\n",
    "### Generated Visualizations\n",
    "1. Enhanced Volcano Plot\n",
    "2. Labeled Volcano Plot (top DMPs)\n",
    "3. Top DMP Heatmap\n",
    "4. Manhattan Plot\n",
    "5. Distribution Plots\n\n",
    "### Files\n",
    "All visualization files are saved in: ", VIS_DIR, "\n"
  )
  
  writeLines(report_text, file.path(VIS_DIR, "visualization_summary.Rmd"))
  
  tryCatch({
    rmarkdown::render(file.path(VIS_DIR, "visualization_summary.Rmd"),
                      output_file = file.path(VIS_DIR, "visualization_summary.html"),
                      quiet = TRUE)
    cat("   Saved: visualization_summary.html\n")
  }, error = function(e) {
    cat("   ⚠️  Could not generate HTML report: ", e$message, "\n")
  })
}

# 保存会话信息
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(VIS_DIR, "visualization_session_info.txt"))
cat("   Saved: visualization_session_info.txt\n")

# ============================================================================
# 8. 完成
# ============================================================================
cat("\n")
cat("Visualization Summary:\n")
cat("- Total probes analyzed:", nrow(dmp_annotated), "\n")
cat("- Significant DMPs (FDR < 0.05):", sum(dmp_annotated$adj.P.Val < 0.05, na.rm = TRUE), "\n")
cat("- Generated visualizations:", length(list.files(VIS_DIR, pattern = "\\.pdf$")), "PDF files\n")
cat("- Output directory:", VIS_DIR, "\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n✓ Visualization step completed successfully!\n")

sink()

# 返回成功状态
quit(status = 0)