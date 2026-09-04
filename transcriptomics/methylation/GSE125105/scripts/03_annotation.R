#!/usr/bin/env Rscript
# ============================================================================
# 03_annotation.R
# GSE125105 DMP结果注释
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("03_annotation.R: GSE125105 DMP Results Annotation\n")
cat("============================================================================\n")
cat("\n")

# 设置绝对路径
RESULTS_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/results"
LOG_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/logs"

# 输出目录
ANNOT_DIR <- file.path(RESULTS_DIR, "03_annotation")
dir.create(ANNOT_DIR, showWarnings = FALSE, recursive = TRUE)

# 日志文件
log_file <- file.path(LOG_DIR, "03_annotation.log")
sink(log_file, append = FALSE, split = TRUE)

cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Annotation Output Directory:", ANNOT_DIR, "\n")
cat("\n")

# ============================================================================
# 1. 加载必要包
# ============================================================================
cat("1. Loading required packages...\n")

# 尝试加载甲基化注释包
annotation_loaded <- FALSE
annotation_method <- ""

if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
  cat("   IlluminaHumanMethylation450kanno.ilmn12.hg19 not available.\n")
  cat("   Trying alternative: IlluminaHumanMethylation450kanno.ilmn12.hg38...\n")
  
  if (!require("IlluminaHumanMethylation450kanno.ilmn12.hg38", quietly = TRUE)) {
    cat("   IlluminaHumanMethylation450kanno.ilmn12.hg38 not available.\n")
    cat("   Trying to install from Bioconductor...\n")
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    tryCatch({
      BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
      library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
      annotation_loaded <- TRUE
      annotation_method <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    }, error = function(e) {
      cat("   Failed to install hg19 annotation. Trying hg38...\n")
      tryCatch({
        BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg38")
        library(IlluminaHumanMethylation450kanno.ilmn12.hg38)
        annotation_loaded <- TRUE
        annotation_method <- "IlluminaHumanMethylation450kanno.ilmn12.hg38"
      }, error = function(e2) {
        cat("   Failed to install annotation packages.\n")
        cat("   Will use minimal annotation from built-in data.\n")
      })
    })
  } else {
    library(IlluminaHumanMethylation450kanno.ilmn12.hg38)
    annotation_loaded <- TRUE
    annotation_method <- "IlluminaHumanMethylation450kanno.ilmn12.hg38"
  }
} else {
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  annotation_loaded <- TRUE
  annotation_method <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
}

if (annotation_loaded) {
  cat("   ✓ Annotation package loaded:", annotation_method, "\n")
} else {
  cat("   ⚠️  Using fallback annotation method\n")
}

# 加载其他必要包
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
# 2. 加载DMP结果
# ============================================================================
cat("2. Loading DMP results...\n")
DMP_DIR <- file.path(RESULTS_DIR, "02_dmp")

# 尝试加载完整结果或显著性结果
dmp_file_all <- file.path(DMP_DIR, "GSE125105_DMP_all_results.csv")
dmp_file_sig <- file.path(DMP_DIR, "GSE125105_DMP_significant.csv")

if (file.exists(dmp_file_all)) {
  dmp_results <- fread(dmp_file_all)
  cat("   Loaded all DMP results:", nrow(dmp_results), "probes\n")
} else if (file.exists(dmp_file_sig)) {
  dmp_results <- fread(dmp_file_sig)
  cat("   Loaded significant DMP results:", nrow(dmp_results), "probes\n")
} else {
  stop("No DMP results found in ", DMP_DIR)
}

# 检查必要的列
required_cols <- c("CpG", "logFC", "P.Value", "adj.P.Val")
missing_cols <- setdiff(required_cols, colnames(dmp_results))
if (length(missing_cols) > 0) {
  stop("Missing required columns in DMP results: ", paste(missing_cols, collapse = ", "))
}

cat("   Columns:", paste(colnames(dmp_results), collapse = ", "), "\n")
cat("   First few CpGs:", paste(head(dmp_results$CpG, 5), collapse = ", "), "\n")
cat("\n")

# ============================================================================
# 3. 注释DMP结果
# ============================================================================
cat("3. Annotating DMP results...\n")

if (annotation_loaded) {
  cat("   Using", annotation_method, "for annotation\n")
  
  # 获取注释数据
  if (annotation_method == "IlluminaHumanMethylation450kanno.ilmn12.hg19") {
    data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else if (annotation_method == "IlluminaHumanMethylation450kanno.ilmn12.hg38") {
    data("IlluminaHumanMethylation450kanno.ilmn12.hg38")
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg38)
  }
  
  # 检查注释数据的结构
  cat("   Annotation dimensions:", dim(anno), "\n")
  cat("   Annotation columns:", paste(colnames(anno)[1:10], collapse=", "), "...\n")
  
  # 选择有用的列
  useful_cols <- c("chr", "pos", "strand", "Name", "UCSC_RefGene_Name", 
                   "UCSC_RefGene_Group", "UCSC_RefGene_Accession",
                   "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island",
                   "Phantom", "DMR", "Enhancer", "HMM_Island", "Regulatory_Feature_Name",
                   "Regulatory_Feature_Group", "DHS", "Methyl450_Loci")
  
  # 保留实际存在的列
  available_cols <- intersect(useful_cols, colnames(anno))
  cat("   Available annotation columns:", length(available_cols), "\n")
  
  # 提取探针的注释信息
  probe_annotation <- anno[rownames(anno) %in% dmp_results$CpG, available_cols, drop = FALSE]
  probe_annotation$CpG <- rownames(probe_annotation)
  
  cat("   Found annotation for", nrow(probe_annotation), "probes\n")
  
  # 合并注释与DMP结果
  annotated_results <- merge(dmp_results, probe_annotation, by = "CpG", all.x = TRUE)
  
} else {
  cat("   ⚠️  Annotation package not available. Using minimal annotation...\n")
  
  # 创建基本注释（仅CpG位置信息）
  # 这里我们解析CpG ID来获取基本位置信息（如果ID格式为cg00000029）
  # 实际上，如果没有注释包，我们只能提供有限的注释
  
  annotated_results <- dmp_results
  
  # 添加占位符列
  annotated_results$chr <- NA
  annotated_results$pos <- NA
  annotated_results$strand <- NA
  annotated_results$UCSC_RefGene_Name <- NA
  annotated_results$UCSC_RefGene_Group <- NA
  annotated_results$Relation_to_UCSC_CpG_Island <- NA
  
  cat("   ⚠️  Only basic placeholder annotation added\n")
}

cat("   Annotated results:", nrow(annotated_results), "rows\n")
cat("   Added columns:", setdiff(colnames(annotated_results), colnames(dmp_results)), "\n")
cat("\n")

# ============================================================================
# 4. 增强注释：添加自定义分类
# ============================================================================
cat("4. Enhancing annotation with custom categories...\n")

# 4.1 按染色体统计
if ("chr" %in% colnames(annotated_results)) {
  chr_stats <- table(annotated_results$chr, useNA = "ifany")
  cat("   Probes per chromosome:\n")
  print(chr_stats)
  
  # 添加染色体分组
  annotated_results$Chromosome_Group <- NA
  annotated_results$Chromosome_Group[annotated_results$chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5")] <- "Large autosomes"
  annotated_results$Chromosome_Group[annotated_results$chr %in% c("chr6", "chr7", "chr8", "chr9", "chr10")] <- "Medium autosomes"
  annotated_results$Chromosome_Group[annotated_results$chr %in% c("chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")] <- "Small autosomes"
  annotated_results$Chromosome_Group[annotated_results$chr %in% c("chrX", "chrY")] <- "Sex chromosomes"
  annotated_results$Chromosome_Group[annotated_results$chr == "chrM"] <- "Mitochondrial"
} else {
  cat("   ⚠️  Chromosome information not available\n")
}

# 4.2 按CpG岛关系分类
if ("Relation_to_UCSC_CpG_Island" %in% colnames(annotated_results)) {
  island_stats <- table(annotated_results$Relation_to_UCSC_CpG_Island, useNA = "ifany")
  cat("   CpG Island relations:\n")
  print(island_stats)
} else {
  cat("   ⚠️  CpG Island relation information not available\n")
}

# 4.3 按基因区域分类
if ("UCSC_RefGene_Group" %in% colnames(annotated_results)) {
  # 简化基因区域分类
  annotated_results$Gene_Region_Simple <- NA
  
  # 解析基因区域（可能有多个用分号分隔）
  gene_groups <- strsplit(as.character(annotated_results$UCSC_RefGene_Group), ";")
  
  # 定义区域优先级
  region_priority <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR", "ExonBnd", "IGR")
  
  for (i in 1:length(gene_groups)) {
    regions <- unique(gene_groups[[i]])
    regions <- regions[!is.na(regions) & regions != ""]
    
    if (length(regions) > 0) {
      # 按优先级选择第一个区域
      for (priority_region in region_priority) {
        if (priority_region %in% regions) {
          annotated_results$Gene_Region_Simple[i] <- priority_region
          break
        }
      }
      # 如果没有匹配的优先级区域，使用第一个区域
      if (is.na(annotated_results$Gene_Region_Simple[i])) {
        annotated_results$Gene_Region_Simple[i] <- regions[1]
      }
    }
  }
  
  gene_region_stats <- table(annotated_results$Gene_Region_Simple, useNA = "ifany")
  cat("   Gene region distribution:\n")
  print(gene_region_stats)
} else {
  cat("   ⚠️  Gene region information not available\n")
}

# 4.4 添加显著性分类
cat("   Adding significance classification...\n")
annotated_results$Significance <- "Not significant"
annotated_results$Significance[annotated_results$adj.P.Val < 0.05] <- "Significant (FDR<0.05)"
annotated_results$Significance[annotated_results$adj.P.Val < 0.05 & 
                                abs(annotated_results$logFC) > 0.1] <- "Significant (FDR<0.05, |logFC|>0.1)"
annotated_results$Significance[annotated_results$adj.P.Val < 0.05 & 
                                annotated_results$logFC > 0.1] <- "Hypermethylated"
annotated_results$Significance[annotated_results$adj.P.Val < 0.05 & 
                                annotated_results$logFC < -0.1] <- "Hypomethylated"

sig_stats <- table(annotated_results$Significance, useNA = "ifany")
cat("   Significance classification:\n")
print(sig_stats)
cat("\n")

# ============================================================================
# 5. 保存注释结果
# ============================================================================
cat("5. Saving annotated results...\n")

# 5.1 完整注释结果
annot_file_all <- file.path(ANNOT_DIR, "GSE125105_DMP_annotated_all.csv")
write.csv(annotated_results, annot_file_all, row.names = FALSE)
cat("   All annotated results:", annot_file_all, "\n")

# 5.2 显著性DMP的注释结果
annotated_sig <- annotated_results[annotated_results$adj.P.Val < 0.05, ]
annot_file_sig <- file.path(ANNOT_DIR, "GSE125105_DMP_annotated_significant.csv")
write.csv(annotated_sig, annot_file_sig, row.names = FALSE)
cat("   Annotated significant DMPs:", annot_file_sig, "\n")

# 5.3 按染色体分割的文件（如果染色体信息可用）
if ("chr" %in% colnames(annotated_results)) {
  chr_list <- unique(annotated_results$chr)
  chr_list <- chr_list[!is.na(chr_list)]
  
  for (chr_name in chr_list) {
    chr_data <- annotated_results[annotated_results$chr == chr_name, ]
    if (nrow(chr_data) > 0) {
      chr_file <- file.path(ANNOT_DIR, paste0("GSE125105_DMP_chr", gsub("chr", "", chr_name), ".csv"))
      write.csv(chr_data, chr_file, row.names = FALSE)
    }
  }
  cat("   Split by chromosome: ", length(chr_list), " files\n")
}

# 5.4 按基因区域分割的文件
if ("Gene_Region_Simple" %in% colnames(annotated_results)) {
  regions <- unique(annotated_results$Gene_Region_Simple)
  regions <- regions[!is.na(regions)]
  
  for (region in regions) {
    region_data <- annotated_results[annotated_results$Gene_Region_Simple == region, ]
    if (nrow(region_data) > 0) {
      region_file <- file.path(ANNOT_DIR, paste0("GSE125105_DMP_region_", gsub("[^A-Za-z0-9]", "_", region), ".csv"))
      write.csv(region_data, region_file, row.names = FALSE)
    }
  }
  cat("   Split by gene region: ", length(regions), " files\n")
}

# 5.5 保存RDS格式
saveRDS(annotated_results, file.path(ANNOT_DIR, "annotated_dmp_results.rds"))
cat("   RDS format:", file.path(ANNOT_DIR, "annotated_dmp_results.rds"), "\n")

# ============================================================================
# 6. 生成注释摘要报告
# ============================================================================
cat("6. Generating annotation summary report...\n")

# 创建摘要数据框
summary_list <- list()

# 基本统计
summary_list[["Total_probes"]] <- nrow(annotated_results)
summary_list[["Significant_probes_FDR_0.05"]] <- sum(annotated_results$adj.P.Val < 0.05)
summary_list[["Significant_probes_FDR_0.01"]] <- sum(annotated_results$adj.P.Val < 0.01)
summary_list[["Hypermethylated_probes"]] <- sum(annotated_results$adj.P.Val < 0.05 & annotated_results$logFC > 0.1)
summary_list[["Hypomethylated_probes"]] <- sum(annotated_results$adj.P.Val < 0.05 & annotated_results$logFC < -0.1)

# 染色体统计（如果有）
if ("chr" %in% colnames(annotated_results)) {
  chr_counts <- table(annotated_results$chr[annotated_results$adj.P.Val < 0.05], useNA = "ifany")
  for (chr_name in names(chr_counts)) {
    summary_list[[paste0("Significant_probes_", chr_name)]] <- chr_counts[chr_name]
  }
}

# 基因区域统计（如果有）
if ("Gene_Region_Simple" %in% colnames(annotated_results)) {
  region_counts <- table(annotated_results$Gene_Region_Simple[annotated_results$adj.P.Val < 0.05], useNA = "ifany")
  for (region_name in names(region_counts)) {
    if (!is.na(region_name)) {
      summary_list[[paste0("Significant_probes_region_", gsub("[^A-Za-z0-9]", "_", region_name))]] <- region_counts[region_name]
    }
  }
}

# CpG岛关系统计（如果有）
if ("Relation_to_UCSC_CpG_Island" %in% colnames(annotated_results)) {
  island_counts <- table(annotated_results$Relation_to_UCSC_CpG_Island[annotated_results$adj.P.Val < 0.05], useNA = "ifany")
  for (island_name in names(island_counts)) {
    if (!is.na(island_name)) {
      summary_list[[paste0("Significant_probes_island_", gsub("[^A-Za-z0-9]", "_", island_name))]] <- island_counts[island_name]
    }
  }
}

# 转换为数据框并保存
summary_df <- data.frame(
  Metric = names(summary_list),
  Value = unlist(summary_list)
)

write.csv(summary_df, file.path(ANNOT_DIR, "annotation_summary.csv"), row.names = FALSE)
cat("   Annotation summary:", file.path(ANNOT_DIR, "annotation_summary.csv"), "\n")

# ============================================================================
# 7. 生成简单的注释可视化
# ============================================================================
cat("7. Generating annotation visualizations...\n")

if (require("ggplot2", quietly = TRUE)) {
  
  # 7.1 染色体分布条形图（如果有染色体信息）
  if ("chr" %in% colnames(annotated_results) && !all(is.na(annotated_results$chr))) {
    chr_data <- as.data.frame(table(annotated_results$chr))
    colnames(chr_data) <- c("Chromosome", "Count")
    
    # 按染色体顺序排序
    chr_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
    chr_data$Chromosome <- factor(chr_data$Chromosome, levels = chr_order)
    chr_data <- chr_data[order(chr_data$Chromosome), ]
    
    pdf(file.path(ANNOT_DIR, "chromosome_distribution.pdf"), width = 12, height = 6)
    p <- ggplot(chr_data, aes(x = Chromosome, y = Count)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      labs(title = "Distribution of DMPs by Chromosome",
           x = "Chromosome", y = "Number of DMPs") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
    dev.off()
    cat("      Saved: chromosome_distribution.pdf\n")
  }
  
  # 7.2 基因区域分布条形图
  if ("Gene_Region_Simple" %in% colnames(annotated_results) && 
      !all(is.na(annotated_results$Gene_Region_Simple))) {
    region_data <- as.data.frame(table(annotated_results$Gene_Region_Simple))
    colnames(region_data) <- c("Gene_Region", "Count")
    region_data <- region_data[region_data$Count > 0, ]
    
    if (nrow(region_data) > 0) {
      pdf(file.path(ANNOT_DIR, "gene_region_distribution.pdf"), width = 10, height = 6)
      p <- ggplot(region_data, aes(x = reorder(Gene_Region, -Count), y = Count)) +
        geom_bar(stat = "identity", fill = "darkgreen") +
        labs(title = "Distribution of DMPs by Gene Region",
             x = "Gene Region", y = "Number of DMPs") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
      dev.off()
      cat("      Saved: gene_region_distribution.pdf\n")
    }
  }
  
  # 7.3 CpG岛关系分布
  if ("Relation_to_UCSC_CpG_Island" %in% colnames(annotated_results) && 
      !all(is.na(annotated_results$Relation_to_UCSC_CpG_Island))) {
    island_data <- as.data.frame(table(annotated_results$Relation_to_UCSC_CpG_Island))
    colnames(island_data) <- c("CpG_Island_Relation", "Count")
    island_data <- island_data[island_data$Count > 0, ]
    
    if (nrow(island_data) > 0) {
      pdf(file.path(ANNOT_DIR, "cpg_island_relation_distribution.pdf"), width = 10, height = 6)
      p <- ggplot(island_data, aes(x = reorder(CpG_Island_Relation, -Count), y = Count)) +
        geom_bar(stat = "identity", fill = "purple") +
        labs(title = "Distribution of DMPs by CpG Island Relation",
             x = "Relation to CpG Island", y = "Number of DMPs") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
      dev.off()
      cat("      Saved: cpg_island_relation_distribution.pdf\n")
    }
  }
}

# ============================================================================
# 8. 完成
# ============================================================================
cat("\n")
cat("Annotation Summary:\n")
print(summary_df)
cat("\n")

cat("Annotation method used:", annotation_method, "\n")
cat("Output files saved to:", ANNOT_DIR, "\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n✓ Annotation step completed successfully!\n")

sink()

# 返回成功状态
quit(status = 0)