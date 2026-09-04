#!/usr/bin/env Rscript
# ============================================================================
# 05_enrichment_and_report.R
# GSE125105 富集分析与最终报告
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("05_enrichment_and_report.R: GSE125105 Enrichment Analysis and Final Report\n")
cat("============================================================================\n")
cat("\n")

# 设置绝对路径
RESULTS_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/results"
LOG_DIR <- "/data03/karama/projects/MDD_analysis/GSE125105/logs"

# 输出目录
ENRICH_DIR <- file.path(RESULTS_DIR, "05_enrichment")
REPORT_DIR <- file.path(RESULTS_DIR, "06_report")
dir.create(ENRICH_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(REPORT_DIR, showWarnings = FALSE, recursive = TRUE)

# 日志文件
log_file <- file.path(LOG_DIR, "05_enrichment_and_report.log")
sink(log_file, append = FALSE, split = TRUE)

cat("Start Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Enrichment Output Directory:", ENRICH_DIR, "\n")
cat("Report Output Directory:", REPORT_DIR, "\n")
cat("\n")

# ============================================================================
# 1. 加载必要包
# ============================================================================
cat("1. Loading required packages...\n")

# 基础包
required_packages <- c("data.table", "dplyr", "ggplot2", "stringr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 富集分析包
enrichment_available <- FALSE
if (!require("clusterProfiler", quietly = TRUE)) {
  cat("   clusterProfiler not available, attempting to install...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  tryCatch({
    BiocManager::install("clusterProfiler")
    library(clusterProfiler)
    enrichment_available <- TRUE
    cat("   ✓ clusterProfiler installed and loaded\n")
  }, error = function(e) {
    cat("   ⚠️  Failed to install clusterProfiler:", e$message, "\n")
    cat("   ⚠️  Enrichment analysis will be limited or skipped\n")
  })
} else {
  library(clusterProfiler)
  enrichment_available <- TRUE
  cat("   ✓ clusterProfiler loaded\n")
}

# 其他生物信息学包
if (enrichment_available) {
  if (!require("org.Hs.eg.db", quietly = TRUE)) {
    tryCatch({
      BiocManager::install("org.Hs.eg.db")
      library(org.Hs.eg.db)
      cat("   ✓ org.Hs.eg.db loaded\n")
    }, error = function(e) {
      cat("   ⚠️  Failed to load org.Hs.eg.db:", e$message, "\n")
    })
  } else {
    library(org.Hs.eg.db)
  }
  
  if (!require("DOSE", quietly = TRUE)) {
    tryCatch({
      BiocManager::install("DOSE")
      library(DOSE)
      cat("   ✓ DOSE loaded\n")
    }, error = function(e) {
      cat("   ⚠️  Failed to load DOSE:", e$message, "\n")
    })
  } else {
    library(DOSE)
  }
}

cat("✓ Packages loaded\n\n")

# ============================================================================
# 2. 加载数据
# ============================================================================
cat("2. Loading data from previous steps...\n")

# 2.1 加载注释的DMP结果
ANNOT_DIR <- file.path(RESULTS_DIR, "03_annotation")
annotated_file <- file.path(ANNOT_DIR, "GSE125105_DMP_annotated_all.csv")

if (file.exists(annotated_file)) {
  dmp_data <- fread(annotated_file)
  cat("   Annotated DMP results:", nrow(dmp_data), "probes\n")
} else {
  # 如果注释文件不存在，加载原始DMP结果
  DMP_DIR <- file.path(RESULTS_DIR, "02_dmp")
  dmp_file <- file.path(DMP_DIR, "GSE125105_DMP_all_results.csv")
  if (file.exists(dmp_file)) {
    dmp_data <- fread(dmp_file)
    cat("   DMP results (not annotated):", nrow(dmp_data), "probes\n")
  } else {
    stop("No DMP results found for enrichment analysis")
  }
}

# 2.2 加载显著性DMP
sig_dmp_file <- file.path(ANNOT_DIR, "GSE125105_DMP_annotated_significant.csv")
if (file.exists(sig_dmp_file)) {
  sig_dmp <- fread(sig_dmp_file)
} else {
  # 从完整结果中提取显著性DMP
  sig_dmp <- dmp_data[dmp_data$adj.P.Val < 0.05, ]
}

cat("   Significant DMPs (FDR < 0.05):", nrow(sig_dmp), "probes\n")
cat("   Hypermethylated DMPs (logFC > 0.1):", 
    sum(sig_dmp$adj.P.Val < 0.05 & sig_dmp$logFC > 0.1, na.rm = TRUE), "\n")
cat("   Hypomethylated DMPs (logFC < -0.1):", 
    sum(sig_dmp$adj.P.Val < 0.05 & sig_dmp$logFC < -0.1, na.rm = TRUE), "\n")

# 2.3 加载样本信息
QC_DIR <- file.path(RESULTS_DIR, "01_qc")
targets <- readRDS(file.path(QC_DIR, "targets.rds"))
cat("   Sample information:", nrow(targets), "samples\n")
cat("   Group distribution:\n")
print(table(targets$Group))
cat("\n")

# ============================================================================
# 3. 准备基因列表用于富集分析
# ============================================================================
cat("3. Preparing gene lists for enrichment analysis...\n")

# 检查是否有基因注释
if ("UCSC_RefGene_Name" %in% colnames(dmp_data)) {
  # 提取基因名称（可能多个基因用分号分隔）
  extract_genes <- function(gene_string) {
    if (is.na(gene_string) || gene_string == "") {
      return(character(0))
    }
    genes <- unlist(strsplit(as.character(gene_string), ";"))
    genes <- unique(trimws(genes))
    genes <- genes[genes != "" & !is.na(genes)]
    return(genes)
  }
  
  # 所有DMP的基因
  all_genes_list <- lapply(dmp_data$UCSC_RefGene_Name, extract_genes)
  all_genes <- unique(unlist(all_genes_list))
  cat("   Total genes in all DMPs:", length(all_genes), "\n")
  
  # 显著性DMP的基因
  sig_genes_list <- lapply(sig_dmp$UCSC_RefGene_Name, extract_genes)
  sig_genes <- unique(unlist(sig_genes_list))
  cat("   Genes in significant DMPs:", length(sig_genes), "\n")
  
  # 超甲基化DMP的基因
  hyper_dmp <- sig_dmp[sig_dmp$logFC > 0.1, ]
  hyper_genes_list <- lapply(hyper_dmp$UCSC_RefGene_Name, extract_genes)
  hyper_genes <- unique(unlist(hyper_genes_list))
  cat("   Genes in hypermethylated DMPs:", length(hyper_genes), "\n")
  
  # 低甲基化DMP的基因
  hypo_dmp <- sig_dmp[sig_dmp$logFC < -0.1, ]
  hypo_genes_list <- lapply(hypo_dmp$UCSC_RefGene_Name, extract_genes)
  hypo_genes <- unique(unlist(hypo_genes_list))
  cat("   Genes in hypomethylated DMPs:", length(hypo_genes), "\n")
  
} else {
  cat("   ⚠️  Gene annotation not available in DMP data\n")
  cat("   ⚠️  Cannot perform gene-based enrichment analysis\n")
  all_genes <- character(0)
  sig_genes <- character(0)
  hyper_genes <- character(0)
  hypo_genes <- character(0)
}

cat("\n")

# ============================================================================
# 4. GO富集分析
# ============================================================================
cat("4. Performing GO enrichment analysis...\n")

if (enrichment_available && length(sig_genes) > 0) {
  
  # 4.1 将基因符号转换为Entrez ID
  cat("   4.1 Converting gene symbols to Entrez IDs...\n")
  
  # 基因符号到Entrez ID的映射
  gene_symbols <- sig_genes
  entrez_map <- tryCatch({
    # 使用org.Hs.eg.db进行映射
    mapIds(org.Hs.eg.db, 
           keys = gene_symbols,
           column = "ENTREZID",
           keytype = "SYMBOL",
           multiVals = "first")
  }, error = function(e) {
    cat("   ⚠️  Gene symbol mapping failed:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(entrez_map)) {
    entrez_ids <- na.omit(entrez_map)
    cat("   Mapped", length(entrez_ids), "of", length(gene_symbols), "genes to Entrez IDs\n")
    
    # 4.2 GO生物过程（BP）富集
    cat("   4.2 GO Biological Process enrichment...\n")
    tryCatch({
      ego_bp <- enrichGO(gene = entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         readable = TRUE)
      
      if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
        cat("   GO BP: Found", nrow(ego_bp), "enriched terms\n")
        
        # 保存结果
        write.csv(as.data.frame(ego_bp), 
                  file.path(ENRICH_DIR, "GO_BP_enrichment.csv"),
                  row.names = FALSE)
        cat("      Saved: GO_BP_enrichment.csv\n")
        
        # 可视化：条形图
        pdf(file.path(ENRICH_DIR, "GO_BP_barplot.pdf"), width = 12, height = 8)
        if (nrow(ego_bp) > 20) {
          show_terms <- 20
        } else {
          show_terms <- nrow(ego_bp)
        }
        print(barplot(ego_bp, showCategory = show_terms, 
                      title = "GO Biological Process Enrichment"))
        dev.off()
        cat("      Saved: GO_BP_barplot.pdf\n")
        
        # 可视化：点图
        pdf(file.path(ENRICH_DIR, "GO_BP_dotplot.pdf"), width = 12, height = 8)
        print(dotplot(ego_bp, showCategory = show_terms,
                      title = "GO Biological Process Enrichment"))
        dev.off()
        cat("      Saved: GO_BP_dotplot.pdf\n")
      } else {
        cat("   GO BP: No significant enrichment found\n")
      }
    }, error = function(e) {
      cat("   ⚠️  GO BP enrichment failed:", e$message, "\n")
    })
    
    # 4.3 GO细胞组分（CC）富集
    cat("   4.3 GO Cellular Component enrichment...\n")
    tryCatch({
      ego_cc <- enrichGO(gene = entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         ont = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         readable = TRUE)
      
      if (!is.null(ego_cc) && nrow(ego_cc) > 0) {
        cat("   GO CC: Found", nrow(ego_cc), "enriched terms\n")
        
        write.csv(as.data.frame(ego_cc), 
                  file.path(ENRICH_DIR, "GO_CC_enrichment.csv"),
                  row.names = FALSE)
        cat("      Saved: GO_CC_enrichment.csv\n")
      } else {
        cat("   GO CC: No significant enrichment found\n")
      }
    }, error = function(e) {
      cat("   ⚠️  GO CC enrichment failed:", e$message, "\n")
    })
    
    # 4.4 GO分子功能（MF）富集
    cat("   4.4 GO Molecular Function enrichment...\n")
    tryCatch({
      ego_mf <- enrichGO(gene = entrez_ids,
                         OrgDb = org.Hs.eg.db,
                         ont = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1,
                         readable = TRUE)
      
      if (!is.null(ego_mf) && nrow(ego_mf) > 0) {
        cat("   GO MF: Found", nrow(ego_mf), "enriched terms\n")
        
        write.csv(as.data.frame(ego_mf), 
                  file.path(ENRICH_DIR, "GO_MF_enrichment.csv"),
                  row.names = FALSE)
        cat("      Saved: GO_MF_enrichment.csv\n")
      } else {
        cat("   GO MF: No significant enrichment found\n")
      }
    }, error = function(e) {
      cat("   ⚠️  GO MF enrichment failed:", e$message, "\n")
    })
    
  } else {
    cat("   ⚠️  Not enough genes mapped for enrichment analysis\n")
  }
  
} else {
  cat("   ⚠️  Enrichment analysis not available or no genes to analyze\n")
}

cat("\n")

# ============================================================================
# 5. KEGG通路富集分析
# ============================================================================
cat("5. Performing KEGG pathway enrichment analysis...\n")

if (enrichment_available && length(sig_genes) > 0 && exists("entrez_ids") && length(entrez_ids) > 0) {
  tryCatch({
    kk <- enrichKEGG(gene = entrez_ids,
                     organism = 'hsa',
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1)
    
    if (!is.null(kk) && nrow(kk) > 0) {
      cat("   KEGG: Found", nrow(kk), "enriched pathways\n")
      
      # 保存结果
      write.csv(as.data.frame(kk), 
                file.path(ENRICH_DIR, "KEGG_enrichment.csv"),
                row.names = FALSE)
      cat("      Saved: KEGG_enrichment.csv\n")
      
      # 可视化
      pdf(file.path(ENRICH_DIR, "KEGG_barplot.pdf"), width = 12, height = 8)
      if (nrow(kk) > 20) {
        show_terms <- 20
      } else {
        show_terms <- nrow(kk)
      }
      print(barplot(kk, showCategory = show_terms, 
                    title = "KEGG Pathway Enrichment"))
      dev.off()
      cat("      Saved: KEGG_barplot.pdf\n")
      
      pdf(file.path(ENRICH_DIR, "KEGG_dotplot.pdf"), width = 12, height = 8)
      print(dotplot(kk, showCategory = show_terms,
                    title = "KEGG Pathway Enrichment"))
      dev.off()
      cat("      Saved: KEGG_dotplot.pdf\n")
      
    } else {
      cat("   KEGG: No significant enrichment found\n")
    }
  }, error = function(e) {
    cat("   ⚠️  KEGG enrichment failed:", e$message, "\n")
  })
} else {
  cat("   ⚠️  KEGG enrichment not available\n")
}

cat("\n")

# ============================================================================
# 6. 疾病本体（DO）富集分析
# ============================================================================
cat("6. Performing Disease Ontology (DO) enrichment analysis...\n")

if (enrichment_available && length(sig_genes) > 0 && exists("entrez_ids") && length(entrez_ids) > 0) {
  if (require("DOSE", quietly = TRUE)) {
    tryCatch({
      do <- enrichDO(gene = entrez_ids,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
      
      if (!is.null(do) && nrow(do) > 0) {
        cat("   DO: Found", nrow(do), "enriched disease terms\n")
        
        write.csv(as.data.frame(do), 
                  file.path(ENRICH_DIR, "DO_enrichment.csv"),
                  row.names = FALSE)
        cat("      Saved: DO_enrichment.csv\n")
        
        # 检查是否与抑郁症相关
        depression_terms <- do[grep("depress|mood|bipolar|mental", 
                                    do$Description, ignore.case = TRUE), ]
        if (nrow(depression_terms) > 0) {
          cat("   Found", nrow(depression_terms), "depression-related terms:\n")
          print(depression_terms$Description)
          write.csv(depression_terms, 
                    file.path(ENRICH_DIR, "depression_related_terms.csv"),
                    row.names = FALSE)
          cat("      Saved: depression_related_terms.csv\n")
        }
        
      } else {
        cat("   DO: No significant enrichment found\n")
      }
    }, error = function(e) {
      cat("   ⚠️  DO enrichment failed:", e$message, "\n")
    })
  } else {
    cat("   ⚠️  DOSE package not available for DO enrichment\n")
  }
} else {
  cat("   ⚠️  DO enrichment not available\n")
}

cat("\n")

# ============================================================================
# 7. 基因组区域富集分析（基于CpG岛关系）
# ============================================================================
cat("7. Performing genomic region enrichment analysis...\n")

if ("Relation_to_UCSC_CpG_Island" %in% colnames(dmp_data)) {
  # 计算所有DMP中CpG岛关系的分布
  all_island_dist <- table(dmp_data$Relation_to_UCSC_CpG_Island, useNA = "ifany")
  all_island_dist <- all_island_dist / sum(all_island_dist) * 100
  
  # 计算显著性DMP中CpG岛关系的分布
  sig_island_dist <- table(sig_dmp$Relation_to_UCSC_CpG_Island, useNA = "ifany")
  if (sum(sig_island_dist) > 0) {
    sig_island_dist <- sig_island_dist / sum(sig_island_dist) * 100
  }
  
  # 创建比较数据框
  island_comparison <- data.frame(
    Relation = names(all_island_dist),
    All_DMPs = as.numeric(all_island_dist),
    Significant_DMPs = as.numeric(sig_island_dist[match(names(all_island_dist), 
                                                        names(sig_island_dist))])
  )
  island_comparison$Significant_DMPs[is.na(island_comparison$Significant_DMPs)] <- 0
  
  write.csv(island_comparison, 
            file.path(ENRICH_DIR, "CpG_island_enrichment.csv"),
            row.names = FALSE)
  cat("      Saved: CpG_island_enrichment.csv\n")
  
  # 可视化
  island_long <- reshape2::melt(island_comparison, id.vars = "Relation")
  colnames(island_long) <- c("Relation", "Group", "Percentage")
  
  pdf(file.path(ENRICH_DIR, "CpG_island_enrichment_plot.pdf"), width = 12, height = 6)
  p <- ggplot(island_long, aes(x = Relation, y = Percentage, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "CpG Island Relation Enrichment",
         subtitle = "Comparison between all DMPs and significant DMPs",
         x = "Relation to CpG Island",
         y = "Percentage (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  dev.off()
  cat("      Saved: CpG_island_enrichment_plot.pdf\n")
  
  cat("   CpG island relation analysis completed\n")
} else {
  cat("   ⚠️  CpG island relation information not available\n")
}

cat("\n")

# ============================================================================
# 8. 生成最终报告
# ============================================================================
cat("8. Generating final report...\n")

# 8.1 收集所有步骤的统计信息
cat("   8.1 Collecting statistics from all steps...\n")

# 从QC步骤
qc_summary_file <- file.path(QC_DIR, "qc_summary.csv")
if (file.exists(qc_summary_file)) {
  qc_stats <- fread(qc_summary_file)
} else {
  qc_stats <- data.frame(Metric = c("QC summary not available"), Value = c("NA"))
}

# 从DMP步骤
dmp_stats_file <- file.path(RESULTS_DIR, "02_dmp", "dmp_statistics_summary.csv")
if (file.exists(dmp_stats_file)) {
  dmp_stats <- fread(dmp_stats_file)
} else {
  dmp_stats <- data.frame(Metric = c("DMP statistics not available"), Value = c("NA"))
}

# 从注释步骤
annot_summary_file <- file.path(ANNOT_DIR, "annotation_summary.csv")
if (file.exists(annot_summary_file)) {
  annot_stats <- fread(annot_summary_file)
} else {
  annot_stats <- data.frame(Metric = c("Annotation summary not available"), Value = c("NA"))
}

# 8.2 创建综合报告
report_text <- paste0(
  "================================================================================\n",
  "                     GSE125105 METHYLATION ANALYSIS REPORT\n",
  "================================================================================\n",
  "\n",
  "Analysis completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "Project directory: /data03/karama/projects/MDD_analysis/GSE125105\n",
  "\n",
  "----------------------------------------------------------------\n",
  "1. DATA OVERVIEW\n",
  "----------------------------------------------------------------\n",
  "\n"
)

# 添加QC统计
report_text <- paste0(report_text, "Quality Control Statistics:\n")
for (i in 1:nrow(qc_stats)) {
  report_text <- paste0(report_text, "  ", qc_stats$Metric[i], ": ", qc_stats$Value[i], "\n")
}

# 添加样本信息
report_text <- paste0(report_text, "\nSample Information:\n")
group_table <- table(targets$Group)
for (group_name in names(group_table)) {
  report_text <- paste0(report_text, "  ", group_name, ": ", group_table[group_name], " samples\n")
}

# 添加DMP统计
report_text <- paste0(report_text, "\n----------------------------------------------------------------\n",
                      "2. DIFFERENTIAL METHYLATION RESULTS\n",
                      "----------------------------------------------------------------\n\n")

for (i in 1:min(10, nrow(dmp_stats))) {
  report_text <- paste0(report_text, "  ", dmp_stats$Metric[i], ": ", dmp_stats$Value[i], "\n")
}

# 添加富集分析摘要
report_text <- paste0(report_text, "\n----------------------------------------------------------------\n",
                      "3. ENRICHMENT ANALYSIS SUMMARY\n",
                      "----------------------------------------------------------------\n\n")

# 检查富集结果文件
enrichment_files <- list.files(ENRICH_DIR, pattern = "_enrichment\\.csv$", full.names = TRUE)
if (length(enrichment_files) > 0) {
  report_text <- paste0(report_text, "Enrichment analyses performed:\n")
  for (file in enrichment_files) {
    analysis_name <- gsub("_enrichment\\.csv$", "", basename(file))
    analysis_name <- gsub("_", " ", analysis_name)
    enrich_data <- fread(file)
    report_text <- paste0(report_text, "  ", analysis_name, ": ", nrow(enrich_data), " significant terms\n")
  }
} else {
  report_text <- paste0(report_text, "No significant enrichment results found or enrichment not performed.\n")
}

# 添加输出文件列表
report_text <- paste0(report_text, "\n----------------------------------------------------------------\n",
                      "4. OUTPUT FILES\n",
                      "----------------------------------------------------------------\n\n")

report_text <- paste0(report_text, "Main output directories:\n")
results_dirs <- list.dirs(RESULTS_DIR, full.names = FALSE, recursive = FALSE)
results_dirs <- results_dirs[results_dirs != ""]
for (dir_name in sort(results_dirs)) {
  dir_path <- file.path(RESULTS_DIR, dir_name)
  file_count <- length(list.files(dir_path))
  report_text <- paste0(report_text, "  ", dir_name, "/: ", file_count, " files\n")
}

# 添加关键文件
report_text <- paste0(report_text, "\nKey result files:\n")
key_files <- c(
  "01_qc/qc_summary.csv",
  "02_dmp/GSE125105_DMP_significant.csv",
  "03_annotation/GSE125105_DMP_annotated_significant.csv",
  "04_visualization/volcano_plot_enhanced.pdf",
  "05_enrichment/GO_BP_enrichment.csv"
)

for (file_path in key_files) {
  full_path <- file.path(RESULTS_DIR, file_path)
  if (file.exists(full_path)) {
    report_text <- paste0(report_text, "  ✓ ", file_path, "\n")
  } else {
    report_text <- paste0(report_text, "  ✗ ", file_path, " (not found)\n")
  }
}

# 添加建议
report_text <- paste0(report_text, "\n----------------------------------------------------------------\n",
                      "5. NEXT STEPS & RECOMMENDATIONS\n",
                      "----------------------------------------------------------------\n\n",
                      "1. Review significant DMPs in '03_annotation/GSE125105_DMP_annotated_significant.csv'\n",
                      "2. Examine enrichment results for biological insights\n",
                      "3. Validate top hits in independent datasets if available\n",
                      "4. Consider functional validation of key genes/pathways\n",
                      "5. Integrate with other omics data for systems-level understanding\n",
                      "\n",
                      "================================================================================\n",
                      "                            END OF REPORT\n",
                      "================================================================================\n"
)

# 8.3 保存报告
report_file <- file.path(REPORT_DIR, "run_summary.md")
writeLines(report_text, report_file)
cat("   Final report saved:", report_file, "\n")

# 8.4 创建HTML版本（如果rmarkdown可用）
if (require("rmarkdown", quietly = TRUE)) {
  html_report <- gsub("\\.md$", ".html", report_file)
  tryCatch({
    rmarkdown::render(report_file, output_file = html_report, quiet = TRUE)
    cat("   HTML report saved:", html_report, "\n")
  }, error = function(e) {
    cat("   ⚠️  Could not generate HTML report:", e$message, "\n")
  })
}

# 8.5 保存会话信息
session_info <- capture.output(sessionInfo())
writeLines(session_info, file.path(REPORT_DIR, "final_session_info.txt"))
cat("   Session info saved: final_session_info.txt\n")

# ============================================================================
# 9. 完成
# ============================================================================
cat("\n")
cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
cat("\n")
cat("Summary:\n")
cat("- Total probes analyzed:", nrow(dmp_data), "\n")
cat("- Significant DMPs (FDR < 0.05):", nrow(sig_dmp), "\n")
cat("- Enrichment analyses performed:", length(enrichment_files), "\n")
cat("- Output directories created:", length(results_dirs), "\n")
cat("\n")
cat("Key output locations:\n")
cat("  QC results:", file.path(RESULTS_DIR, "01_qc"), "\n")
cat("  DMP results:", file.path(RESULTS_DIR, "02_dmp"), "\n")
cat("  Annotation:", file.path(RESULTS_DIR, "03_annotation"), "\n")
cat("  Visualizations:", file.path(RESULTS_DIR, "04_visualization"), "\n")
cat("  Enrichment:", file.path(RESULTS_DIR, "05_enrichment"), "\n")
cat("  Final report:", file.path(RESULTS_DIR, "06_report"), "\n")
cat("\n")
cat("End Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("\n✓✓✓ GSE125105 analysis pipeline completed successfully! ✓✓✓\n")

sink()

# 返回成功状态
quit(status = 0)