# ============================================================================
# 1. 加载包
# ============================================================================
library(GEOquery)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(stringr)
library(DESeq2)
library(apeglm)
library(dplyr)

# ============================================================================
# 2. 读取并预处理 GSE80655 计数矩阵（Entrez ID -> Gene Symbol）
# ============================================================================
file_path <- "/Users/tan/Developer/projects/r_project/Study/Psychology_study/datasets/GSE80655_raw_counts_GRCh38.p13_NCBI.tsv"
counts_raw <- read.delim(file_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 转换 ID
symbols <- mapIds(org.Hs.eg.db, keys = rownames(counts_raw), column = "SYMBOL",
                  keytype = "ENTREZID", multiVals = "first")
valid <- !is.na(symbols)
counts_valid <- counts_raw[valid, ]
symbols_valid <- symbols[valid]

# 处理重复 Symbol（保留总计数最高的行）
gene_sums <- rowSums(counts_valid)
order_idx <- order(gene_sums, decreasing = TRUE)
symbols_ordered <- symbols_valid[order_idx]
counts_ordered <- counts_valid[order_idx, ]
unique_idx <- !duplicated(symbols_ordered)
final_counts <- counts_ordered[unique_idx, ]
rownames(final_counts) <- symbols_ordered[unique_idx]

cat("最终矩阵维度:", dim(final_counts), "\n")

# ============================================================================
# 3. 获取表型数据（通过 GEOquery）
# ============================================================================
gse <- getGEO("GSE80655", GSEMatrix = TRUE)[[1]]
pheno <- pData(gse)

# 提取所需列
pheno_sub <- pheno[, c("geo_accession", "brain region:ch1", "clinical diagnosis:ch1")]
colnames(pheno_sub) <- c("GSM_ID", "Brain_Region_Raw", "Diagnosis_Raw")

# 筛选 MDD 和 Control
keep <- pheno_sub$Diagnosis_Raw %in% c("Control", "Major Depression")
meta <- pheno_sub[keep, ]

# 清理名称
meta$Brain_Region <- make.names(gsub(" ", "_", meta$Brain_Region_Raw))
meta$Diagnosis <- make.names(gsub(" ", "_", meta$Diagnosis_Raw))
meta$Group <- factor(paste(meta$Brain_Region, meta$Diagnosis, sep = "_"))

# 按脑区分组
regions <- unique(meta$Brain_Region)
region_gsm <- lapply(regions, function(r) meta$GSM_ID[meta$Brain_Region == r])
names(region_gsm) <- regions

# ============================================================================
# 4. 将表达矩阵与元数据对齐（只保留匹配样本）
# ============================================================================
common_gsm <- intersect(meta$GSM_ID, colnames(final_counts))
expr <- final_counts[, common_gsm]
meta <- meta[match(colnames(expr), meta$GSM_ID), ]

# 转为整数（DESeq2 要求）
expr_int <- round(expr)
storage.mode(expr_int) <- "integer"

# ============================================================================
# 5. 对每个脑区运行 DESeq2（仅含 Diagnosis 协变量）
# ============================================================================
deg_list <- list()

for (r in regions) {
  cat("\n处理脑区:", r, "\n")
  gsm_r <- intersect(region_gsm[[r]], colnames(expr_int))
  if (length(gsm_r) < 2) { warning("样本不足，跳过"); next }
  
  counts_r <- expr_int[, gsm_r]
  meta_r <- meta[meta$GSM_ID %in% gsm_r, ]
  meta_r <- meta_r[match(colnames(counts_r), meta_r$GSM_ID), ]
  
  # 过滤零表达基因
  counts_r <- counts_r[rowSums(counts_r) > 0, ]
  
  dds <- DESeqDataSetFromMatrix(counts_r, colData = meta_r, design = ~ Diagnosis)
  dds$Diagnosis <- relevel(dds$Diagnosis, ref = "Control")
  dds <- DESeq(dds, quiet = TRUE)
  
  # 提取结果（使用 apeglm 收缩）
  res <- tryCatch(
    lfcShrink(dds, contrast = c("Diagnosis", "Major_Depression", "Control"), type = "apeglm"),
    error = function(e) results(dds, contrast = c("Diagnosis", "Major_Depression", "Control"))
  )
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("GeneSymbol") %>%
    filter(pvalue < 0.01) %>%
    rename(logFC = log2FoldChange, P.Value = pvalue, adj.P.Val = padj) %>%
    mutate(Region = r) %>%
    select(GeneSymbol, logFC, P.Value, adj.P.Val, Region)
  
  deg_list[[r]] <- res_df
  cat("  - DEG 数量 (P<0.01):", nrow(res_df), "\n")
}

# ============================================================================
# 6. 合并结果并保存
# ============================================================================
all_degs <- do.call(rbind, deg_list)
write.csv(all_degs, "GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv", row.names = FALSE)
cat("\n结果已保存。总 DEG 数:", nrow(all_degs), "\n")
