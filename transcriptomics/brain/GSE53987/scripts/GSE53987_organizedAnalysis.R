# ============================================================================
# 1. 加载必要的包
# ============================================================================
library(GEOquery)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(hgu133plus2.db)
library(stringr)  

# ============================================================================
# 2. 读取原始表达矩阵（跳过注释行，保留ID_REF作为行名）
# ============================================================================
file_path <- "/Users/tan/Developer/projects/r_project/Study/Psychology_study/datasets/GSE53987_series_matrix.txt"

raw_data <- read.delim(
  file = file_path,
  skip = 66,
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  check.names = FALSE,  # 保持GSM ID原始格式
  as.is = TRUE          # 防止ID_REF转换为因子
)

# 将ID_REF列设为行名，并移除该列
if ("ID_REF" %in% colnames(raw_data)) {
  rownames(raw_data) <- raw_data[, "ID_REF"]
  expr_matrix <- as.matrix(raw_data[, -which(colnames(raw_data) == "ID_REF")])
  storage.mode(expr_matrix) <- "numeric"
} else {
  stop("未找到 'ID_REF' 列，请检查 skip 行数。")
}

cat("表达矩阵维度:", dim(expr_matrix), "\n")

# ============================================================================
# 3. 从文件注释行提取样本分组信息（!Sample_source_name_ch1）
# ============================================================================
geo_lines <- readLines(file_path)
phenotype_line <- grep("^!Sample_source_name_ch1", geo_lines, value = TRUE)
if (length(phenotype_line) == 0) {
  stop("未找到以 '!Sample_source_name_ch1' 开头的行。")
}

# 解析标签：去掉首列标签和双引号
labels_split <- strsplit(phenotype_line, "\t")[[1]]
labels_cleaned <- gsub("\"", "", labels_split[-1])

# ============================================================================
# 4. 筛选 MDD 和 Control 样本
# ============================================================================
target_pattern <- "major depressive disorder|control"
selected_idx <- grepl(target_pattern, labels_cleaned)
selected_labels <- labels_cleaned[selected_idx]

# 将标签拆分为脑区和疾病状态
split_mat <- stringr::str_split_fixed(selected_labels, ", ", n = 2)
brain_area <- split_mat[, 1]
disease <- split_mat[, 2]

# 对于简单的 "hippocampus" 标签，手动补全
is_simple <- (disease == "")
brain_area[is_simple] <- "hippocampus"
disease[is_simple] <- split_mat[is_simple, 1]

# 创建最终分组因子 (脑区_疾病)
final_group <- factor(paste(brain_area, disease, sep = "_"))
levels(final_group) <- gsub(" ", "_", levels(final_group))
levels(final_group) <- gsub("\\(|\\)", "", levels(final_group))

cat("分组样本计数:\n")
print(table(final_group))

# ============================================================================
# 5. 构建样本元数据（GSM ID、脑区、疾病）
# ============================================================================
all_gsm <- colnames(expr_matrix)
selected_gsm <- all_gsm[selected_idx]

# 规范化脑区名称（与后续 regions 列表保持一致）
brain_area_clean <- gsub(" ", "_", brain_area)
brain_area_clean <- gsub("\\(|\\)", "", brain_area_clean)
brain_area_clean <- gsub("_BA46", "_cortex_BA46", brain_area_clean)
brain_area_clean[grepl("Pre-frontal_cortex", brain_area_clean)] <- "Pre-frontal_cortex_BA46"
brain_area_clean[grepl("Associative_striatum", brain_area_clean)] <- "Associative_striatum"
brain_area_clean[grepl("hippocampus", brain_area_clean)] <- "hippocampus"

sample_metadata <- data.frame(
  GSM_ID = selected_gsm,
  Group = final_group,
  Brain_Area = brain_area_clean,
  Disease = disease,
  stringsAsFactors = FALSE
)

# 按脑区拆分 GSM ID 列表
regions <- c("hippocampus", "Pre-frontal_cortex_BA46", "Associative_striatum")
regions_gsm_list <- list()
for (r in regions) {
  regions_gsm_list[[r]] <- sample_metadata$GSM_ID[sample_metadata$Brain_Area == r]
  cat(r, "样本数:", length(regions_gsm_list[[r]]), "\n")
}

# ============================================================================
# 6. 表达矩阵预处理（去除NA、低表达探针）
# ============================================================================
expr_matrix <- expr_matrix[rowSums(is.na(expr_matrix)) == 0, ]
low_cutoff <- median(expr_matrix)
expr_matrix <- expr_matrix[rowMeans(expr_matrix) > low_cutoff, ]

# 子集化到选中的样本（MDD和Control）
expr_sub <- expr_matrix[, selected_idx]

# ============================================================================
# 7. 按脑区进行差异表达分析（limma）
# ============================================================================
deg_list <- list()

for (region in regions) {
  cat("\n处理脑区:", region, "\n")
  
  # 获取当前脑区的样本
  gsm_ids <- regions_gsm_list[[region]]
  expr_region <- expr_sub[, colnames(expr_sub) %in% gsm_ids, drop = FALSE]
  meta_region <- sample_metadata[sample_metadata$GSM_ID %in% gsm_ids, ]
  
  # 确保列顺序与元数据一致
  expr_region <- expr_region[, meta_region$GSM_ID, drop = FALSE]
  
  # 构造分组因子（使用make.names处理组名，避免特殊字符）
  group <- factor(meta_region$Group, levels = unique(meta_region$Group))
  group_valid <- factor(group, labels = make.names(levels(group)))
  
  # 设计矩阵
  design <- model.matrix(~0 + group_valid)
  colnames(design) <- gsub("group_valid", "", colnames(design))
  
  # 拟合线性模型
  fit <- lmFit(expr_region, design)
  
  # 构建对比：MDD vs Control
  mdd_name <- make.names(paste0(region, "_major_depressive_disorder"))
  ctrl_name <- make.names(paste0(region, "_control"))
  
  if (!all(c(mdd_name, ctrl_name) %in% colnames(design))) {
    warning("在", region, "中未找到MDD或Control组，跳过。")
    next
  }
  
  contrast_mat <- makeContrasts(contrasts = paste(mdd_name, ctrl_name, sep = " - "),
                                levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast_mat))
  
  # 提取结果（P < 0.01，未调整）
  tt <- topTable(fit2, coef = 1, number = Inf, adjust.method = "fdr")
  tt_sig <- tt[tt$P.Value < 0.01, ]
  cat("  - 显著探针数 (P<0.01):", nrow(tt_sig), "\n")
  
  if (nrow(tt_sig) == 0) {
    deg_list[[region]] <- data.frame()
    next
  }
  
  # 探针ID转基因符号（hgu133plus2.db）
  probe_ids <- rownames(tt_sig)
  symbols <- suppressWarnings(
    mapIds(hgu133plus2.db, keys = probe_ids, column = "SYMBOL",
           keytype = "PROBEID", multiVals = "first")
  )
  
  # 若映射失败较多，尝试转为大写再映射
  if (sum(!is.na(symbols)) < 10) {
    symbols_upper <- suppressWarnings(
      mapIds(hgu133plus2.db, keys = toupper(probe_ids), column = "SYMBOL",
             keytype = "PROBEID", multiVals = "first")
    )
    if (sum(!is.na(symbols_upper)) > sum(!is.na(symbols))) {
      symbols <- symbols_upper
    }
  }
  
  # 整理结果
  res <- tt_sig[, c("logFC", "P.Value", "adj.P.Val")]
  res$ProbeID <- probe_ids
  res$GeneSymbol <- symbols
  res$Region <- region
  res <- res[!is.na(res$GeneSymbol), ]   # 移除无基因符号的探针
  
  deg_list[[region]] <- res
  cat("  - 最终有基因符号的DEG数:", nrow(res), "\n")
}

# ============================================================================
# 8. 合并所有脑区结果并保存
# ============================================================================
all_degs <- do.call(rbind, deg_list)
cat("\n总DEG数（全部脑区）:", nrow(all_degs), "\n")

output_file <- "GSE53987_MDD_vs_Control_BrainRegion_DEGs.csv"
write.csv(all_degs, file = output_file, row.names = FALSE)
cat("结果已保存至:", output_file, "\n")
