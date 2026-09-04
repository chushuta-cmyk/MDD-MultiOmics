# GSE53987 — 脑后转录组 (Microarray)

## 1. GEO 信息

| 项目 | 内容 |
|------|------|
| GEO Accession | GSE53987 |
| 标题 | Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia: a comparison of prefrontal cortex, striatum, and hippocampus among matched tetrads of controls with subjects diagnosed with schizophrenia, bipolar or major depressive disorder |
| 物种 | Homo sapiens |
| 平台 | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| 技术 | Microarray |
| 参考 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987 |
| 原始文献 | Lanz TA et al. Transl Psychiatry. 2019;9(1):151. doi:10.1038/s41398-019-0492-8 |

## 2. 样本设计

- **总样本数**：105（仅选取 MDD 与 Control 样本，原始研究为精神分裂症/双相/MDD/对照四联体设计）
- **脑区**：
  - 海马 (hippocampus)
  - 前额叶皮层 BA46 (Pre-frontal cortex / DLPFC)
  - 联合纹状体 (Associative striatum)
- **分组**：MDD vs Control（按脑区分别比较）

## 3. 在本项目中的角色

- **角色**：主要发现脑转录组数据集 (Primary discovery)
- **手稿章节**：§2.2 Brain transcriptomic analysis、§3.1 Brain transcriptomic alterations、Table 2、Figure 2
- **关键结果**：
  - 联合纹状体：303 DEGs（194 上调 / 109 下调）
  - BA46：611 DEGs（286 上调 / 325 下调）
  - 海马：215 DEGs（129 上调 / 86 下调）
  - 跨脑区 DEG 重叠极小，支持基因层面异质性结论

## 4. 分析流程

1. **数据加载**：从 series matrix（skip=66 行）读取表达矩阵，`ID_REF` 设为行名
2. **表型提取**：从 `!Sample_source_name_ch1` 行解析脑区与诊断标签，筛选 MDD/Control
3. **探针→基因映射**：`AnnotationDbi::mapIds` + `hgu133plus2.db`（PROBEID → SYMBOL，multiVals="first"）
4. **预处理过滤**：移除含 NA 的行；移除低表达探针（rowMeans > 中位数）
5. **差异表达**：按脑区分别运行 `limma`（lmFit → makeContrasts → contrasts.fit → eBayes），筛选 P.Value < 0.01
6. **结果整合**：按脑区合并 DEG 结果，输出 `GSE53987_MDD_vs_Control_BrainRegion_DEGs.csv`

## 5. 主要脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `GSE53987_analysis.R` | `complete_projects/GSE53987/scripts/` | 主分析脚本（加载→表型→limma DEG→探针映射→输出） |
| `brain_session_study.R` | `complete_projects/GSE53987/scripts/` | 辅助/探索性脚本 |

## 6. 输出文件

- `GSE53987_MDD_vs_Control_BrainRegion_DEGs.csv`：按脑区整合的 DEG 列表（GeneSymbol, logFC, P.Value, adj.P.Val, Region）
