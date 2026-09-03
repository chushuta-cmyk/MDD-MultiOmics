# GSE54563 — 脑转录组验证集 (Microarray)

## 1. GEO 信息

| 项目 | 内容 |
|------|------|
| GEO Accession | GSE54563 |
| 标题 | Gene expression in anterior cingulate cortex in major depressive disorder |
| 物种 | Homo sapiens |
| 平台 | Illumina HumanHT-12 V3.0 expression beadchip (GPL6947) |
| 技术 | Microarray |
| 参考 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54563 |
| 原始文献 | Chang LC et al. PLoS One. 2014;9(3):e90980. doi:10.1371/journal.pone.0090980 |

## 2. 样本设计

- **总样本数**：50
- **MDD**：25
- **Control**：25
- **设计**：MDD-Control 匹配对
- **脑区**：前扣带回皮层 (Anterior Cingulate Cortex, ACC)
- **探针数**：48,803

## 3. 在本项目中的角色

- **角色**：独立验证数据集 (Independent validation)
- **手稿章节**：§2.2（验证）、§3.1（基因层面验证）、Table 3、Figure 2C、Supplementary Table S3
- **关键结果**：
  - p < 0.05：757 个探针水平信号
  - p < 0.01：92 个
  - p < 0.001：4 个
  - 平均 logFC：0.000；中位 logFC：−0.0177
  - 平均 Cohen's d：0.181；中位 Cohen's d：0.149
  - 与主要发现数据集（GSE53987, GSE80655）在基因层面重叠有限
  - 功能富集突出髓鞘化 (myelination)、少突胶质细胞功能、神经递质信号通路

## 4. 分析流程

1. **数据加载**：从 series matrix 读取表达矩阵与表型
2. **差异表达**：limma（MDD vs Control），计算 logFC、t 统计量、Cohen's d
3. **基因选择策略**：
   - 方案 1：仅 P 值（P < 0.05 / 0.01 / 0.001）
   - 方案 2：logFC + P 值双指标（|logFC| > 0.5 & P < 0.05）
   - 方案 3：宽松标准（|logFC| > 0.3 & P < 0.05）
   - 方案 4：仅 logFC（|logFC| > 1.0，因数据 logFC 较大）
   - **最终采用**：|logFC| > 1.0（因平均 logFC 大，adj.P 接近 1，多重检验校正过于严格）
4. **基因注释**：`illuminaHumanv3.db`（PROBEID → SYMBOL）
5. **GO 富集**：`clusterProfiler::enrichGO`（BP，pvalueCutoff=0.05，BH 校正）
6. **血液 vs 脑对比**：与 GSE98793 血液数据进行统计对比（样本量、平台、显著基因数、效应量）
7. **MOFA 导出**：独立脚本导出标准化对象

## 5. 主要脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `GSE54563__GSE54563_FinalAnalysis_PvalueFC_Integration.R` | `complete_projects/GSE54563/scripts/` | 最终分析脚本（P+logFC 双指标选择、GO 富集、血液-脑对比、综合报告） |
| `GSE54563__GSE54563_BugFix_Complete_Analysis.R` | `complete_projects/GSE54563/scripts/` | Bug 修复版完整分析 |
| `GSE54563__GSE54563_Correct_Phenotyping_Analysis.R` | `complete_projects/GSE54563/scripts/` | 表型校正版分析 |
| `GSE54563__Visualization_Code.R` | `complete_projects/GSE54563/scripts/` | 可视化代码 |
| `GSE54563_standalone_MOFA_v2.R` | `transcriptomics/GSE54563/scripts/` | 独立 MOFA 输入导出（不依赖先前 RStudio 对象） |

## 6. 输出文件

- `GSE54563_Significant_DEGs_P005_FC10.csv`：显著 DEG 列表（|logFC|>1.0）
- `GSE54563_GO_Enrichment.csv`：GO 富集结果
- `Blood_vs_Brain_Comparison.csv`：GSE98793（血液）vs GSE54563（脑）统计对比
- `GSE54563_Analysis_Status.csv`：分析状态汇总
- `Final_Statistics_Comparison.csv`：最终统计对比
- `next_steps/MOFA/preprocessing/GSE54563/`：标准化对象（expression.rds, metadata.rds, feature_annotation.rds, validation_report.md）
