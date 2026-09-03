# GSE88890 — 脑 DNA 甲基化 (探索性)

## 1. GEO 信息

| 项目 | 内容 |
|------|------|
| GEO Accession | GSE88890 |
| 标题 | Methylomic profiling of cortex samples from completed suicide cases implicates a role for PSORS1C3 in major depression and suicide |
| 物种 | Homo sapiens |
| 平台 | Illumina HumanMethylation450 BeadChip |
| 技术 | DNA 甲基化微阵列 |
| 参考 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88890 |
| 原始文献 | Murphy TM et al. Transl Psychiatry. 2017;7(1):e989. doi:10.1038/tp.2016.249 |

## 2. 样本设计

- **总样本数**：75
- **MDD（自杀）**：37
- **Control**：38
- **脑区**：
  - BA11（眶额叶皮层）：40 样本
  - BA25（膝下前扣带回）：35 样本
- **探针数**：Illumina 450K（约 485,000 探针）

## 3. 在本项目中的角色

- **角色**：探索性脑甲基化 (Exploratory brain methylation)
- **手稿章节**：§2.6（探索性分析，因信号弱未纳入主分析，仅作补充解释）
- **关键结果**：
  - 显著 DMPs（FDR < 0.05）：**14 个**
  - 超甲基化：10 个（71.4%）
  - 低甲基化：4 个（28.6%）
  - **信号弱**：DMR 分析、通路富集未完成/跳过
  - 因显著 DMP 极少，脑甲基化结果未纳入主分析，仅作补充

## 4. 分析流程

1. **数据加载与 QC**：CHAMP 流程（`GSE88890_With_GEO2R_Groups_CHAMP.*.R` 多个版本）
2. **DMP 分析**：limma，FDR < 0.05
3. **注释**：探针→基因映射
4. **可视化**：QC 密度图、火山图、Top30 DMP 热图
5. **MOFA 导出**：基因水平 beta 矩阵导出（`GSE88890_standalone_MOFA_v2.R`）

## 5. 主要脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `GSE88890__GSE88890_*.R`（13+ 个版本） | `complete_projects/GSE88890/scripts/` | CHAMP 完整分析（含 Auto_Grouping、Complete_Analysis、Fixed_ChAMP、Improved、Offline、Ultimate_Fixed、With_Actual_groupid、With_GEO2R_Groups、With_SeriesMatrix_Metadata 等版本） |
| `GSE88890_With_GEO2R_Groups_CHAMP.fixed.R` | `GSE88890_full/` | 修复版 CHAMP 分析 |
| `GSE88890_With_GEO2R_Groups_CHAMP.fullrerun.R` | `GSE88890_full/` | 完整重跑版 |
| `GSE88890_With_GEO2R_Groups_CHAMP.server.R` | `GSE88890_full/` | 服务器版 |
| `GSE88890_standalone_MOFA_v2.R` | `transcriptomics/GSE88890/scripts/` | 独立 MOFA 输入导出（基因水平 beta 矩阵） |

## 6. 输出文件

- `GSE88890_full/FINAL_REPORT.txt`：最终分析报告
- `GSE88890_full/results/`：DMP 结果、注释、火山图、QC 图
- `GSE88890_full/data/`：数据文件
- `GSE88890_full/logs/`：日志
- `next_steps/MOFA/preprocessing/GSE88890/`：标准化对象（expression.rds, metadata.rds, feature_annotation.rds, validation_report.md, dataset_summary.json）
