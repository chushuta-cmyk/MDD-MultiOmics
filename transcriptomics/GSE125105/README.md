# GSE125105 — 外周血 DNA 甲基化 (主分析)

## 1. GEO 信息

| 项目 | 内容 |
|------|------|
| GEO Accession | GSE125105 |
| 标题 | DeepWAS: multivariate genotype-phenotype associations by directly integrating regulatory information using deep learning（关联 MDD 血液甲基化） |
| 物种 | Homo sapiens |
| 平台 | Illumina HumanMethylation450 BeadChip |
| 技术 | DNA 甲基化微阵列 |
| 参考 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125105 |
| 原始文献 | Arloth J et al. PLoS Comput Biol. 2020;16(2):e1007616. doi:10.1371/journal.pcbi.1007616 |

## 2. 样本设计

- **组织**：全血 (Whole blood)
- **分组**：MDD vs Control
- **甲基化相关基因数**：4,488（DMP 映射到基因后）

## 3. 在本项目中的角色

- **角色**：DNA 甲基化主分析 (Primary methylation analysis)
- **手稿章节**：§2.6 Peripheral blood DNA methylation analysis、§3.7 Epigenetic regulation preferentially targets network hub genes、Figure 7、Supplementary Table S8–S9
- **关键结果**：
  - 甲基化相关基因：4,488 个
  - 与跨脑区收敛 DEG（≥3 脑区，n=10）重叠：**仅 1 个基因（COX19）**
  - 与 WGCNA hub 基因（n=616）重叠：**125 个基因**
  - 重叠 hub 基因包括：SYNE1、NTRK2、STAT1、ITGB1、PICALM、PRKCB、DNAH10、RAB2A、FLOT1、PRKAR2A、AGBL2、CCDC33 等
  - **核心结论**：表观遗传调控优先靶向网络枢纽基因，而非持续差异表达的基因

## 4. 分析流程

1. **数据加载与 QC**：
   - `01_load_and_qc.R` / `01_load_and_qc_original.R`：加载甲基化矩阵，质量控制
   - 内存优化脚本（大矩阵处理）
2. **DMP 分析**：
   - `02_dmp_analysis.R` / `02_dmp_analysis_simple.R`：limma 差异甲基化位点分析
   - `GSE125105_DMP_Limma_Only_fixed.R`：仅 limma 的 DMP 分析（修复版）
   - 阈值：adjusted p-value < 0.05
3. **注释**：
   - `03_annotation.R`：DMP 探针→基因注释映射
4. **可视化**：
   - `04_visualization.R`：火山图、热图等
5. **富集与报告**：
   - `05_enrichment_and_report.R`：功能富集与综合报告
6. **跨组学重叠分析**：
   - `futher_analysis/03_methylation_overlap/03_methylation_overlap.R`：甲基化基因与收敛 DEG / WGCNA hub 基因的重叠分析

## 5. 主要脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `01_load_and_qc.R` / `01_load_and_qc_original.R` | `complete_projects/GSE125105/pipeline/` | 数据加载与 QC |
| `02_dmp_analysis.R` / `02_dmp_analysis_simple.R` | `complete_projects/GSE125105/pipeline/` | DMP 分析 |
| `GSE125105_DMP_Limma_Only_fixed.R` | `complete_projects/GSE125105/pipeline/` | 仅 limma DMP（修复版） |
| `03_annotation.R` | `complete_projects/GSE125105/pipeline/` | 探针→基因注释 |
| `04_visualization.R` | `complete_projects/GSE125105/pipeline/` | 可视化 |
| `05_enrichment_and_report.R` | `complete_projects/GSE125105/pipeline/` | 富集与报告 |
| `GSE125105__GSE125105_CHAMP_*.R`（多个） | `complete_projects/GSE125105/scripts/` | CHAMP 流程脚本 |
| `03_methylation_overlap.R` | `futher_analysis/03_methylation_overlap/` | 跨组学重叠分析 |

## 6. 输出文件

- `futher_analysis/03_methylation_overlap/methylation_gene_list.txt`：甲基化相关基因列表（4,488）
- `futher_analysis/03_methylation_overlap/methylation_overlap_convergent.csv`：甲基化×收敛 DEG 重叠（COX19）
- `futher_analysis/03_methylation_overlap/methylation_overlap_hub.csv`：甲基化×WGCNA hub 重叠（125 基因）
- `futher_analysis/03_methylation_overlap/methylation_overlap_summary.txt`：重叠分析总结
- `next_steps/MOFA/inputs/GSE88890_methylation_gene_level.csv/rds`：MOFA 甲基化输入（注：文件名标为 GSE88890，实际需确认）
