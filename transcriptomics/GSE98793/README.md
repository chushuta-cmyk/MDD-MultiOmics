# GSE98793 — 外周血转录组 & 免疫分析 (主)

## 1. GEO 信息

| 项目 | 内容 |
|------|------|
| GEO Accession | GSE98793 |
| 标题 | Gene expression analysis in whole blood samples obtained from donors diagnosed with Major Depressive Disorder compared with healthy controls |
| 物种 | Homo sapiens |
| 平台 | Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) |
| 技术 | Microarray |
| 参考 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98793 |
| 原始文献 | Leday GGR et al. Biol Psychiatry. 2018;83(1):70-80. doi:10.1016/j.biopsych.2017.01.021 |

## 2. 样本设计

- **总样本数**：192
- **MDD**：128
- **Control**：64
- **组织**：全血 (Whole blood)
- **探针数**：54,675
- **表达值**：已 log2 GC-RMA 标准化（无需重新标准化）
- **可用协变量**：Age、Sex、Batch、临床变量（含抗抑郁药暴露状态）

## 3. 在本项目中的角色

- **角色**：外周免疫分析主要数据集 (Peripheral immune profiling, primary)
- **手稿章节**：§2.5 Peripheral immune-related analysis、§3.5 Peripheral immune signatures、§3.6 Cross-tissue integration、Table 4、Figure 6、Supplementary Table S7
- **关键结果**（MCP-counter 反卷积，MDD vs Control）：
  - CD8+ T 细胞：Δ = −0.181，FDR = 0.0453（显著下调）
  - Cytotoxicity score：Δ = −0.167，FDR = 0.0453（显著下调）
  - 总 T 细胞：Δ = −0.107，FDR = 0.0478（显著下调）
  - 内皮细胞：Δ = +0.058，FDR = 0.0388（显著上调）
  - B 细胞：FDR = 0.165（无显著差异，FC ≈ 0.985，提示保留而非耗竭）
  - 免疫检查点：TIGIT 显著下调，PDCD1 不变（非典型耗竭表型）
  - AQP4–B 细胞：负相关（支持 BBB / 转运功能障碍假说）

## 4. 分析流程

1. **数据加载**：从 series matrix 读取（已 GC-RMA 标准化）
2. **预处理**：
   - 探针→基因映射（GPL570，`hgu133plus2.db`）
   - 重复探针保留最高平均表达
   - 移除无基因 Symbol 的探针
   - 移除低表达探针
3. **差异表达**：limma（MDD vs Control）
4. **免疫细胞反卷积**：
   - **MCP-counter**（主要方法，手稿报告）
   - **quanTIseq**（交叉验证）
   - 通过 `immunedeconv` R 包调用
5. **统计比较**：Wilcoxon rank-sum 检验，FDR（Benjamini-Hochberg）校正
6. **免疫检查点/耗竭标记分析**：TIGIT、PDCD1 等
7. **AQP4–B 细胞相关性分析**
8. **MOFA 导出**：标准化对象导出

## 5. 主要脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `GSE98793__GSE98793_Immune_Infiltration_Analysis.R` | `complete_projects/GSE98793/scripts/` | 免疫浸润分析主脚本（immunedeconv 多算法、统计比较、可视化） |
| `GSE98793__GSE98793_Complete_Analysis_Final.R` | `complete_projects/GSE98793/scripts/` | 完整分析流程（参考） |
| `GSE98793__GSE98793_Expression_Analysis_Complete.R` | `complete_projects/GSE98793/scripts/` | 表达预处理（参考） |
| `GSE98793__GSE98793_Analysis_Correct_Grouping.R` | `complete_projects/GSE98793/scripts/` | 分组校正（已合并） |
| `GSE98793_preprocessing_FIXED_v3.R` | `transcriptomics/GSE98793/scripts/` | 预处理修复版 v3（库路径、注释包、探针映射、select() 冲突修复） |
| `GSE98793_DATASET.md` | `complete_projects/GSE98793/` | 已有数据集说明文档 |

## 6. 输出文件

- `immune_results/Immune_Cell_Comparison.csv`：免疫细胞 MDD vs Control 统计比较
- `immune_results/Immune_Cell_Boxplot.pdf`：免疫细胞丰度箱线图
- `immune_results/Immune_Cell_Heatmap.pdf`：免疫细胞丰度热图
- `immune_results/Immune_Cell_FC_Plot.pdf`：log2(FC) 条形图
- `immune_results/Immune_Analysis_Summary.txt`：分析总结
- `futher_analysis/02_immune_validation/`：免疫验证结果（Table_S6_Immune_Infiltration.csv、exhaustion_marker_summary.csv、aqp4_bcell_correlation.csv、immune_key_findings.txt）
- `next_steps/MOFA/preprocessing/GSE98793/`：标准化对象（expression.rds, metadata.rds, feature_annotation.rds, validation_report.md）
