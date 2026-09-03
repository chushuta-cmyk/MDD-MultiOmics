# GSE39653 — 外周血转录组 & 免疫分析 (次/验证)

## 1. GEO 信息

| 项目 | 内容 |
|------|------|
| GEO Accession | GSE39653 |
| 标题 | Inflammation and neurological disease-related genes are differentially expressed in depressed patients with mood disorders and correlate with morphometric and functional imaging abnormalities |
| 物种 | Homo sapiens |
| 平台 | Illumina HumanHT-12 V4.0 expression beadchip (GPL10558) |
| 技术 | Microarray |
| 参考 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39653 |
| 原始文献 | Savitz J et al. Brain Behav Immun. 2013;31:161-171. doi:10.1016/j.bbi.2012.10.007 |

## 2. 样本设计

- **组织**：PBMC / 全血
- **原始分组**：MDD、Control、BD（Bipolar Disorder，双相情感障碍）
- **本项目使用**：仅选取 MDD 与 Control（排除 BD）
- **样本 ID 格式**：title 列为 "MDD-XXX"、"HC-XXX"、"BD2-XXX" 格式，需从前缀提取诊断
- **注释文件**：需手动准备 `GPL10558_annotation.txt`（含 ProbeID 与 SYMBOL 列）

## 3. 在本项目中的角色

- **角色**：外周免疫分析次要/验证数据集 (Peripheral immune profiling, secondary/validation)
- **手稿章节**：§2.5（方法中列出）；**§3.5 结果中目前仅报告 GSE98793，GSE39653 待显式补充**（Tier 1 修订第 7 条）
- **关键结果**：
  - B 细胞验证：与脑组织（GSE80655，B cell FC = 0.759, p = 0.019）进行方向一致性评估
  - MCP-counter 反卷积：所有免疫细胞类型 MDD vs Control 比较（含 FDR 校正）
  - 作为 GSE98793 外周免疫发现的独立验证

## 4. 分析流程

1. **数据加载**：从 series matrix 读取（手动指定行范围：表型 32–62 行，表达矩阵 64–43180 行）
2. **表型解析**：从 `title` 列前缀提取诊断（MDD / HC / BD），筛选 MDD vs Control（排除 BD）
3. **表达矩阵清洗**：过滤坏行（列数不匹配的行），`ID_REF` 设为行名
4. **探针→基因映射**：
   - 使用 GPL10558 注释文件（非 Bioconductor 包，因 `illuminaHumanv4.db` 依赖问题）
   - 自动检测 ProbeID 列（ID / ID_REF / probe_id 等）与 SYMBOL 列
   - 重复基因取平均
5. **免疫细胞反卷积**：`immunedeconv::deconvolute(method = "mcp_counter")`
6. **B 细胞比较**：
   - 提取 B 细胞评分（自动检测行名 "B cell" / "B_cell" 等）
   - MDD vs Control Wilcoxon 检验
   - 计算 FC、均值、中位数、SD
7. **与脑组织结果对比**：与 GSE80655 脑 B 细胞结果（FC=0.759, p=0.019）进行方向/幅度/显著性一致性评估
8. **所有细胞类型比较**：遍历所有免疫细胞类型，Wilcoxon 检验 + BH FDR 校正

## 5. 主要脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `GSE39653__GSE39653_Priority1_Final.R` | `complete_projects/GSE39653/scripts/` | 主分析脚本 v2.0（B 细胞验证完整流程，无 illuminaHumanv4.db 依赖，使用 GPL 注释文件） |
| `GSE39653__GSE39653_draft.R` | `complete_projects/GSE39653/scripts/` | 草稿版 |
| `Future_Work_Complete_Pipeline.R` | `complete_projects/GSE39653/scripts/` | 未来工作完整管线 |

## 6. 输出文件

- `GSE39653_Blood_Immune_Scores.csv`：所有免疫细胞类型评分
- `Blood_BCell_Validation_Result.csv`：B 细胞 MDD vs Control 比较结果（含与脑结果对比）
- `GSE39653_Blood_AllCellTypes_Comparison.csv`：所有细胞类型统计比较（含 FDR）
- `GSE39653_Blood_BCells_Boxplot.pdf`：B 细胞箱线图
