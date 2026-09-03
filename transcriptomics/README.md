# Transcriptomics & Methylation Datasets

本目录汇总 MDD 多组学项目中使用的全部 GEO 数据集，按 GSE accession 分子目录组织。每个子目录包含该数据集的预处理脚本、MOFA 输入导出脚本，以及一份 README 说明其在项目中的角色、样本设计、分析流程和关键结果。

## 数据集总览

| GSE | 数据类型 | 平台 | 组织/脑区 | 样本数 (MDD/Ctrl) | 项目角色 | 手稿章节 |
|-----|---------|------|----------|-------------------|---------|---------|
| [GSE53987](./GSE53987/) | 脑转录组 (microarray) | Affymetrix HG-U133 Plus 2.0 (GPL570) | 海马、BA46、联合纹状体 | 105 | 主要发现 | §2.2, §3.1, Table 2 |
| [GSE80655](./GSE80655/) | 脑转录组 (RNA-seq) | NCBI GRCh38.p13 raw counts | AnCg、DLPFC、nAcc | 281 | 主要发现 | §2.2, §3.1, Table 2 |
| [GSE54563](./GSE54563/) | 脑转录组 (microarray) | Illumina HumanHT-12 V3 (GPL6947) | 前扣带回皮层 (ACC) | 50 (25/25) | 独立验证 | §2.2, §3.1, Table 3 |
| [GSE98793](./GSE98793/) | 外周血转录组 | Affymetrix HG-U133 Plus 2.0 (GPL570) | 全血 | 192 (128/64) | 外周免疫分析 (主) | §2.5, §3.5, Table 4 |
| [GSE39653](./GSE39653/) | 外周血转录组 | Illumina HumanHT-12 V4 (GPL10558) | PBMC/全血 | MDD+Ctrl+BD (BD排除) | 外周免疫分析 (次/验证) | §2.5 (结果待补) |
| [GSE125105](./GSE125105/) | DNA 甲基化 | Illumina HumanMethylation450K | 全血 | — | 甲基化主分析 | §2.6, §3.7 |
| [GSE88890](./GSE88890/) | DNA 甲基化 | Illumina HumanMethylation450K | 脑 (BA11, BA25) | 75 (37/38) | 探索性脑甲基化 | §2.6 (探索性) |
| [GSE113725](./GSE113725/) | DNA 甲基化 | — | 全血 | — | 评估后排除 (统计功效不足) | §2.6 (排除) |

## 通用分析流程

### 脑转录组 (GSE53987, GSE80655, GSE54563)
1. 数据加载（series matrix 或 raw counts）
2. 探针/基因 ID 映射 → Gene Symbol（去重，保留最高表达）
3. 按脑区分割样本
4. 差异表达分析（limma 用于 microarray，DESeq2 用于 RNA-seq）
5. GO / Reactome 富集（clusterProfiler, ReactomePA）
6. WGCNA 共表达网络分析
7. 跨脑区收敛基因识别

### 外周血转录组 & 免疫 (GSE98793, GSE39653)
1. 数据加载与预处理
2. 免疫细胞反卷积（MCP-counter, quanTIseq via immunedeconv）
3. MDD vs Control 统计比较（Wilcoxon, FDR 校正）
4. 免疫检查点/耗竭标记分析
5. AQP4–B 细胞相关性分析

### DNA 甲基化 (GSE125105, GSE88890)
1. 数据加载与 QC（CHAMP / minfi 流程）
2. DMP 鉴定（limma，adjusted p < 0.05）
3. 探针→基因注释映射
4. 与转录组收敛基因 / WGCNA hub 基因的重叠分析

## MOFA 整合

所有数据集均通过 `*_standalone_MOFA_v*.R` 或 `*_to_MOFA_*.R` 脚本导出为标准化对象（expression.rds / metadata.rds / feature_annotation.rds），存放于 `next_steps/MOFA/inputs/`，供 MOFA+ 多组学因子分析使用。
