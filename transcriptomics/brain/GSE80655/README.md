# GSE80655 — 脑转录组 (RNA-seq)

## 1. GEO 信息

| 项目 | 内容 |
|------|------|
| GEO Accession | GSE80655 |
| 标题 | Post-mortem molecular profiling of three psychiatric disorders |
| 物种 | Homo sapiens |
| 平台 | RNA-seq (raw counts, NCBI GRCh38.p13) |
| 技术 | RNA-seq |
| 参考 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655 |
| 原始文献 | Ramaker RC et al. Genome Med. 2017;9:72. doi:10.1186/s13073-017-0458-5 |

## 2. 样本设计

- **总样本数**：281（MOFA 预处理验证：39,376 基因 × 281 样本）
- **脑区**：
  - 前扣带回皮层 (AnCg / anterior cingulate cortex)
  - 背外侧前额叶皮层 (DLPFC)
  - 伏隔核 (nAcc / nucleus accumbens)
- **分组**：MDD (Major Depression) vs Control

## 3. 在本项目中的角色

- **角色**：主要发现脑转录组数据集 (Primary discovery)
- **手稿章节**：§2.2 Brain transcriptomic analysis、§3.1 Brain transcriptomic alterations、Table 2、Figure 2
- **关键结果**：
  - AnCg：621 DEGs（293 上调 / 328 下调）
  - DLPFC：293 DEGs（101 上调 / 192 下调）
  - nAcc：372 DEGs（217 上调 / 155 下调）
  - WGCNA：DLPFC 多个显著模块（Modules 18, 25, 43, 44, 45），hub 基因包括 ATM、CLPTM1、ACER3、FUS、MAP1LC3A

## 4. 分析流程

1. **数据加载**：读取 raw counts TSV（`GSE80655_raw_counts_GRCh38.p13_NCBI.tsv`），GeneID 设为行名
2. **ID 映射**：`AnnotationDbi::mapIds` + `org.Hs.eg.db`（ENTREZID → SYMBOL），移除未映射基因
3. **基因去重**：重复 Symbol 保留总计数最高的探针/行
4. **表型获取**：`GEOquery::getGEO` 获取 pData，筛选 `clinical diagnosis:ch1` 为 Control / Major Depression
5. **样本匹配**：表达矩阵列名与 GSM ID 对齐，按脑区分割
6. **差异表达**：按脑区分别运行 `DESeq2`（DESeqDataSetFromMatrix → DESeq → lfcShrink with apeglm），筛选 pvalue < 0.01
7. **结果整合**：输出 `GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv`
8. **MOFA 导出**：VST 标准化后导出为 `GSE80655_VST_normalized.rds`

## 5. 主要脚本

| 脚本 | 位置 | 用途 |
|------|------|------|
| `GSE80655__GSE80655.R` | `complete_projects/GSE80655/scripts/` | 主分析脚本（raw counts→DESeq2 DEG→输出） |
| `01_loader.R` / `02_preprocess.R` / `03_normalization.R` / `04_quality_control.R` | `complete_projects/GSE80655/` | 模块化预处理脚本 |
| `GSE80655_to_MOFA_noGEO.R` | `transcriptomics/GSE80655/scripts/` | 独立 MOFA 输入导出（VST 标准化，不依赖 GEOquery） |

## 6. 输出文件

- `GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv`：按脑区整合的 DEG 列表
- `GSE80655_VST_normalized.rds`：VST 标准化表达矩阵（MOFA 输入）
- `next_steps/MOFA/preprocessing/GSE80655/`：标准化对象（expression.rds, metadata.rds, feature_annotation.rds）
