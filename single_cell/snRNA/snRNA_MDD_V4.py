import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scrublet as scr
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.io as sio
from harmonypy import run_harmony
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pertpy as pt

# ==============================================================================
# SECTION 1: Path Configurations & Environment Setup
# ==============================================================================
data_dir = "/root/bioinformatics/GSE213982/data/"
plot_dir = os.path.join(data_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)

mtx_path = os.path.join(data_dir, "GSE213982_combined_counts_matrix.mtx.gz")
genes_path = os.path.join(data_dir, "GSE213982_combined_counts_matrix_genes_rows.csv.gz")
cells_path = os.path.join(data_dir, "GSE213982_combined_counts_matrix_cells_columns.csv.gz")

# 配置 Scanpy 绘图参数
sc.settings.set_figure_params(dpi=150, facecolor='white', color_map='viridis')

# ==============================================================================
# SECTION 2: Sparse Data Matrix Loading & Object Initialization
# ==============================================================================
print(">>> Loading sparse matrix...")
X = sio.mmread(mtx_path).T.tocsr() 

genes = pd.read_csv(genes_path, index_col=0) 
cells = pd.read_csv(cells_path, index_col=0)

print(">>> Constructing AnnData...")
adata = ad.AnnData(X=X, obs=cells, var=genes)
adata.var_names_make_unique()

# ==============================================================================
# SECTION 3: Metadata Structuring, Quote Stripping & Pattern Inference (REVISED)
# ==============================================================================
print(">>> Reconstructing metadata and stripping quotes...")

# 1. Filter out the header leftover row '"x"' or 'x'
valid_cells = ~adata.obs_names.isin(['"x"', 'x'])
adata = adata[valid_cells].copy()

# 2. Completely strip double quotes from barcodes (e.g., '"F1.AAAC...' -> 'F1.AAAC...')
adata.obs_names = adata.obs_names.map(lambda bar: str(bar).replace('"', ''))

# --- NEW: Fix the pandas ambiguity error ---
# If a column named 'x' exists, drop it
if 'x' in adata.obs.columns:
    adata.obs = adata.obs.drop(columns=['x'])

# Rename the index to something explicit so Seaborn doesn't get confused
adata.obs.index.name = 'cell_barcode'
# -------------------------------------------

# 3. Tokenize the cleaned barcode to extract sample_id
# Structure: sample_id.barcode.cell_type_major
adata.obs['sample_id'] = adata.obs_names.str.split('.').str[0]
adata.obs['orig_ident'] = adata.obs['sample_id']

# 4. Infer sex from the donor prefix (F -> Female, M -> Male)
adata.obs['sex'] = adata.obs['sample_id'].map(lambda s_id: 'Female' if str(s_id).startswith('F') else 'Male')

# 5. Map diagnosis based on GSE213982 actual MDD sample list
mdd_samples = [
    'F3', 'F5', 'F7', 'F8', 'F12', 'F14', 'F17', 'F25', 'F33', 'F35', 
    'M3', 'M6', 'M11', 'M15', 'M17', 'M20', 'M23', 'M31', 'M34'
] 

adata.obs['diagnosis'] = adata.obs['sample_id'].map(lambda s_id: 'MDD' if s_id in mdd_samples else 'Control')
adata.obs['diagnosis'] = adata.obs['diagnosis'].astype('category')

print("=== Metadata Reconstruction Verification ===")
print(adata.obs[['sample_id', 'sex', 'diagnosis']].head())
print(f"\nCell count per Sex group:\n{adata.obs['sex'].value_counts()}")
print(f"Cell count per Diagnosis group:\n{adata.obs['diagnosis'].value_counts()}\n")

# ==============================================================================
# SECTION 4: QC Processing & Graphical Metrics Export
# ==============================================================================
print(">>> Quality Control Processing...")
adata.var['mt'] = adata.var_names.str.contains('^mt-', case=False, regex=True)
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Graphics Output 1: Pre-filtering Violin Diagnostic Plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "01_qc_metrics_pre_filter.png"), bbox_inches='tight')
plt.close()

print(f"   Before QC Filter Count: {adata.n_obs} cells")
adata = adata[
    (adata.obs['n_genes_by_counts'] > 500) &
    (adata.obs['n_genes_by_counts'] < 8000) &
    (adata.obs['total_counts'] > 1000) &
    (adata.obs['pct_counts_mt'] < 10), :
].copy()
print(f"   After QC Filter Count: {adata.n_obs} cells")

# ==============================================================================
# SECTION 5: Mask-Safe Doublet Detection (Per Sample)
# ==============================================================================
print(">>> Scrublet Doublet Detection...")
adata.obs['doublet_score'] = 0.0
adata.obs['is_doublet'] = False

for sample in adata.obs['sample_id'].unique():
    sample_mask = adata.obs['sample_id'] == sample
    sample_adata = adata[sample_mask]
    
    scrub = scr.Scrublet(sample_adata.X)
    d_scores, d_preds = scrub.scrub_doublets()
    
    # 使用显式内部掩码安全定位分配，防止行序列错位
    adata.obs.loc[sample_mask, 'doublet_score'] = d_scores
    adata.obs.loc[sample_mask, 'is_doublet'] = d_preds

# Graphics Output 2: Doublet Score Group Distribution across Batches
plt.figure(figsize=(12, 4))
sns.violinplot(data=adata.obs, x='sample_id', y='doublet_score', hue='is_doublet', split=True)
plt.xticks(rotation=90)
plt.title("Doublet Score Stratification per Sample Batch")
plt.savefig(os.path.join(plot_dir, "02_scrublet_distribution.png"), bbox_inches='tight')
plt.close()

adata = adata[~adata.obs['is_doublet']].copy()
print(f"   After Doublet Removal Count: {adata.n_obs} cells")

# ==============================================================================
# CRITICAL FIX: Store raw counts BEFORE any normalization
# ==============================================================================
print("\n>>> DIAGNOSTIC: Checking raw count matrix before normalization...")
print(f"   X matrix type: {type(adata.X)}")
print(f"   X matrix shape: {adata.X.shape}")
print(f"   X matrix dtype: {adata.X.dtype}")

# Check a few values
if hasattr(adata.X, 'toarray'):
    sample_vals = adata.X[:5, :5].toarray()
else:
    sample_vals = adata.X[:5, :5]
print(f"   Sample values from X (first 5x5):\n{sample_vals}")
print(f"   Min value in X: {adata.X.min()}")
print(f"   Max value in X: {adata.X.max()}")
print(f"   Mean value in X: {adata.X.mean()}")

# Store raw counts - this MUST happen before normalization
adata.layers['counts'] = adata.X.copy()

print(f"\n   Raw counts stored in layers['counts']")
print(f"   layers['counts'] type: {type(adata.layers['counts'])}")

# ==============================================================================
# SECTION 6: Log Transformation, Scaleless PCA, & Direct Harmony Core Fix
# ==============================================================================
print("\n>>> Normalization and Dimensional Transformation...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

print(f"\n>>> DIAGNOSTIC: After normalization...")
print(f"   X matrix is now normalized/log-transformed")
if hasattr(adata.X, 'toarray'):
    sample_vals_norm = adata.X[:5, :5].toarray()
else:
    sample_vals_norm = adata.X[:5, :5]
print(f"   Sample normalized values (first 5x5):\n{sample_vals_norm}")

# Verify raw counts are still preserved
if hasattr(adata.layers['counts'], 'toarray'):
    sample_vals_raw = adata.layers['counts'][:5, :5].toarray()
else:
    sample_vals_raw = adata.layers['counts'][:5, :5]
print(f"   Sample RAW count values from layers['counts'] (first 5x5):\n{sample_vals_raw}")

sc.pp.scale(adata, max_value=10, zero_center=False)
sc.tl.pca(adata, svd_solver='arpack')

# Graphics Output 3: Pre-harmony Structural Batch Variance PCA Map
sc.pl.pca(adata, color='sample_id', show=False)
plt.savefig(os.path.join(plot_dir, "03_pca_pre_harmony.png"), bbox_inches='tight')
plt.close()

print(">>> Executing Explicit Harmony Batch Integration...")
pca_matrix = adata.obsm['X_pca']
ho = run_harmony(pca_matrix, adata.obs, ['sample_id'])

# -------------------------------------------------------------------------
# DYNAMIC SHAPE FIX: VERSION-PROOF MATRIX ALIGNMENT
# -------------------------------------------------------------------------
harmony_embedding = ho.Z_corr
if harmony_embedding.shape[0] != adata.n_obs:
    print("   [Auto-Fix] Transposing Harmony output to match Cell dimensions...")
    harmony_embedding = harmony_embedding.T
else:
    print("   [Auto-Fix] Harmony output already aligns with Cell dimensions. Skipping transpose.")

adata.obsm['X_pca_harmony'] = harmony_embedding
print(f"   Harmony Alignment Completed. Dimension Shape Matrix: {adata.obsm['X_pca_harmony'].shape}")
# -------------------------------------------------------------------------

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40, use_rep='X_pca_harmony')
sc.tl.umap(adata)

# Graphics Output 4: Post-Integration Batch Mixing Verification Plots
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
sc.pl.umap(adata, color='sample_id', ax=axes[0], show=False, title="UMAP: Structured Sample Batches")
sc.pl.umap(adata, color='sex', ax=axes[1], show=False, title="UMAP: Stratified Sex Variables")
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "04_umap_post_harmony_alignment.png"), bbox_inches='tight')
plt.close()

# ==============================================================================
# SECTION 7 & 8: Graph Resolution Partitioning & Multi-Metadata Mapping (FIXED)
# ==============================================================================
print(">>> Topographical Clustering Analysis...")
sc.tl.leiden(adata, resolution=0.8)

cluster_to_celltype = {
    "0": "ExN", "1": "ExN", "2": "Oligo", "3": "ExN", "4": "InN",
    "5": "Astro", "6": "OPC", "7": "Micro", "8": "InN_SST", 
    "9": "InN_PVALB", "10": "Endo", "11": "InN_VIP", "12": "ExN", 
    "13": "Oligo", "14": "Pericyte"
}
adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_celltype).astype('category')

# Graphics Output 5: Cell Type Topological Atlas Maps
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sc.pl.umap(adata, color='cell_type', legend_loc='on data', ax=axes[0], show=False, title="dlPFC Cell Typing Annotations")
sc.pl.umap(adata, color='diagnosis', ax=axes[1], show=False, title="Stratified Phenotype (MDD vs Control)")
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "05_umap_annotated_topology.png"), bbox_inches='tight')
plt.close()

adata.write(os.path.join(data_dir, "combined_annotated.h5ad"))

# ==============================================================================
# SECTION 9: Aggregated Cellular Matrix Processing & Linear Mixed Modeling (PyDESeq2)
# ==============================================================================
print("\n>>> Running Pseudobulk Processing Sequence (Pure Pandas/NumPy Engine)...")
de_results = []

# 获取所有唯一的细胞类型
cell_types = adata.obs['cell_type'].dropna().unique()

for ct in cell_types:
    print(f"\n" + "="*60)
    print(f"--- Processing Cell Type: {ct} ---")
    print("="*60)
    
    # 提取当前细胞类型的单细胞子集
    ct_adata = adata[adata.obs['cell_type'] == ct]
    
    # -------------------------------------------------------------------------
    # 检查点 1: 单细胞子集规模验证
    # -------------------------------------------------------------------------
    print(f"  [CHECK 1] Single-Cell Subset Info:")
    print(f"    - Cell Count (n_obs): {ct_adata.n_obs}")
    print(f"    - Unique Samples (n_samples): {ct_adata.obs['sample_id'].nunique()}")
    
    if ct_adata.n_obs < 50 or ct_adata.obs['sample_id'].nunique() < 6:
        print(f"    [SKIP] Insufficient data footprint profiles for {ct}. Skipping...")
        continue
    
    # -------------------------------------------------------------------------
    # 检查点 2: 原始计数矩阵提取与类型检查
    # -------------------------------------------------------------------------
    print("  [CHECK 2] Extracting Raw Integer Counts from layers['counts']...")
    raw_counts = ct_adata.layers['counts']
    
    # 如果是稀疏矩阵，转换为稠密数组
    if hasattr(raw_counts, 'toarray'):
        raw_counts = raw_counts.toarray()
        
    print(f"    - Matrix Type: {type(raw_counts)}")
    print(f"    - Matrix Shape: {raw_counts.shape}")
    print(f"    - Value Range: Min={raw_counts.min()}, Max={raw_counts.max()}, Mean={raw_counts.mean():.4f}")
    
    # 潜在风险校验：确保数据是未做过 Log 转换的原始整数 counts
    if raw_counts.max() < 50 and raw_counts.mean() < 1:
        print("    [WARNING] Max value is suspiciously low. Ensure this is RAW count, not normalized data!")

    # -------------------------------------------------------------------------
    # 检查点 3: 伪块 (Pseudobulk) 样本合并矩阵检查
    # -------------------------------------------------------------------------
    print("  [CHECK 3] Executing Native Pandas Pseudobulk Aggregation (Summing rows by sample_id)...")
    temp_df = pd.DataFrame(
        raw_counts, 
        index=ct_adata.obs['sample_id'], 
        columns=ct_adata.var_names
    )
    
    # 按 sample_id 进行分组求和
    pseudobulk_counts = temp_df.groupby(level=0).sum()
    # 强制转换为整型，并将空值填 0
    pseudobulk_counts = pseudobulk_counts.fillna(0).astype(int)
    
    print(f"    - Aggregated Matrix Shape (Samples x Genes): {pseudobulk_counts.shape}")
    print("    - First 3 rows and 5 columns of Aggregated Matrix:")
    print(pseudobulk_counts.iloc[:3, :5])
    
    # -------------------------------------------------------------------------
    # 检查点 4: 低表达/无变异基因过滤检查
    # -------------------------------------------------------------------------
    print("  [CHECK 4] Filtering low-expression genes (Total counts across all samples >= 10)...")
    gene_sums = pseudobulk_counts.sum(axis=0)
    filtered_counts = pseudobulk_counts.loc[:, gene_sums >= 10]
    
    print(f"    - Before Filtering Gene Count: {pseudobulk_counts.shape[1]}")
    print(f"    - After Filtering Gene Count: {filtered_counts.shape[1]}")
    print(f"    - Filtered Matrix Sparsity: {(filtered_counts == 0).sum().sum() / filtered_counts.size:.2%}")
    
    if filtered_counts.shape[1] == 0:
        print(f"    [SKIP] No genes passed filtering threshold for {ct}. Skipping...")
        continue

    # -------------------------------------------------------------------------
    # 检查点 5: 元数据对齐与临床实验设计验证
    # -------------------------------------------------------------------------
    print("  [CHECK 5] Aligning Clinical Metadata with Aggregated Matrix...")
    # 提取唯一的样本临床信息
    sample_info = ct_adata.obs[['sample_id', 'diagnosis', 'sex']].drop_duplicates().set_index('sample_id')
    # 严格按照 counts 矩阵的样本顺序进行切片对齐
    aligned_meta = sample_info.loc[filtered_counts.index]
    
    print("    - Aligned Design Matrix Metadata Table:")
    print(aligned_meta)
    
    # 校验临床对比组是否存在
    diagnosis_counts = aligned_meta['diagnosis'].value_counts()
    print(f"    - Group Distribution: {dict(diagnosis_counts)}")
    
    if len(diagnosis_counts) < 2:
        print("    [CRITICAL SKIP] Missing one of the comparison groups (Control or MDD) in this subset. Skipping...")
        continue

    # -------------------------------------------------------------------------
    # 执行 PyDESeq2 统计模型建模（完全替换后的现代化稳健代码）
    # -------------------------------------------------------------------------
    print("  [MODEL] Initializing PyDESeq2 DeseqDataSet...")
    try:
        # 1. 使用现代化设计公式，彻底消除 DeprecationWarning
        dds = DeseqDataSet(
            counts=filtered_counts,
            metadata=aligned_meta,
            design="~diagnosis"
        )
        
        print("  [MODEL] Fitting DESeq2 Model Parameters (Size factors, Dispersions, LFCs)...")
        dds.deseq2()
        
        print("  [MODEL] Extracting Statistical Testing Metrics via DeseqStats...")
        stat_res = DeseqStats(dds, contrast=["diagnosis", "MDD", "Control"])
        
        # 2. 核心修复：必须显式调用 summary() 激活内部计算，否则没有 results_df 属性
        stat_res.summary()
        
        # 3. 此时读取结果矩阵，绝对不会再触发 AttributeError
        res_df = stat_res.results_df.copy()
        
        # 4. 提取、保存与流转数据
        res_df['cell_type'] = ct
        res_df['gene'] = res_df.index
        de_results.append(res_df)
        
        # 打印单群成功状态
        sig_counts = (res_df['padj'] < 0.05).sum()
        print(f"    [SUCCESS] DE pipeline finished for {ct}. Significantly altered genes (padj < 0.05): {sig_counts}")
            
    except Exception as e:
        # 捕获异常，确保某一细胞类型收敛失败时（如 Astro 报拟合不收敛），自动跳过并处理下一个，不让整个脚本退出
        print(f"    [CRITICAL ERROR] DE Pipeline Failure on {ct}: {e}")

# ==============================================================================
# SECTION 9 后续：结果合并与可视化
# ==============================================================================
if de_results:
    final_de = pd.concat(de_results, ignore_index=True) # 改变 1: 强制重置索引，防止基因名索引冲突导致合并错位
    final_de.to_csv(os.path.join(data_dir, "pseudobulk_de_results.csv"), index=False)
    print(f"\n>>> Global Pseudobulk DE Analysis complete. Total results saved: {final_de.shape[0]} rows.")
else:
    final_de = pd.DataFrame() # 保证后面画图引用时变量存在
    print("\n>>> WARNING: No DE results were generated across any cell type.")

# ==============================================================================
# SECTION 10: Cellular Abundance Analysis
# ==============================================================================
print("\n>>> Cellular Abundance Analysis...")
props_df = pd.crosstab(adata.obs['sample_id'], adata.obs['cell_type'], normalize='index')

sample_diagnosis = adata.obs[['sample_id', 'diagnosis']].drop_duplicates().set_index('sample_id')
props_df = props_df.reset_index()
props_df['diagnosis'] = props_df['sample_id'].map(sample_diagnosis['diagnosis'])

props_melted = props_df.melt(
    id_vars=['sample_id', 'diagnosis'], 
    var_name='cell_type', 
    value_name='proportion'
)

plt.figure(figsize=(14, 6))
sns.boxplot(data=props_melted, x='cell_type', y='proportion', hue='diagnosis')
plt.xticks(rotation=45, ha='right')
plt.title("Cell Type Proportions: MDD vs Control")
plt.tight_layout()
plt.savefig(os.path.join(plot_dir, "10_celltype_proportions.png"), bbox_inches='tight')
plt.close()

# ==============================================================================
# SECTION 11: Main Composite Panel Production Block
# ==============================================================================
print(">>> Building Multi-Panel Consolidated Figure File...")
fig = plt.figure(figsize=(24, 10)) # 优化：略微调大画布物理长宽，给右侧堆叠图的图例留出空间
gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.1, 1.2]) # 优化：平衡子图比例

# Ax1: UMAP
ax1 = fig.add_subplot(gs[0])
sc.pl.umap(adata, color='cell_type', legend_loc='on data', title="Annotated Single-Nuclei Atlas", ax=ax1, show=False)

# Ax2: DE 点图
ax2 = fig.add_subplot(gs[1])
if de_results and not final_de.empty:
    # 改变 2: 显式过滤非 NaN 的 padj 值，避免由于统计未收敛产生的空值导致 groupby 报错
    valid_de = final_de[final_de['padj'].notna()].copy()
    sig_de = valid_de[valid_de['padj'] < 0.05].copy()
    
    if not sig_de.empty:
        # 改变 3: 现代化 Pandas 语法适配。旧版 groupby().apply() 会引入多重索引，导致后续画图找不到列名
        top_de = sig_de.sort_values('padj').groupby('cell_type').head(3).reset_index(drop=True)
        top_de['-log10(padj)'] = -np.log10(top_de['padj'])
        
        sns.scatterplot(data=top_de, x='cell_type', y='gene', size='-log10(padj)', hue='log2FoldChange', palette='coolwarm', ax=ax2)
        ax2.set_title("Top Transcriptomic Alterations (MDD vs Ctrl)")
        ax2.tick_params(axis='x', rotation=45)
        # 将点图的图例移到子图下方，防止挤压右侧
        ax2.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)
    else:
        ax2.text(0.5, 0.5, 'No Statistically Significant DE Results\n(padj < 0.05)', ha='center', va='center')
else:
    ax2.text(0.5, 0.5, 'DE Data Incomplete', ha='center', va='center')

# Ax3: 样本群构成堆叠图
ax3 = fig.add_subplot(gs[2])
props_summary = adata.obs.groupby(['diagnosis', 'cell_type'], observed=False).size().unstack()
props_summary = props_summary.div(props_summary.sum(axis=1), axis=0)
props_summary.plot(kind='bar', stacked=True, ax=ax3, colormap='tab20')
ax3.set_title("Cohort Composition Fractions")
ax3.set_ylabel("Fractional Cell Layout")
ax3.tick_params(axis='x', rotation=0)
# 改变 4: 强制将长条图的细胞类型图例（Legend）移到画布最右侧外部，彻底避免其覆盖长条图本身
ax3.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), title="Cell Types", frameon=False)

# 强制校准整体排版边界，不截断任何文字
plt.tight_layout()

# 保存
fig.savefig(os.path.join(plot_dir, "MDD_dlPFC_multipanel.pdf"), bbox_inches='tight', dpi=300)
plt.close(fig)
print("\n>>> Process Pipeline Finalized and Complete.")