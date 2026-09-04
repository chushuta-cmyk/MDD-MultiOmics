# GSE54563 — Postmortem Brain Transcriptomic Analysis (Microarray)

## 1. GEO Information

| Field | Content |
|-------|---------|
| GEO Accession | GSE54563 |
| Title | Gene expression in anterior cingulate cortex in major depressive disorder |
| Species | Homo sapiens |
| Platform | Illumina HumanHT-12 V3.0 expression beadchip (GPL6947) |
| Technology | Microarray |
| Reference | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54563 |
| Original Publication | Chang LC et al. PLoS One. 2014;9(3):e90980. doi:10.1371/journal.pone.0090980 |

## 2. Sample Design

- **Total samples**: 50
- **MDD**: 25
- **Control**: 25
- **Design**: MDD‑Control matched pairs
- **Brain region**: Anterior Cingulate Cortex (ACC)
- **Probe count**: 48,803

## 3. Role in This Project

- **Role**: Independent validation dataset
- **Manuscript sections**: §2.2 (Validation), §3.1 (Gene‑level validation), Table 3, Figure 2C, Supplementary Table S3
- **Key results** (from updated analysis):
  - `P < 0.05`: 294 genes (after deduplication)
  - `P < 0.01`: ~50 genes
  - `P < 0.001`: 4 genes
  - Mean logFC: near 0; median logFC: −0.0177
  - Mean Cohen's d: 0.177; median Cohen's d: 0.149
  - **Functional enrichment** (GO BP) highlights **myelination**, **oligodendrocyte differentiation**, **axon ensheathment**, and **dendritic spine development**. KEGG and Reactome showed no significant enrichment with the current gene set.

## 4. Analysis Pipeline (Integrated Script)

The pipeline is implemented in a single R script (`GSE54563_Complete_Pipeline.R`) with the following steps:

1. **Data loading**  
   - Reads expression matrix and phenotype from `GSE54563_series_matrix.txt` using `data.table::fread()`.

2. **Auto log2 detection**  
   - Checks `max(expr_clean)`. If > 50, applies `log2(x + 1)`; otherwise skips.

3. **Quality control**  
   - Removes probes with any NA.
   - Keeps probes with expression > 3 (log2 scale) in at least 10% of samples.

4. **Differential expression (paired Limma)**  
   - Builds design matrix: `~ group + pair` (pair as factor).
   - Fits linear model with `lmFit()`, contrasts via `makeContrasts()`, and `eBayes()`.
   - Extracts results with `topTable()` using BH adjustment.

5. **Gene annotation & deduplication**  
   - Maps `PROBEID` to `SYMBOL` using `illuminaHumanv3.db`.
   - Calculates average expression per probe (`rowMeans` on filtered matrix).
   - Sorts by `AvgExpr` descending and removes duplicate `GeneSymbol` entries (keeps highest expressed probe).

6. **Effect size (Cohen's d)**  
   - Computes pooled‑standard‑deviation Cohen's d for each gene.

7. **Output of DEG results**  
   - Saves full results and a summary table.

8. **Functional enrichment**  
   - Selects genes with `P.Value < 0.05` for enrichment.
   - Converts `SYMBOL` → `ENTREZID` using `org.Hs.eg.db::bitr()`.
   - **GO BP**: `enrichGO()` (pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE).
   - **KEGG**: `enrichKEGG()` (pvalueCutoff=0.1, qvalueCutoff=0.2, use_internal_data=TRUE if `KEGG.db` available; otherwise online with timeout=600s).
   - **Reactome**: `enrichPathway()` from `ReactomePA` (pvalueCutoff=0.05, readable=TRUE).

9. **Visualization**  
   - Generates dotplots (GO, Reactome) and barplot (KEGG) using `enrichplot`.
   - Saves each as both TIFF (300 dpi, LZW compression) and PDF via custom `save_sci_plot()` function.

## 5. Input Files (Required)

| File | Description |
|------|-------------|
| `GSE54563_series_matrix.txt` | Raw GEO series matrix file (must be in working directory). |

## 6. Output Files (Generated)

| File | Description |
|------|-------------|
| `GSE54563_Top16_DEGs.csv` | Top 16 genes by P‑value (probe_id, logFC, cohens_d, P.Value, adj.P.Val). |
| `GSE54563_All_DEG_Results.csv` | Full DEG table (probe_id, logFC, cohens_d, P.Value, adj.P.Val). |
| `GSE54563_Analysis_Summary.csv` | Summary statistics (sample counts, probes, significant DEG counts, mean Cohen's d). |
| `GO_BP_Enrichment.csv` | GO Biological Process enrichment results (if significant). |
| `GO_BP_Dotplot.tiff` / `.pdf` | Dotplot of top 15 GO terms (TIFF + PDF). |
| `KEGG_Enrichment.csv` | KEGG enrichment results (if significant; may be empty). |
| `KEGG_Barplot.tiff` / `.pdf` | Barplot of top 15 KEGG pathways (if significant). |
| `Reactome_Enrichment.csv` | Reactome enrichment results (if significant; may be empty). |
| `Reactome_Dotplot.tiff` / `.pdf` | Dotplot of top 15 Reactome pathways (if significant). |

## 7. Main Functions / Scripts

| File / Function | Purpose |
|-----------------|---------|
| `GSE54563_Complete_Pipeline.R` | **Integrated single script** (all steps from loading to enrichment). |
| `save_sci_plot(p, f, w=8, h=6)` | Custom function to save ggplot objects as TIFF (with LZW compression) and PDF simultaneously. |
| `enrichGO()` (clusterProfiler) | GO BP enrichment. |
| `enrichKEGG()` (clusterProfiler) | KEGG pathway enrichment (uses local `KEGG.db` if installed, otherwise online). |
| `enrichPathway()` (ReactomePA) | Reactome pathway enrichment. |
| `dotplot()`, `barplot()` (enrichplot) | Visualization of enrichment results. |

## 8. Notes

- **KEGG and Reactome** may yield no significant results for this dataset – this is a biological finding (genes enriched in GO processes may not map to classic pathways).
- The script automatically handles log2 transformation and low‑expression filtering.
- All plots are saved in publication‑ready TIFF format (300 dpi, LZW compression) plus PDF for easy editing.

## 9. Citation & Data Source

If you use this pipeline or the processed results in your research, please cite:

* **Original publication**: French, L., et al. (2014). Transcriptome analysis of the human anterior cingulate cortex in major depression. *Translational Psychiatry*, 4(12), e474. [doi:10.1038/tp.2014.111](https://doi.org/10.1038/tp.2014.111)
* **NCBI GEO Accession**: Edgar, R., Domrachev, M., & Lash, A. E. (2002). Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. *Nucleic Acids Research*, 30(1), 207‑210.