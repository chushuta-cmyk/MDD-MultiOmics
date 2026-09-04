# GSE53987: Postmortem Brain Transcriptomic Analysis in MDD

This directory contains the complete analysis pipeline for processing the GSE53987 microarray dataset, evaluating transcriptomic dysregulation across multiple postmortem brain regions in Major Depressive Disorder (MDD). All required scripts and output files are provided.

---

## Dataset Overview

-   **Platform**: Affymetrix Human Genome U133 Plus 2.0 Array (GPL570)
    
-   **Brain Regions**: Hippocampus (HIP), Pre‑frontal Cortex BA46 (DLPFC), Associative Striatum (STR)
    
-   **Groups**: MDD vs. age/sex‑matched Controls
    
-   **Total Samples**: 105 (MDD and Control only; original study used a tetrad design including schizophrenia and bipolar disorder)
    

---

## Pipeline Workflow

The analysis is implemented in `GSE53987_organizedAnalysis.R` and performs the following steps:

1.  **Data Loading**  
    Reads the local series matrix file (`skip = 66`), sets `ID_REF` as row names.
    
2.  **Automatic Log2 Transformation**  
    If raw expression values exceed 50, the matrix is automatically log2‑transformed to meet limma assumptions.
    
3.  **Quantile Normalization**  
    Applies `limma::normalizeBetweenArrays()` to reduce technical variation.
    
4.  **Quality Control**  
    Generates diagnostic plots:
    
    -   Boxplots before/after normalization
        
    -   Density plot of expression distributions
        
5.  **Expression Filtering**
    
    -   Removes probes with any missing values
        
    -   Removes low‑expression probes (rowMeans < median)
        
    -   Removes low‑variance probes (below the 25th percentile)
        
6.  **Phenotype Parsing**  
    Extracts brain region and diagnosis from `!Sample_source_name_ch1`, retaining only MDD and Control samples.
    
7.  **Region‑Stratified Differential Expression**  
    For each brain region:
    
    -   Fits a linear model using `limma`
        
    -   Defines contrast: MDD vs Control
        
    -   Applies empirical Bayes moderation (`eBayes`)
        
    -   Filters results by `P.Value < 0.01` (raw p‑value, not FDR‑adjusted)
        
8.  **Probe‑to‑Gene Annotation**  
    Maps Affymetrix probe IDs to Gene Symbols using `hgu133plus2.db`, with an uppercase fallback for unmatched probes.
    
9.  **Deduplication**  
    For genes with multiple probes, retains the probe with the highest absolute logFC.
    
10.  **Visualization**  
     Generates volcano plots and heatmaps for each region and for the entire dataset (uses `plot_generator.R`).
     
11.  **Output**  
     Merges all regional DEG tables into a single CSV file.
     

---

## Scripts

| Script | Description |
| --- | --- |
| `GSE53987_organizedAnalysis.R` | Main analysis script (loading → normalization → DEA → annotation → visualization) |
| `plot_generator.R` | Helper script for volcano and heatmap generation (required) |

---

## Dependencies

### R Packages

```
library(stringr)
library(GEOquery)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(hgu133plus2.db)
```


### External Script

-   `plot_generator.R` – must be placed in the same directory as the main script.
    

---

## Usage

1.  Place `GSE53987_series_matrix.txt` in a `/data/` directory (relative to the script location) or adjust the `file_path` variable.
    
2.  Ensure `plot_generator.R` is in the same folder.
    
3.  Run the script:

```
Rscript GSE53987_reorganized.R
```    

All outputs will be written to the current working directory (or to `/results/` if you modify the path).

---

## Output Files

| File | Description |
| --- | --- |
| `GSE53987_MDD_vs_Control_BrainRegion_DEGs.csv` | Merged DEG results (GeneSymbol, logFC, P.Value, adj.P.Val, Region) |
| `QC_Boxplot_BeforeNorm.pdf` | Boxplot of raw expression values |
| `QC_Boxplot_AfterNorm.pdf` | Boxplot after quantile normalization |
| `QC_DensityPlot.pdf` | Density plot of normalized expression |
| `volcano_<region>_volcano.pdf/.tiff` | Volcano plots per brain region and for all regions combined |
| `heatmap_<region>_heatmap.pdf/.tiff` | Heatmaps of top 40 DEGs per region (with sample annotations) |
| `session_info.txt` | Full R session information for reproducibility |

> **Note**: Both PDF and high‑resolution TIFF (600 dpi, LZW compression) versions are generated for all figures.

---

## Key Findings

The following numbers of unique genes were identified as differentially expressed (raw P < 0.01) in each brain region:

| Brain Region | DEGs (unique genes) |
| --- | --- |
| Hippocampus | 40 |
| Pre‑frontal Cortex BA46 | 267 |
| Associative Striatum | 206 |

> **Interpretation**: The minimal overlap across regions supports a gene‑level heterogeneity model, where distinct sets of genes are dysregulated in different brain areas, rather than a common core signature.

---

## Citation

If you use this pipeline or the processed results in your research, please cite:

-   **Original Study**: Lanz, T. A., et al. (2019). Postmortem transcriptional profiling reveals widespread increase in inflammation in schizophrenia: a comparison of prefrontal cortex, striatum, and hippocampus among matched tetrads of controls with subjects diagnosed with schizophrenia, bipolar or major depressive disorder. *Translational Psychiatry*, 9(1), 151. [doi:10.1038/s41398-019-0492-8](https://doi.org/10.1038/s41398-019-0492-8)
    
-   **NCBI GEO Accession**: [GSE53987](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987)