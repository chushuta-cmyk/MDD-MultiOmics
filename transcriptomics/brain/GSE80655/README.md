# GSE80655: Postmortem Brain RNA-seq Analysis in MDD

This directory contains the complete analysis pipeline for processing the **GSE80655** RNA-seq dataset, focusing on molecular dysregulation in the **Anterior Cingulate Cortex (AnCg)**, **Nucleus Accumbens (nAcc)**, and **Dorsolateral Prefrontal Cortex (DLPFC)** of subjects with Major Depressive Disorder (MDD).

---

## Dataset Overview

| Attribute | Details |
|-----------|---------|
| **Dataset ID** | [GSE80655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655) |
| **Platform** | Illumina HiSeq 2000 (RNA-seq) |
| **Tissue Type** | Postmortem brain: AnCg, nAcc, DLPFC |
| **Experimental Groups** | MDD vs. Matched Controls (139 total samples, 69 MDD / 70 Control) |
| **Reference** | Ramaker RC et al. *Genome Med*. 2017;9:72. [doi:10.1186/s13073-017-0458-5](https://doi.org/10.1186/s13073-017-0458-5) |

---

## Pipeline Workflow

| Step | Script | Input | Output |
|:---:|:---|:---|:---|
| 1 | `01_loader.R` | Raw count matrix (`GSE80655_raw_counts_GRCh38.p13_NCBI.tsv`) | Deduplicated count matrix: `GSE80655_raw_counts_dedup.rds` |
| 2 | `02_preprocess.R` | Deduplicated counts + GEO phenotype | Integer count matrix (`GSE80655_expr_counts_int.rds`), metadata (`GSE80655_metadata_clean.rds`), region list (`GSE80655_regions_list.rds`), QC plots (PCA + correlation heatmap), Cohen's d per region |
| 3 | `03_normalization.R` | Integer counts + metadata | VST‑transformed matrix (`GSE80655_expression_vst.rds`), feature annotation (`GSE80655_feature_annotation.rds`), post‑normalization PCA |
| 4 | `04_quality_control.R` | Integer counts, metadata, region list, VST matrix | Merged DEG table (`GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv`), GO/KEGG/Reactome enrichment (`*_Enrichment.csv`), volcano plot, MA plot, enrichment plots |

---

## Key Findings

| Region | DEGs (P < 0.01) | Mean \|Cohen's d\| | Top Enriched Pathways |
|--------|-----------------|-------------------|------------------------|
| **AnCg** | 620 | 0.2544 | Oxidative phosphorylation, PI3K-Akt, Mitochondrial translation |
| **nAcc** | 372 | 0.2459 | PI3K-Akt, FoxO, Elastic fibre formation |
| **DLPFC** | 293 | 0.2323 | PI3K-Akt, Translation, Elastic fibre formation |

**Biological Interpretation**:
- AnCg shows the strongest transcriptional signal, consistent with GSE54563 validation.
- Enrichment of **mitochondrial translation** and **oxidative phosphorylation** suggests energy metabolism disruption.
- **Extracellular matrix / elastic fibre** pathways indicate structural remodeling in MDD.

---

## Output Files

| File | Description |
|------|-------------|
| `GSE80655_raw_counts_dedup.rds` | Deduplicated raw count matrix |
| `GSE80655_expr_counts_int.rds` | Integer count matrix (DESeq2-ready) |
| `GSE80655_metadata_clean.rds` | Cleaned sample metadata |
| `GSE80655_regions_list.rds` | GSM IDs per brain region |
| `GSE80655_expression_vst.rds` | VST‑transformed matrix (MOFA-ready) |
| `GSE80655_feature_annotation.rds` | Gene annotation file |
| `GSE80655_MDD_vs_Control_BrainRegion_DEGs.csv` | Merged DEG table (all regions, P < 0.01) |
| `GSE80655_KEGG_Enrichment.csv` | KEGG pathway enrichment results |
| `GSE80655_Reactome_Enrichment.csv` | Reactome pathway enrichment results |
| `QC_Sample_Correlation_Heatmap.tiff` / `.pdf` | Sample correlation heatmap |
| `QC_PCA.tiff` / `.pdf` | PCA plot (pre‑VST) |
| `QC_PCA_VST.tiff` / `.pdf` | PCA plot (post‑VST) |
| `Volcano_GSE80655.tiff` / `.pdf` | Volcano plot |
| `MA_GSE80655.tiff` / `.pdf` | MA plot |
| `GSE80655_KEGG_Barplot.tiff` / `.pdf` | KEGG pathway barplot |
| `GSE80655_Reactome_Dotplot.tiff` / `.pdf` | Reactome pathway dotplot |

---

## Usage

### Prerequisites
- R (≥ 4.3)
- Required packages: `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`, `ggplot2`, `pheatmap`, `dplyr`, `tibble`, `GEOquery`

Install missing packages:
```r
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db",
                       "ReactomePA", "ggplot2", "pheatmap", "dplyr",
                       "tibble", "GEOquery"))
```

### Run the Pipeline
Place `GSE80655_raw_counts_GRCh38.p13_NCBI.tsv` in the working directory and run scripts in order:
```bash
Rscript 01_loader.R
Rscript 02_preprocess.R
Rscript 03_normalization.R
Rscript 04_quality_control.R
```

---

## Citation & Data Source

If you use this pipeline or the processed results in your research, please cite:

- **Original study**: Ramaker RC, Bowling KM, Lasseigne BN, et al. Post‑mortem molecular profiling of three psychiatric disorders. *Genome Med*. 2017;9:72. [doi:10.1186/s13073-017-0458-5](https://doi.org/10.1186/s13073-017-0458-5)
- **NCBI GEO Accession**: [GSE80655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80655)