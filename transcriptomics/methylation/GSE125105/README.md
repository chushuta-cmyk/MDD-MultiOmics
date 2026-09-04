# GSE125105 — Peripheral Blood DNA Methylation Analysis

This pipeline processes Illumina HumanMethylation450 BeadChip data from GSE125105 (MDD vs control whole blood) to identify differentially methylated probes (DMPs), annotate them to genes, perform functional enrichment, and integrate with brain expression results.

## 1. GEO Information

- **Accession**: [GSE125105](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125105)
- **Title**: DeepWAS: multivariate genotype-phenotype associations by directly integrating regulatory information using deep learning
- **Platform**: Illumina HumanMethylation450 BeadChip
- **Tissue**: Whole blood
- **Groups**: MDD vs Control
- **Reference**: Arloth J et al., *PLoS Comput Biol* 2020;16(2):e1007616. [doi:10.1371/journal.pcbi.1007616](https://doi.org/10.1371/journal.pcbi.1007616)

## 2. Sample Design

- Tissue: Whole blood
- Groups: MDD vs Control
- Methylation‑associated genes (mapped from DMPs): **4,488**

## 3. Role in This Project

- **Chapter**: §2.6 Peripheral blood DNA methylation analysis, §3.7 Epigenetic regulation preferentially targets network hub genes, Figure 7, Supplementary Table S8–S9
- **Key findings**:
  - Methylation genes: 4,488
  - Overlap with convergent brain DEGs (≥3 regions, n=10): **only 1 gene (*COX19*)**
  - Overlap with WGCNA hub genes (n=616): **125 genes**
  - Selected hub genes: *SYNE1*, *NTRK2*, *STAT1*, *ITGB1*, *PICALM*, *PRKCB*, *DNAH10*, *RAB2A*, *FLOT1*, *PRKAR2A*, *AGBL2*, *CCDC33*
  - **Core conclusion**: Epigenetic regulation preferentially targets network hub genes rather than consistently differentially expressed genes.

## 4. Pipeline Overview

1. **Data loading & QC**  
   `01_load_and_qc.R` / `01_load_and_qc_original.R` – load methylation matrix, QC (missingness, distribution plots, MDS), save RDS objects.

2. **DMP analysis**  
   `02_dmp_analysis.R` – chunked limma (10,000 probes per chunk) to avoid memory overflow.  
   `02_dmp_analysis_simple.R` – quick test on first 10,000 probes.  
   `GSE125105_DMP_Limma_Only_fixed.R` – standalone limma fallback.  
   Threshold: adjusted p‑value < 0.05.

3. **Annotation**  
   `03_annotation.R` – map probes to genes, chromosome, CpG islands, gene regions (using IlluminaHumanMethylation450kanno.ilmn12.hg19).

4. **Visualization**  
   `04_visualization.R` – volcano plot, heatmap, Manhattan plot, distribution plots (all PDF).

5. **Enrichment & report**  
   `05_enrichment_and_report.R` – GO (BP/CC/MF), KEGG, Disease Ontology; compiles Markdown/HTML report.

6. **Cross‑omics overlap**  
   `futher_analysis/03_methylation_overlap/03_methylation_overlap.R` – overlaps methylation genes with convergent DEGs and WGCNA hub genes.

## 5. Main Scripts

| Script | Location | Purpose |
|--------|----------|---------|
| `01_load_and_qc.R` / `_original.R` | `complete_projects/GSE125105/pipeline/` | Data loading and QC |
| `02_dmp_analysis.R` / `_simple.R` | `complete_projects/GSE125105/pipeline/` | DMP analysis |
| `GSE125105_DMP_Limma_Only_fixed.R` | `complete_projects/GSE125105/pipeline/` | Standalone limma (fallback) |
| `03_annotation.R` | `complete_projects/GSE125105/pipeline/` | Probe→gene annotation |
| `04_visualization.R` | `complete_projects/GSE125105/pipeline/` | Visualisation |
| `05_enrichment_and_report.R` | `complete_projects/GSE125105/pipeline/` | Enrichment and report |
| `GSE125105__GSE125105_CHAMP_*.R` | `complete_projects/GSE125105/scripts/` | CHAMP workflow scripts |
| `03_methylation_overlap.R` | `futher_analysis/03_methylation_overlap/` | Cross‑omics overlap |

## 6. Key Outputs

- `GSE125105_DMP_all_results.csv` – full DMP table  
- `GSE125105_DMP_significant.csv` – DMPs with FDR < 0.05  
- `GSE125105_DMP_annotated_all.csv` – annotated full table  
- `GSE125105_DMP_annotated_significant.csv` – annotated significant DMPs  
- `methylation_gene_list.txt` – unique gene symbols (4,488)  
- `methylation_overlap_convergent.csv` – overlap with convergent DEGs (only *COX19*)  
- `methylation_overlap_hub.csv` – overlap with WGCNA hubs (125 genes)  
- `methylation_overlap_summary.txt` – overlap statistics  
- `volcano_plot_enhanced.pdf`, `top_dmp_heatmap.pdf` – figures  
- `GO_BP_enrichment.csv`, `KEGG_enrichment.csv` – enrichment results

## 7. Dependencies

Main R packages: `data.table`, `limma`, `ggplot2`, `reshape2`, `dplyr`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`, `clusterProfiler`, `org.Hs.eg.db`, `DOSE`, `pheatmap`, `ggrepel`, `RColorBrewer`.

Install Bioconductor packages:
```r
BiocManager::install(c("limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       "clusterProfiler", "org.Hs.eg.db", "DOSE"))

## 8. Usage

1. Place data files in `DATA_DIR` (defined in each script):  
   - `GSE125105_series_matrix.txt.gz`  
   - `GSE125105_matrix_normalized.txt.gz`  
   - `GSE125105_CHAMP_targets.csv` (sample ID and group)

2. Run the pipeline sequentially:

```bash
Rscript 01_load_and_qc.R
Rscript 02_dmp_analysis.R
Rscript 03_annotation.R
Rscript 04_visualization.R
Rscript 05_enrichment_and_report.R
```

3. After completion, run the overlap analysis:

```bash
cd futher_analysis/03_methylation_overlap/
Rscript 03_methylation_overlap.R
```

## 9. Main Findings

- **4,488** methylation‑associated genes (mapped from significant DMPs).  
- Overlap with **convergent brain DEGs** (≥3 regions): only **1 gene** (*COX19*).  
- Overlap with **WGCNA hub genes** (n=616): **125 genes**, including important hubs like *SYNE1*, *NTRK2*, *STAT1*, *ITGB1*, *PICALM*, *PRKCB*, *DNAH10*, etc.  
- **Conclusion**: Epigenetic regulation preferentially targets network hub genes rather than genes that show consistent expression changes across brain regions. This supports a model where DNA methylation shapes the core regulatory architecture of MDD.

## 10. Notes

- All scripts are designed to run on a Linux server with sufficient memory (16 GB+ recommended).  
- The `simple` version (`02_dmp_analysis_simple.R`) is for testing only and does not process the full dataset.  
- For integration with other omics (e.g., MOFA), methylation gene‑level data are provided in `next_steps/MOFA/inputs/` 
