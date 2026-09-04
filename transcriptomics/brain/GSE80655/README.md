# GSE88890 — Brain DNA Methylation (Exploratory)

This directory contains scripts for analyzing the Illumina 450K DNA methylation dataset [GSE88890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88890), which profiles post‑mortem brain tissue (BA11 and BA25) from MDD suicide cases and non‑psychiatric controls.

> **Note on role in the project:** This dataset is **exploratory**. The methylation signal is weak (only 14 significant DMPs at FDR < 0.05). Therefore, it was **not included in the main analysis** of the paper but is retained as supplementary material.

---

## 📊 Dataset Overview

| Item | Description |
| --- | --- |
| **GEO Accession** | [GSE88890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88890) |
| **Title** | Methylomic profiling of cortex samples from completed suicide cases implicates a role for PSORS1C3 in major depression and suicide |
| **Platform** | Illumina HumanMethylation450 BeadChip |
| **Total samples** | 75 (37 MDD suicide, 38 controls) |
| **Brain regions** | BA11 (40 samples), BA25 (35 samples) |
| **Reference** | Murphy TM et al. *Transl Psychiatry*. 2017;7(1):e989. [doi:10.1038/tp.2016.249](https://doi.org/10.1038/tp.2016.249) |

### Key findings in this dataset

-   Significant DMPs at FDR < 0.05: **14**
    
    -   Hypermethylated: 10 (71.4%)
        
    -   Hypomethylated: 4 (28.6%)
        
-   DMR analysis and pathway enrichment were **not performed** due to the weak signal.
    
-   This dataset is **not part of the main manuscript** but is used for supplementary exploration and as input for MOFA integration.
    

---

## 📁 Retained Scripts (5 Core Scripts)

After reviewing all available versions, we have selected **5 functional scripts** that cover the full workflow (loading, QC, DMP, annotation, visualisation, and optional continuation).

| Script | Original File | Purpose | When to use |
| --- | --- | --- | --- |
| **01_main_analysis.R** | `GSE88890_Complete_Analysis_Final.R` | **End‑to‑end pipeline**   Extracts metadata from `series_matrix.txt` (lines 50–53), runs QC → DMP (limma) → annotation → volcano/heatmap → report. | **First run** (recommended). Works locally or on server. |
| **02_improved_analysis.R** | `GSE88890_Improved_Complete_Analysis.R` | **Enhanced version** with SNP‑probe removal, covariate adjustment (Region + Gender), and a lightweight DMR scan. | If you need stricter covariate control or want to explore DMRs. |
| **03_resume_analysis.R** | `GSE88890_With_GEO2R_Groups_CHAMP.fixed.R` | **Resume / patch‑run** script.   Checks existing files in `./data` and `./results`, skips completed steps, and only runs missing DMR or pathway enrichment. | After you have DMP results and want to add DMR or GO/KEGG without recomputing all probes. |
| **04_diagnose_metadata.R** | `GSE88890_Diagnostic_Script.R` | **Metadata diagnostic tool**   Prints all sample‑related lines from `series_matrix.txt` to verify line numbers (50, 52, 53). | **Only if 01_main_analysis.R fails** due to missing groupid or region info. |
| **05_next_steps_guide.R** | `GSE88890_Next_Steps_Guide.R` | **Planning guide** (no computation)   Provides code suggestions for downstream analyses (DMR, GO/KEGG, brain‑blood integration, MOFA export). | Read after completing the main analysis to plan next steps. |

---

## 📂 Recommended Directory Structure

```
GSE88890/
├── data/
│   ├── GSE88890\_series\_matrix.txt          # required
│   └── GSE88890\_normalisedBetas.csv        # required
├── scripts/
│   ├── 01\_main\_analysis.R
│   ├── 02\_improved\_analysis.R
│   ├── 03\_resume\_analysis.R
│   ├── diagnostics/
│   │   └── 04\_diagnose\_metadata.R
│   └── guides/
│       └── 05\_next\_steps\_guide.R
├── results/                                 # auto‑generated
│   ├── 01\_DMP\_Global\_Results.csv
│   ├── 02\_DMP\_Annotated\_Significant.csv
│   ├── 03\_VolcanoPlot.pdf
│   ├── 04\_Top30\_DMP\_Heatmap.pdf
│   └── FINAL\_REPORT.txt
└── README.md                               # this file
```


> **Note:** If you plan to export gene‑level beta matrices for MOFA, use the standalone script `GSE88890_standalone_MOFA_v2.R` (located in `transcriptomics/GSE88890/scripts/`), which is not part of this core set.

---

## 🔧 Dependencies

Key R packages: `data.table`, `limma`, `ggplot2`, `pheatmap`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`.  
Optional for advanced steps: `bumphunter`, `missMethyl`, `clusterProfiler`, `org.Hs.eg.db`.

Install Bioconductor packages:

```
BiocManager::install(c("limma", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       "bumphunter", "missMethyl", "clusterProfiler", "org.Hs.eg.db"))
```

---

## 🚀 Recommended Execution Order

### 1\. Run main analysis

```
cd /path/to/GSE88890
Rscript scripts/01\_main\_analysis.R
```

If you encounter errors (e.g., “groupid not found”), run:


```
Rscript scripts/diagnostics/04\_diagnose\_metadata.R
```

Then adjust the line numbers (`50`, `52`, `53`) in `01_main_analysis.R` and re‑run.

### 2\. (Optional) Improved version


```
Rscript scripts/02\_improved\_analysis.R
```

### 3\. (Optional) Add DMR / pathway enrichment


```
Rscript scripts/03\_resume\_analysis.R
```

### 4\. Review results

-   `results/FINAL_REPORT.txt` – full statistics.
    
-   `results/03_VolcanoPlot.pdf` and `04_Top30_DMP_Heatmap.pdf` – visualisation.
    

---

## 📤 Output Files

| File | Description |
| --- | --- |
| `01_DMP_Global_Results.csv` | Full DMP table (all CpGs). |
| `02_DMP_Annotated_Significant.csv` | Annotated significant DMPs (FDR < 0.05). |
| `03_VolcanoPlot.pdf` | Volcano plot highlighting significant DMPs. |
| `04_Top30_DMP_Heatmap.pdf` | Heatmap of top 30 DMPs. |
| `FINAL_REPORT.txt` | Complete summary report. |
| *(optional)* `05_DMR_Results.csv` | DMR results (if run). |
| *(optional)* `06_GO_Enrichment.csv` | GO enrichment (if run). |

---

## ⚠️ Important Notes

-   **Weak signal**: This dataset produced only 14 significant DMPs; it is **exploratory** and **not used as primary evidence** in the manuscript.
    
-   **Memory**: 450K arrays require **8–16 GB RAM**.
    
-   **Line numbers**: Verify `series_matrix.txt` line numbers with the diagnostic script if parsing fails.
    
-   **MOFA export**: If you need gene‑level beta matrices for MOFA integration, use the dedicated script `GSE88890_standalone_MOFA_v2.R` (not included here).
    

---

## 📖 Citation

If you use this dataset or pipeline, please cite:

-   **GEO accession**: [GSE88890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE88890)
    
-   **Original study**: Murphy TM et al. (2017). *Transl Psychiatry* 7(1):e989. [doi:10.1038/tp.2016.249](https://doi.org/10.1038/tp.2016.249)
    

---

*Last updated: 2026‑09‑04*