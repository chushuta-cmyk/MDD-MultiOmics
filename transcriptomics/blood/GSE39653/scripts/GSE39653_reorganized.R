# ============================================================================
# GSE39653 Blood B Cell Validation - Clean Pipeline
# Objective: Compare B cell abundance (MCP-counter) between MDD and Control
#           and validate against brain tissue findings.
# ============================================================================

cat("\n=== GSE39653 Blood B Cell Validation Started ===\n")

# ----------------------------------------------------------------------------
# 1. Setup
# ----------------------------------------------------------------------------
library(immunedeconv)
library(ggplot2)

# ----------------------------------------------------------------------------
# 2. Load data from local file
# ----------------------------------------------------------------------------
file_path <- "datasets/GSE39653_series_matrix.txt"
all_lines <- readLines(file_path)

# 2.1 Phenotype (lines 32-62)
pheno_lines <- all_lines[32:62]
pheno_split <- strsplit(pheno_lines, "\t")
pheno_names <- sub("^!Sample_", "", sapply(pheno_split, `[`, 1))
pheno_matrix <- do.call(rbind, lapply(pheno_split, `[`, -1))
meta <- as.data.frame(t(pheno_matrix), stringsAsFactors = FALSE)
colnames(meta) <- pheno_names
meta[] <- lapply(meta, function(x) gsub('"', '', x))

# 2.2 Expression matrix (lines 64-43180)
mat_lines <- all_lines[64:43180]
tmp <- tempfile()
writeLines(mat_lines, tmp)
expr_raw <- read.table(tmp, header = TRUE, sep = "\t", check.names = FALSE)
unlink(tmp)

rownames(expr_raw) <- expr_raw[[1]]
expr <- as.matrix(expr_raw[, -1])
storage.mode(expr) <- "numeric"

cat(sprintf("Data loaded: %d probes x %d samples\n", nrow(expr), ncol(expr)))

# ----------------------------------------------------------------------------
# 3. Sample grouping (MDD vs Control, exclude BD)
# ----------------------------------------------------------------------------
meta$Diagnosis <- ifelse(
  grepl("^MDD", meta$title), "MDD",
  ifelse(grepl("^HC", meta$title), "Control", NA)
)

keep <- !is.na(meta$Diagnosis)
meta <- meta[keep, ]
expr <- expr[, keep]

cat(sprintf("Samples: MDD=%d, Control=%d\n",
            sum(meta$Diagnosis == "MDD"),
            sum(meta$Diagnosis == "Control")))

# ----------------------------------------------------------------------------
# 4. Probe → Gene Symbol mapping (GPL10558 annotation)
# ----------------------------------------------------------------------------
gpl_file <- "datasets/GPL10558_annotation.txt"
if (!file.exists(gpl_file)) stop("GPL annotation file not found")

gpl <- read.delim(gpl_file, stringsAsFactors = FALSE)

# Auto-detect column names
probe_col <- intersect(c("ID", "ID_REF", "probe_id", "ProbeID"), colnames(gpl))[1]
symbol_col <- intersect(c("SYMBOL", "Gene_Symbol", "gene_symbol", "Symbol"), colnames(gpl))[1]

if (is.na(probe_col) || is.na(symbol_col)) stop("Column detection failed")

# Build mapping and aggregate duplicates
probe_map <- gpl[!is.na(gpl[[symbol_col]]) & gpl[[symbol_col]] != "", c(probe_col, symbol_col)]
colnames(probe_map) <- c("ProbeID", "Symbol")
probe_map <- probe_map[!duplicated(probe_map$ProbeID), ]

common <- intersect(rownames(expr), probe_map$ProbeID)
expr_mapped <- expr[common, ]
rownames(expr_mapped) <- probe_map[match(common, probe_map$ProbeID), "Symbol"]

# Average duplicate genes
expr_final <- aggregate(expr_mapped, by = list(Gene = rownames(expr_mapped)), FUN = mean)
rownames(expr_final) <- expr_final$Gene
expr_final$Gene <- NULL
expr_final <- as.matrix(expr_final)

cat(sprintf("Final expression: %d genes x %d samples\n", nrow(expr_final), ncol(expr_final)))

# ----------------------------------------------------------------------------
# 5. MCP-counter immune deconvolution
# ----------------------------------------------------------------------------
immune_scores <- deconvolute(expr_final, method = "mcp_counter")
rownames(immune_scores) <- immune_scores$cell_type
immune_scores$cell_type <- NULL
immune_scores <- as.matrix(immune_scores)

cat(sprintf("Immune scores: %d cell types x %d samples\n",
            nrow(immune_scores), ncol(immune_scores)))

# ----------------------------------------------------------------------------
# 6. B cell comparison: MDD vs Control
# ----------------------------------------------------------------------------
# Align sample order
meta_ordered <- meta[match(colnames(immune_scores), meta$geo_accession), ]

# Find B cell row
bcell_row <- rownames(immune_scores)[grepl("B.cell", rownames(immune_scores), ignore.case = TRUE)][1]
if (is.na(bcell_row)) stop("B cell row not found")

bcell <- as.numeric(immune_scores[bcell_row, ])
mdd_b <- bcell[meta_ordered$Diagnosis == "MDD"]
ctrl_b <- bcell[meta_ordered$Diagnosis == "Control"]

# Statistical test
wt <- wilcox.test(mdd_b, ctrl_b)
fc <- mean(mdd_b) / mean(ctrl_b)

cat(sprintf("\nB cell: FC=%.3f, p=%.4f\n", fc, wt$p.value))

# ----------------------------------------------------------------------------
# 7. All cell types comparison
# ----------------------------------------------------------------------------
all_res <- data.frame()
for (ct in rownames(immune_scores)) {
  s <- as.numeric(immune_scores[ct, ])
  m <- s[meta_ordered$Diagnosis == "MDD"]
  c <- s[meta_ordered$Diagnosis == "Control"]
  w <- wilcox.test(m, c)
  all_res <- rbind(all_res, data.frame(
    CellType = ct, FC = mean(m)/mean(c), Pvalue = w$p.value
  ))
}
all_res$FDR <- p.adjust(all_res$Pvalue, method = "BH")
all_res <- all_res[order(all_res$Pvalue), ]

# ----------------------------------------------------------------------------
# 8. Compare with brain tissue results
# ----------------------------------------------------------------------------
brain_fc <- 0.759
brain_p <- 0.019

consistency <- (fc < 1 & brain_fc < 1)  # both decreased
cat(sprintf("Brain tissue: FC=%.3f, p=%.4f\n", brain_fc, brain_p))
cat(sprintf("Consistency (both decreased): %s\n", consistency))

# ----------------------------------------------------------------------------
# 9. Visualization
# ----------------------------------------------------------------------------
plot_df <- data.frame(
  Sample = colnames(immune_scores),
  BCells = bcell,
  Diagnosis = meta_ordered$Diagnosis
)

p <- ggplot(plot_df, aes(x = Diagnosis, y = BCells, fill = Diagnosis)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c("Control" = "#377EB8", "MDD" = "#E41A1C")) +
  labs(title = "B Cell Abundance: MDD vs Control",
       subtitle = paste0("FC = ", round(fc, 3), ", p = ", signif(wt$p.value, 4)),
       y = "MCP-counter Score") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("GSE39653_BCell_Boxplot.pdf", p, width = 7, height = 6)

# ----------------------------------------------------------------------------
# 10. Export results
# ----------------------------------------------------------------------------
write.csv(immune_scores, "GSE39653_Immune_Scores.csv")
write.csv(all_res, "GSE39653_All_CellTypes.csv", row.names = FALSE)

cat("\nOutputs saved.\n")
