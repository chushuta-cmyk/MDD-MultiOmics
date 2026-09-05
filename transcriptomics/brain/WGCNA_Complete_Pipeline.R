# ==============================================================================
# WGCNA_Complete_Pipeline.R
# Complete Weighted Gene Co-expression Network Analysis (WGCNA) pipeline
# for three brain regions (hippocampus, associative striatum, PFC BA46)
# from GSE53987 and GSE80655, following standard protocols.
# ==============================================================================

cat("\n")
cat(strrep("=", 80), "\n")
cat("🧬 COMPLETE WGCNA PIPELINE: Multi-region MDD Transcriptomics\n")
cat(strrep("=", 80), "\n\n")

# ------------------------------------------------------------------------------
# 0. Load required packages and set up environment
# ------------------------------------------------------------------------------
cat("STEP 0: Loading required packages...\n")

required_packages <- c(
  "WGCNA", "dplyr", "tidyr", "ggplot2", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "ReactomePA", "RColorBrewer",
  "tibble", "stringr"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "ReactomePA", "WGCNA")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg, quiet = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Enable multithreading (adjust to your machine)
enableWGCNAThreads(nThreads = 6)

# Set directories
input_dir <- "./WGCNA_Input_Files/"
output_dir <- "./WGCNA_Results/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define regions (must match file names in input_dir)
regions <- c("hippocampus", "Associative_striatum", "Pre-frontal_cortex_BA46")

cat("✓ Packages loaded.\n")
cat("   Regions:", paste(regions, collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# 1. Load expression and clinical data for each region
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 1: Loading expression and clinical data\n")
cat(strrep("=", 80), "\n\n")

expr_list <- list()
trait_list <- list()

for (region in regions) {
  expr_file <- paste0(input_dir, region, "_expression.csv")
  clin_file <- paste0(input_dir, region, "_clinical.csv")
  
  if (!file.exists(expr_file) || !file.exists(clin_file)) {
    stop(paste0("Missing input files for region: ", region))
  }
  
  # Load expression (rows = genes, columns = samples)
  expr_data <- read.csv(expr_file, row.names = 1, stringsAsFactors = FALSE)
  # Transpose to samples x genes for WGCNA
  datExpr <- as.matrix(t(expr_data))
  storage.mode(datExpr) <- "numeric"
  
  # Load clinical traits
  datTraits <- read.csv(clin_file, row.names = 1, stringsAsFactors = FALSE)
  storage.mode(datTraits$MDD_Diagnosis) <- "numeric"
  
  # Match samples
  sample_names <- rownames(datExpr)
  datTraits <- datTraits[match(sample_names, rownames(datTraits)), , drop = FALSE]
  
  # Quality check
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) {
    datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
    datTraits <- datTraits[gsg$goodSamples, , drop = FALSE]
    cat(paste0("   ⚠️  ", region, ": removed ", sum(!gsg$goodSamples), " samples and ", sum(!gsg$goodGenes), " genes.\n"))
  }
  
  expr_list[[region]] <- datExpr
  trait_list[[region]] <- datTraits
  
  cat(paste0("   ✅ ", region, ": ", nrow(datExpr), " samples, ", ncol(datExpr), " genes.\n"))
}
cat("\n")

# ------------------------------------------------------------------------------
# 2. Soft-thresholding power selection for each region
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 2: Soft-thresholding power selection\n")
cat(strrep("=", 80), "\n\n")

power_results <- list()
power_plots <- list()

for (region in regions) {
  datExpr <- expr_list[[region]]
  
  # Choose a set of powers
  powers <- c(1:20)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
  
  fit_indices <- sft$fitIndices
  # Pick the first power that yields R² > 0.8; fallback to 12
  suggested_power <- fit_indices$Power[which(fit_indices$SFT.R.sq > 0.8)][1]
  if (is.na(suggested_power)) suggested_power <- 12
  
  power_results[[region]] <- list(
    suggested = suggested_power,
    fit_indices = fit_indices
  )
  
  cat(paste0("   ✅ ", region, ": suggested power = ", suggested_power, "\n"))
  
  # Save plot
  p <- ggplot(fit_indices, aes(x = Power, y = SFT.R.sq)) +
    geom_point(size = 2) +
    geom_line() +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    labs(title = paste0(region, ": Scale-free topology fit"),
         x = "Soft-threshold power", y = "SFT R²") +
    theme_minimal()
  power_plots[[region]] <- p
}

# Save combined plot
pdf(file.path(output_dir, "Power_Selection_All_Regions.pdf"), width = 12, height = 8)
for (region in regions) {
  print(power_plots[[region]])
}
dev.off()
cat("   ✅ Power selection plots saved.\n\n")

# Use a unified power for consistency (we'll use the maximum suggested or a fixed value)
# Here we use 12 as in your original code, but you could use region-specific.
power_used <- 12
cat(paste0("   Using unified power = ", power_used, " for all regions.\n\n"))

# ------------------------------------------------------------------------------
# 3. Network construction and module detection (per region)
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 3: Network construction and module detection\n")
cat(strrep("=", 80), "\n\n")

wgcna_results <- list()

for (region in regions) {
  cat(paste0("--- ", region, " ---\n"))
  datExpr <- expr_list[[region]]
  datTraits <- trait_list[[region]]
  
  # Build network using blockwiseModules
  net <- blockwiseModules(
    datExpr,
    power = power_used,
    maxBlockSize = 5000,
    TOMType = "unsigned",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE,
    verbose = 0
  )
  
  # Extract module eigengenes and order them
  MEs <- net$MEs
  MEs_renamed <- orderMEs(MEs)
  
  # Module-trait correlations
  moduleTraitCor <- cor(MEs_renamed, datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  # Identify MDD-associated modules (p < 0.05)
  mdd_modules <- rownames(moduleTraitCor)[moduleTraitPvalue[, 1] < 0.05]
  
  cat(paste0("   Total modules: ", length(unique(net$colors)), "\n"))
  cat(paste0("   MDD-associated modules: ", length(mdd_modules), "\n"))
  
  # Extract hub genes (kME > 0.7) for each MDD module
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs_renamed, use = "p"))
  hub_list <- list()
  
  for (mod in mdd_modules) {
    color <- gsub("ME", "", mod)
    hub_indices <- which(geneModuleMembership[[mod]] > 0.7)
    hub_genes <- rownames(geneModuleMembership)[hub_indices]
    
    if (length(hub_genes) > 0) {
      hub_df <- data.frame(
        GeneSymbol = hub_genes,
        kME = geneModuleMembership[hub_genes, mod],
        Direction = moduleTraitCor[mod, 1],  # correlation with MDD
        Region = region,
        Module_Color = color
      ) %>% arrange(desc(kME))
      hub_list[[color]] <- hub_df
    }
  }
  
  if (length(hub_list) > 0) {
    region_hubs <- bind_rows(hub_list)
    cat(paste0("   Hub genes (kME > 0.7): ", nrow(region_hubs), "\n"))
  } else {
    region_hubs <- data.frame(
      GeneSymbol = character(0), kME = numeric(0),
      Direction = numeric(0), Region = character(0), Module_Color = character(0)
    )
    cat("   No hub genes found (kME > 0.7).\n")
  }
  
  # Save per-region hub genes
  write.csv(region_hubs,
            file.path(output_dir, paste0(region, "_Hub_Genes.csv")),
            row.names = FALSE)
  
  # Store results
  wgcna_results[[region]] <- list(
    net = net,
    MEs = MEs_renamed,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    mdd_modules = mdd_modules,
    hub_genes = region_hubs,
    datExpr = datExpr,
    datTraits = datTraits
  )
  
  cat("\n")
}

# Combine all hub genes
all_hubs <- bind_rows(lapply(wgcna_results, function(x) x$hub_genes))
write.csv(all_hubs,
          file.path(output_dir, "All_Regions_Hub_Genes.csv"),
          row.names = FALSE)
cat("✅ Hub genes combined: ", nrow(all_hubs), " total.\n\n")

# ------------------------------------------------------------------------------
# 4. Module-trait correlation heatmaps
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 4: Module-trait correlation heatmaps\n")
cat(strrep("=", 80), "\n\n")

for (region in regions) {
  res <- wgcna_results[[region]]
  if (is.null(res) || nrow(res$MEs) == 0) next
  
  pdf(file.path(output_dir, paste0(region, "_Module_Trait_Heatmap.pdf")),
      width = 8, height = 10)
  plotEigengeneNetworks(
    res$MEs,
    setLabels = region,
    marDendro = c(0, 4, 1, 2),
    marHeatmap = c(3, 4, 1, 2),
    cex.lab = 0.8
  )
  dev.off()
}
cat("✅ Module-trait heatmaps saved.\n\n")

# ------------------------------------------------------------------------------
# 5. Cross-region module preservation analysis
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 5: Cross-region module preservation\n")
cat(strrep("=", 80), "\n\n")

# Use the first region as reference (e.g., associative striatum)
ref_region <- "Associative_striatum"
ref_res <- wgcna_results[[ref_region]]
ref_network <- list(data = t(ref_res$datExpr))
ref_colors <- ref_res$net$colors

# Test against other regions
test_regions <- setdiff(regions, ref_region)
preservation_results <- list()

for (test in test_regions) {
  test_res <- wgcna_results[[test]]
  test_network <- list(data = t(test_res$datExpr))
  test_colors <- test_res$net$colors
  
  cat(paste0("   ", ref_region, " vs ", test, "\n"))
  
  mp <- modulePreservation(
    multiData = list(Ref = ref_network, Test = test_network),
    multiColor = list(Ref = ref_colors, Test = test_colors),
    referenceNetworks = 1,
    nPermutations = 200,
    verbose = 0
  )
  
  mp_stats <- mp$preservation$Z$Ref$Test[, c("moduleSize", "Zsummary.pres")]
  mp_stats$Module <- rownames(mp_stats)
  preservation_results[[test]] <- mp_stats
  
  write.csv(mp_stats,
            file.path(output_dir, paste0("Preservation_", ref_region, "_vs_", test, ".csv")),
            row.names = FALSE)
}
cat("✅ Preservation analysis completed.\n\n")

# ------------------------------------------------------------------------------
# 6. Functional enrichment of hub genes (GO and Reactome)
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 6: Functional enrichment of hub genes\n")
cat(strrep("=", 80), "\n\n")

# Use all unique hub genes
hub_symbols <- unique(all_hubs$GeneSymbol)
if (length(hub_symbols) > 10) {
  # Map to Entrez IDs
  entrez_ids <- suppressMessages(
    mapIds(org.Hs.eg.db,
           keys = hub_symbols,
           column = "ENTREZID",
           keytype = "SYMBOL",
           multiVals = "first")
  )
  entrez_clean <- entrez_ids[!is.na(entrez_ids)]
  
  # GO BP enrichment
  go_res <- enrichGO(
    gene = entrez_clean,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    go_df <- as.data.frame(go_res@result)
    write.csv(go_df,
              file.path(output_dir, "HubGenes_GO_BP_Enrichment.csv"),
              row.names = FALSE)
    
    # Dotplot
    p_go <- dotplot(go_res, showCategory = 15)
    ggsave(file.path(output_dir, "HubGenes_GO_BP_Dotplot.pdf"),
           p_go, width = 10, height = 8)
    cat("✅ GO enrichment completed.\n")
  }
  
  # Reactome enrichment
  reactome_res <- enrichPathway(
    gene = entrez_clean,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  if (!is.null(reactome_res) && nrow(reactome_res@result) > 0) {
    reactome_df <- as.data.frame(reactome_res@result)
    write.csv(reactome_df,
              file.path(output_dir, "HubGenes_Reactome_Enrichment.csv"),
              row.names = FALSE)
    cat("✅ Reactome enrichment completed.\n")
  }
} else {
  cat("⚠️  Not enough hub genes for enrichment analysis.\n")
}
cat("\n")

# ------------------------------------------------------------------------------
# 7. Cross-region hub gene overlap analysis
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 7: Cross-region hub gene overlap\n")
cat(strrep("=", 80), "\n\n")

region_gene_lists <- lapply(regions, function(r) {
  res <- wgcna_results[[r]]
  if (!is.null(res) && nrow(res$hub_genes) > 0) {
    return(unique(res$hub_genes$GeneSymbol))
  } else {
    return(character(0))
  }
})
names(region_gene_lists) <- regions

# Pairwise overlaps
pairwise_overlap <- data.frame()
for (i in 1:length(regions)) {
  for (j in (i+1):length(regions)) {
    if (j <= length(regions)) {
      overlap <- intersect(region_gene_lists[[regions[i]]],
                           region_gene_lists[[regions[j]]])
      pairwise_overlap <- rbind(pairwise_overlap,
                                data.frame(
                                  Region1 = regions[i],
                                  Region2 = regions[j],
                                  Overlap = length(overlap),
                                  Pct1 = round(100 * length(overlap) / length(region_gene_lists[[regions[i]]]), 1),
                                  Pct2 = round(100 * length(overlap) / length(region_gene_lists[[regions[j]]]), 1)
                                ))
    }
  }
}

write.csv(pairwise_overlap,
          file.path(output_dir, "HubGene_CrossRegion_Overlap.csv"),
          row.names = FALSE)

# Common to all three
common_all <- Reduce(intersect, region_gene_lists)
cat(paste0("   Genes shared across all 3 regions: ", length(common_all), "\n"))
if (length(common_all) > 0) {
  cat("   ", paste(common_all, collapse = ", "), "\n")
}
cat("\n")

# ------------------------------------------------------------------------------
# 8. Generate summary report
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("STEP 8: Generating summary report\n")
cat(strrep("=", 80), "\n\n")

sink(file.path(output_dir, "WGCNA_Analysis_Report.txt"))

cat("========================================\n")
cat("WGCNA Analysis Report\n")
cat("Date:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("========================================\n\n")

cat("Regions analyzed:\n")
for (r in regions) cat(paste0("  - ", r, "\n"))
cat("\n")

cat("WGCNA parameters:\n")
cat(paste0("  - Power: ", power_used, "\n"))
cat("  - Min module size: 30\n")
cat("  - Merge cut height: 0.25\n")
cat("  - Hub gene kME cutoff: 0.7\n\n")

cat("Summary per region:\n")
for (region in regions) {
  res <- wgcna_results[[region]]
  if (!is.null(res)) {
    cat(paste0("\n--- ", region, " ---\n"))
    cat(paste0("  Total modules: ", length(unique(res$net$colors)), "\n"))
    cat(paste0("  MDD-associated modules: ", length(res$mdd_modules), "\n"))
    cat(paste0("  Hub genes: ", nrow(res$hub_genes), "\n"))
    if (nrow(res$hub_genes) > 0) {
      cat("  Top 5 hub genes:\n")
      top5 <- res$hub_genes %>% arrange(desc(kME)) %>% head(5)
      for (i in 1:nrow(top5)) {
        cat(paste0("    ", top5$GeneSymbol[i], " (kME = ", round(top5$kME[i], 3), ")\n"))
      }
    }
  }
}
cat("\n")

cat("Cross-region preservation (Z-summary):\n")
for (test in test_regions) {
  mp <- preservation_results[[test]]
  if (!is.null(mp)) {
    cat(paste0("\n  ", ref_region, " vs ", test, ":\n"))
    mp_sig <- mp[mp$Zsummary.pres > 2, ]
    if (nrow(mp_sig) > 0) {
      cat(paste0("    Modules with Z > 2: ", nrow(mp_sig), "\n"))
      for (i in 1:min(5, nrow(mp_sig))) {
        cat(paste0("      Module ", mp_sig$Module[i], ": Z = ", round(mp_sig$Zsummary.pres[i], 2), "\n"))
      }
    } else {
      cat("    No modules with Z > 2.\n")
    }
  }
}
cat("\n")

cat("Cross-region hub gene overlap:\n")
if (length(common_all) > 0) {
  cat(paste0("  Genes shared by all 3 regions: ", length(common_all), "\n"))
  cat(paste0("  ", paste(common_all, collapse = ", "), "\n"))
} else {
  cat("  No genes shared by all 3 regions.\n")
}
cat("\n")

cat("Output files:\n")
cat(paste0("  - ", output_dir, "*_Hub_Genes.csv\n"))
cat(paste0("  - ", output_dir, "All_Regions_Hub_Genes.csv\n"))
cat(paste0("  - ", output_dir, "*_Module_Trait_Heatmap.pdf\n"))
cat(paste0("  - ", output_dir, "Preservation_*.csv\n"))
cat(paste0("  - ", output_dir, "HubGenes_*_Enrichment.csv\n"))
cat(paste0("  - ", output_dir, "HubGene_CrossRegion_Overlap.csv\n"))
cat("\n========================================\n")
cat("Analysis complete.\n")
cat("========================================\n")

sink()

cat("✅ Summary report saved.\n\n")

# ------------------------------------------------------------------------------
# 9. Final message
# ------------------------------------------------------------------------------
cat(strrep("=", 80), "\n")
cat("✅ WGCNA COMPLETE PIPELINE FINISHED SUCCESSFULLY\n")
cat(strrep("=", 80), "\n\n")
cat("📁 Results saved in:", output_dir, "\n")
cat("📄 See WGCNA_Analysis_Report.txt for summary.\n")
cat("📊 Hub genes are in All_Regions_Hub_Genes.csv.\n\n")
cat(strrep("=", 80), "\n")