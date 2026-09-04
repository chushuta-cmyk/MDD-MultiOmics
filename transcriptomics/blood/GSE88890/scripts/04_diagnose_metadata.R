# ============================================================================
# GSE88890 Metadata Diagnostic Tool
# Checks series_matrix.txt to verify line numbers and metadata extraction
# ============================================================================

cat("\n=== GSE88890 Metadata Diagnostic Tool ===\n\n")

# ----------------------------------------------------------------------------
# 1. Locate Sample_title line
# ----------------------------------------------------------------------------
cat("--- Searching for Sample_title line ---\n\n")

series_matrix_file <- "GSE88890_series_matrix.txt"
all_lines <- readLines(series_matrix_file)

for (i in 1:70) {
  if (grepl("Sample_title", all_lines[i])) {
    cat(sprintf("Line %d: %s\n", i, substr(all_lines[i], 1, 150)))
  }
}

# ----------------------------------------------------------------------------
# 2. Extract and validate sample names
# ----------------------------------------------------------------------------
cat("\n--- Extracting sample names ---\n")

title_line <- all_lines[grep("!Sample_title", all_lines)]

if (length(title_line) > 0) {
  title_parts <- unlist(strsplit(title_line[1], "\t"))
  sample_names <- gsub('"', '', title_parts[-1])
  
  cat(sprintf("✓ Found %d samples\n", length(sample_names)))
  cat("\nFirst 20 sample names:\n")
  print(sample_names[1:min(20, length(sample_names))])
  
  # Check for duplicates
  cat("\nChecking for duplicates...\n")
  dup_names <- sample_names[duplicated(sample_names)]
  if (length(dup_names) > 0) {
    cat(sprintf("⚠️ Found %d duplicate names:\n", length(dup_names)))
    print(table(dup_names))
  } else {
    cat("✓ No duplicates found\n")
  }
  
  # Check for unexpected content
  cat("\nChecking sample name formats...\n")
  if (any(grepl("tissue|Brodmann|gender", sample_names))) {
    cat("⚠️ WARNING: Sample names contain unexpected fields!\n")
    bad_idx <- which(grepl("tissue|Brodmann|gender", sample_names))
    print(sample_names[bad_idx])
  } else {
    cat("✓ Sample names format OK\n")
  }
} else {
  cat("✗ ERROR: Cannot find !Sample_title line\n")
}

# ----------------------------------------------------------------------------
# 3. Extract all characteristics lines
# ----------------------------------------------------------------------------
cat("\n--- Extracting characteristics lines ---\n")

all_char <- grep("^!Sample_characteristics_ch1", all_lines, value = TRUE)
cat(sprintf("✓ Found %d characteristics lines\n\n", length(all_char)))

cat("First 5 characteristics lines:\n")
for (i in 1:min(5, length(all_char))) {
  content <- substr(all_char[i], 1, 100)
  cat(sprintf("  [%d] %s...\n", i, content))
}

# ----------------------------------------------------------------------------
# 4. Extract individual features
# ----------------------------------------------------------------------------
cat("\n--- Extracting features ---\n")

# groupid
groupid_line <- all_char[grep("groupid", all_char)]
if (length(groupid_line) > 0) {
  groupid_parts <- unlist(strsplit(groupid_line[1], "\t"))[-1]
  groupid_values <- gsub('"', '', groupid_parts)
  groupid_values <- gsub('.*groupid: ', '', groupid_values)
  
  cat(sprintf("✓ groupid: %d values\n", length(groupid_values)))
  cat(sprintf("  Values: %s\n", paste(unique(groupid_values), collapse = ", ")))
}

# Region / Brodmann area
region_line <- all_char[grep("tissue|Brodmann", all_char)]
if (length(region_line) > 0) {
  region_parts <- unlist(strsplit(region_line[1], "\t"))[-1]
  region_values <- gsub('"', '', region_parts)
  brodmann_numbers <- gsub('.*Brodmann area ([0-9]+).*', '\\1', region_values)
  
  cat(sprintf("✓ region: %d values\n", length(region_values)))
  cat(sprintf("  Values: %s\n", paste(unique(brodmann_numbers), collapse = ", ")))
}

# Gender
gender_line <- all_char[grep("gender", all_char)]
if (length(gender_line) > 0) {
  gender_parts <- unlist(strsplit(gender_line[1], "\t"))[-1]
  gender_values <- gsub('"', '', gender_parts)
  gender_values <- gsub('.*gender: ', '', gender_values)
  
  cat(sprintf("✓ gender: %d values\n", length(gender_values)))
  cat(sprintf("  Values: %s\n", paste(unique(gender_values), collapse = ", ")))
}

# ----------------------------------------------------------------------------
# 5. Consistency check
# ----------------------------------------------------------------------------
cat("\n--- Data consistency check ---\n")

cat(sprintf("Sample names count: %d\n", length(sample_names)))
cat(sprintf("groupid count: %d\n", length(groupid_values)))
cat(sprintf("Region count: %d\n", length(region_values)))
cat(sprintf("Gender count: %d\n", length(gender_values)))

if (length(sample_names) == length(groupid_values) && 
    length(sample_names) == length(region_values) &&
    length(sample_names) == length(gender_values)) {
  cat("\n✓ All counts match!\n")
} else {
  cat("\n⚠️ Count mismatch! Check data extraction.\n")
}

# ----------------------------------------------------------------------------
# 6. Build and display corrected targets table
# ----------------------------------------------------------------------------
cat("\n--- Creating corrected targets table ---\n")

sample_id <- sprintf("Sample_%03d", 1:length(sample_names))

phenotype_values <- ifelse(
  grepl("MDD suicide", groupid_values),
  "MDD_suicide",
  "Control"
)

region_clean <- gsub('.*Brodmann area ([0-9]+).*', 'BA\\1', region_values)

targets <- data.frame(
  Sample_ID = sample_id,
  Sample_Name = sample_names,
  Phenotype = as.factor(phenotype_values),
  Region = as.factor(region_clean),
  Gender = as.factor(gender_values),
  stringsAsFactors = FALSE
)

cat("✓ targets table created\n\n")
cat("Targets summary:\n")
print(head(targets, 10))

cat("\nSample distribution:\n")
print(table(targets$Phenotype, targets$Region))

cat("\n✓ Diagnostic check completed!\n")

# ----------------------------------------------------------------------------
# 7. Fix code snippet (reusable solution)
# ----------------------------------------------------------------------------
cat("\n--- Fix code snippet (copy this into main script) ---\n\n")

cat("
# Correct metadata extraction code

series_matrix_file <- 'GSE88890_series_matrix.txt'
all_lines <- readLines(series_matrix_file)

# Extract sample names
title_line <- all_lines[grep('^!Sample_title', all_lines)]
title_parts <- unlist(strsplit(title_line, '\\t'))
sample_names <- gsub('\"', '', title_parts[-1])

# Extract all characteristics
all_char <- grep('^!Sample_characteristics_ch1', all_lines, value = TRUE)

# groupid
groupid_line <- all_char[grep('groupid', all_char)]
groupid_raw <- gsub('.*groupid: |\"', '', unlist(strsplit(groupid_line, '\\t'))[-1])
phenotype <- ifelse(grepl('MDD suicide', groupid_raw), 'MDD_suicide', 'Control')

# Region
region_line <- all_char[grep('tissue|Brodmann', all_char)]
region_raw <- gsub('\"', '', unlist(strsplit(region_line, '\\t'))[-1])
region <- gsub('.*Brodmann area ([0-9]+).*', 'BA\\\\1', region_raw)

# Gender
gender_line <- all_char[grep('gender', all_char)]
gender <- gsub('.*gender: |\"', '', unlist(strsplit(gender_line, '\\t'))[-1])

# Build targets table
targets <- data.frame(
  Sample_ID = sprintf('Sample_%03d', 1:length(sample_names)),
  Sample_Name = sample_names,
  Phenotype = as.factor(phenotype),
  Region = as.factor(region),
  Gender = as.factor(gender),
  stringsAsFactors = FALSE
)
rownames(targets) <- targets$Sample_ID

# Verify
print(table(targets$Phenotype, targets$Region))
")