#!/usr/bin/env Rscript
# Prepare Mime-compatible datasets from RNA-seq pipeline output
# This script creates datasets for Mime machine learning framework
#
# Output format for Mime:
# - First column: ID (sample ID)
# - Survival analysis: OS.time (survival time), OS (status: 0=censored, 1=event)
# - Response prediction: Var (Y=response, N=no response)
# - Remaining columns: gene expression (log2(x+1) transformed)

# Setup logging
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# Get parameters from Snakemake
samples_file <- snakemake@input[["samples"]]
counts_file <- snakemake@input[["expression"]]
output_file <- snakemake@output[["mime_rds"]]


dataset <- snakemake@params[["dataset"]]
project <- snakemake@params[["project"]]
response_col <- snakemake@params[["response_col"]]

# Ensure response_col is a character string
if (is.list(response_col)) {
    response_col <- unlist(response_col)[1]
}
response_col <- as.character(response_col)


cat("==============================================================\n")
cat("Mime Dataset Preparation\n")
cat("==============================================================\n\n")

cat("Parameters:\n")
cat("  Samples file:", samples_file, "\n")
cat("  Project:", project, "\n")
cat("  Count matrix:", counts_file, "\n")
cat("  Output:", output_file, "\n")




# Read samples metadata
cat("Reading samples metadata...\n")
samples_df <- read.delim(samples_file, header=TRUE, stringsAsFactors = FALSE)
samples_df <- samples_df %>% dplyr::filter(project_id == project, dataset_id == dataset)
cat("  Loaded", nrow(samples_df), "samples for project", project, "\n")

# Read count matrix
cat("Reading count matrix...\n")
counts_df <- read.delim(counts_file, header=TRUE, row.names = 1, check.names = FALSE)
cat("  Dimensions:", nrow(counts_df), "genes x", ncol(counts_df), "samples\n")

# Transpose counts: genes as columns, samples as rows (Mime format)
cat("Transforming count matrix...\n")
expression_df <- as.data.frame(t(counts_df))
expression_df$ID <- rownames(expression_df)
cat("  Transposed to:", nrow(expression_df), "samples x", ncol(expression_df), "columns\n")

# Apply log2(x+1) transformation (standard for RNA-seq)
cat("Applying log2(x+1) transformation...\n")
gene_cols <- setdiff(colnames(expression_df), "ID")
# for (col in gene_cols) {
#     expression_df[[col]] <- log2(expression_df[[col]] + 1)
# }
expression_df[, gene_cols] <- log2(expression_df[, gene_cols] + 1)  
cat("  Transformation applied to", length(gene_cols), "genes\n")

# Filter samples to those in count matrix
expression_df <- expression_df[expression_df$ID %in% samples_df$sample_name, ]
cat("  Matched", nrow(expression_df), "samples between count matrix and metadata\n")


cat("\nPreparing response prediction dataset...\n")
cat("  Response column:", response_col, "\n")
cat("  Response column type:", class(response_col), "\n")

# Validate required columns
if (is.null(response_col) || response_col == "" || is.na(response_col)) {
    stop("ERROR: response_col is required for response prediction")
}

if (!response_col %in% colnames(samples_df)) {
    cat("  Available columns in samples_df:\n")
    print(colnames(samples_df))
    stop("ERROR: Response column '", response_col, "' not found in samples metadata")
}

# Prepare columns to select from samples_df
cols_to_select <- c("sample_name", response_col)
cat("  Columns to select from samples_df:", paste(cols_to_select, collapse = ", "), "\n")

# Create a safe subset of samples_df
samples_subset <- samples_df[, cols_to_select, drop = FALSE]
cat("  samples_subset dimensions:", nrow(samples_subset), "x", ncol(samples_subset), "\n")

# Merge response information
cat("  Merging expression data with response information...\n")
result_df <- merge(
    expression_df,
    samples_subset,
    by.x = "ID",
    by.y = "sample_name",
    all.x = TRUE
)
cat("  Merged result dimensions:", nrow(result_df), "x", ncol(result_df), "\n")

# Rename column to Mime format
colnames(result_df)[colnames(result_df) == response_col] <- "Var"

# Convert to Y/N format
if (is.numeric(result_df$Var)) {
    cat("  Converting numeric response to Y/N format...\n")
    result_df$Var <- ifelse(result_df$Var > 0, "Y", "N")
} else if (is.logical(result_df$Var)) {
    cat("  Converting logical response to Y/N format...\n")
    result_df$Var <- ifelse(result_df$Var, "Y", "N")
} else if (is.character(result_df$Var) || is.factor(result_df$Var)) {
    cat("  Converting character/factor response to Y/N format...\n")
    result_df$Var <- ifelse(toupper(as.character(result_df$Var)) %in% c("Y", "YES", "1", "TRUE","TUMOR"), "Y", "N")
} else {
    stop("ERROR: Unsupported data type for response variable 'Var'")
}

# Remove samples with missing response
n_before <- nrow(result_df)
result_df <- result_df[!is.na(result_df$Var), ]
if (nrow(result_df) < n_before) {
    cat("  Removed", n_before - nrow(result_df), "samples with missing response data\n")
}


cat("\nResponse data summary:\n")
cat("  Response (Y):", sum(result_df$Var == "Y"), "samples\n")
cat("  Non-response (N):", sum(result_df$Var == "N"), "samples\n")
cat("  Response rate:", mean(result_df$Var == "Y"), "\n")


# Save to RDS file
cat("\nSaving Mime-compatible dataset to:", output_file, "\n")
saveRDS(result_df, file = output_file)

cat("\n==============================================================\n")
cat("Dataset preparation complete!\n")
cat("==============================================================\n\n")


# Close logging
sink()
sink(type="message")
close(log)