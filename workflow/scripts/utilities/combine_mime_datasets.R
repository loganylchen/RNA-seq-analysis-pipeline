#!/usr/bin/env Rscript
# Combine individual Mime dataset RDS files into one list for Mime analysis
# Determines training dataset from dataset_type == "discovery" in samples.tsv

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

cat("==============================================================\n")
cat("Combine Mime Datasets\n")
cat("==============================================================\n\n")

# Get parameters from Snakemake
samples_file <- snakemake@params[["samples"]]
project <- snakemake@params[["project"]]

# Get individual RDS files from input
rds_files <- snakemake@input[["rds_files"]]
output_file <- snakemake@output[["combined_rds"]]

cat("Parameters:\n")
cat("  Samples file:", samples_file, "\n")
cat("  Project:", project, "\n")
cat("  Number of RDS files:", length(rds_files), "\n")
cat("  Output file:", output_file, "\n\n")

# Read samples metadata
cat("Reading samples metadata...\n")
samples_df <- read.delim(samples_file, header=TRUE, stringsAsFactors = FALSE)
samples_df <- samples_df[samples_df$project_id == project, ]
cat("  Loaded", nrow(samples_df), "samples for project", project, "\n")

# Show unique dataset_types and dataset_ids
if ("dataset_type" %in% colnames(samples_df)) {
    cat("  Available dataset_types:", paste(unique(samples_df$dataset_type), collapse=", "), "\n")
} else {
    cat("  WARNING: dataset_type column not found in samples file\n")
}
if ("dataset_id" %in% colnames(samples_df)) {
    cat("  Available dataset_ids:", paste(unique(samples_df$dataset_id), collapse=", "), "\n")
} else {
    cat("  WARNING: dataset_id column not found in samples file\n")
}
cat("\n")

# Function to extract dataset name from RDS file path
extract_dataset_name <- function(rds_path) {
    # Extract filename without extension and path
    filename <- basename(rds_path)
    # Remove "_response_dataset.rds" suffix
    dataset_name <- sub("_response_dataset\\.rds$", "", filename)
    return(dataset_name)
}

# Load all individual RDS files
cat("\nLoading individual RDS files...\n")
all_datasets <- list()

for (i in seq_along(rds_files)) {
    rds_file <- rds_files[i]
    dataset_name <- extract_dataset_name(rds_file)

    cat("  Loading:", rds_file, "\n")
    cat("    Dataset name:", dataset_name, "\n")

    # Load the RDS file
    dataset_data <- readRDS(rds_file)

    # Show what we loaded
    if (is.data.frame(dataset_data)) {
        cat("    Loaded data.frame with", nrow(dataset_data), "samples x", ncol(dataset_data), "columns\n")
        cat("    Columns:", paste(colnames(dataset_data), collapse=", "), "\n")
    } else if (is.list(dataset_data)) {
        cat("    Loaded list with", length(dataset_data), "elements\n")
        cat("    Element names:", paste(names(dataset_data), collapse=", "), "\n")
    }

    all_datasets[[dataset_name]] <- dataset_data
}

cat("\nLoaded", length(all_datasets), "datasets:\n")
for (name in names(all_datasets)) {
    cat("  ", name, "\n")
}
cat("\n")

# Determine training and validation datasets
# Filter samples where dataset_type == "discovery" and extract unique dataset_id
training_dataset_ids <- samples_df %>%
    filter(dataset_type == "discovery") %>%
    pull(dataset_id) %>%
    unique()

if (length(training_dataset_ids) == 0) {
    stop("ERROR: No samples found with dataset_type 'discovery'")
}

if (length(training_dataset_ids) > 1) {
    stop("ERROR: Multiple dataset_ids found with dataset_type 'discovery': ",
         paste(training_dataset_ids, collapse=", "))
}

training_dataset_id <- training_dataset_ids[1]
cat("  Training dataset_id (from dataset_type='discovery'):", training_dataset_id, "\n")

# Find the matching dataset name in RDS files
# Dataset names in RDS files are the dataset_id values
training_name <- training_dataset_id

if (!training_name %in% names(all_datasets)) {
    stop("ERROR: Training dataset '", training_name, "' not found in RDS files.\n",
         "Available datasets: ", paste(names(all_datasets), collapse=", "))
}

validation_names <- setdiff(names(all_datasets), training_name)

cat("\nTraining dataset (Dataset1):", training_name, "\n")
if (length(validation_names) > 0) {
    cat("Validation datasets:", paste(validation_names, collapse=", "), "\n")
} else {
    cat("No validation datasets\n")
}
cat("\n")

# Create Mime-compatible list
# Dataset1 = training (discovery)
# Dataset2, Dataset3, ... = validation datasets

mime_datasets <- list()

# Add training dataset as Dataset1
cat("\nCreating Mime dataset structure:\n")
cat("  Dataset1 (training):", training_name, "\n")
mime_datasets[["Dataset1"]] <- all_datasets[[training_name]]
cat("    Samples:", nrow(mime_datasets[["Dataset1"]]), "\n")

# Add validation datasets
if (length(validation_names) > 0) {
    for (i in seq_along(validation_names)) {
        dataset_name <- validation_names[i]
        mime_dataset_name <- paste0("Dataset", i + 1)
        cat("  ", mime_dataset_name, " (validation):", dataset_name, "\n")
        mime_datasets[[mime_dataset_name]] <- all_datasets[[dataset_name]]
        cat("    Samples:", nrow(mime_datasets[[mime_dataset_name]]), "\n")
    }
}

# Save combined RDS file
cat("\nSaving combined Mime dataset to:", output_file, "\n")
saveRDS(mime_datasets, file = output_file)

cat("\n==============================================================\n")
cat("Combined Mime dataset created successfully!\n")
cat("==============================================================\n\n")

cat("Output summary:\n")
for (name in names(mime_datasets)) {
    data <- mime_datasets[[name]]
    if (is.data.frame(data)) {
        cat("  ", name, ":", nrow(data), "samples x", ncol(data), "columns\n")
    } else if (is.matrix(data)) {
        cat("  ", name, ":", nrow(data), "samples x", ncol(data), "columns\n")
    } else {
        cat("  ", name, ":", class(data)[1], "\n")
    }
}
cat("\n")

cat("To use in Mime:\n")
cat("  load('", output_file, "')\n", sep="")
cat("  # Access datasets: list_train_vali_Data$Dataset1 (training),\n")
cat("  #                 list_train_vali_Data$Dataset2, etc. (validation)\n\n")

# Close logging
sink()
sink(type="message")
close(log)