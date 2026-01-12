#!/usr/bin/env Rscript

# Comprehensive QC Summary from MultiQC
# Reads MultiQC parsed data and generates a summary table with sample metadata

# Load required libraries
suppressPackageStartupMessages({
    library(jsonlite)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(yaml)
})

# Logging
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat("==============================================================\n")
cat("QC Summary from MultiQC Data\n")
cat("==============================================================\n\n")

# Get input parameters
samples_file <- snakemake@params$samples
project <- snakemake@params$project
multiqc_dir <- snakemake@input$multiqc_dir

# Get output paths
output_summary <- snakemake@output$summary

cat("Project:", project, "\n")
cat("Samples file:", samples_file, "\n")
cat("MultiQC directory:", multiqc_dir, "\n")
cat("Output summary:", output_summary, "\n\n")

# Load sample metadata
cat("Loading sample metadata...\n")
samples_df <- read.delim(samples_file, comment.char = "#", stringsAsFactors = FALSE)
samples_df <- samples_df[samples_df$project_id == project, ]
cat("  Loaded", nrow(samples_df), "samples\n")
cat("  Columns:", paste(colnames(samples_df), collapse = ", "), "\n\n")

# MultiQC data directory
multiqc_data_dir <- file.path(multiqc_dir, "multiqc_data")
cat("Looking for MultiQC data in:", multiqc_data_dir, "\n\n")

# Check if multiqc_data directory exists
if (!dir.exists(multiqc_data_dir)) {
    cat("ERROR: MultiQC data directory not found:", multiqc_data_dir, "\n")
    cat("Please ensure MultiQC has been run first.\n")
    quit(status = 1)
}

# List all files in multiqc_data directory
multiqc_files <- list.files(multiqc_data_dir, pattern = "\\.txt$|\\.json$", full.names = TRUE)
cat("Found", length(multiqc_files), "MultiQC data files:\n")
for (f in multiqc_files) {
    cat("  ", basename(f), "\n")
}
cat("\n")

# Initialize list to store QC metrics for all samples
qc_metrics_list <- list()

# Find all sample names from samples_df
sample_names <- samples_df$sample_name

cat("Processing samples from MultiQC data...\n")
cat(paste0(rep("-", 78), collapse = ""), "\n\n")

# Function to safely read TSV files
safe_read_tsv <- function(file_path) {
    tryCatch({
        data <- read.delim(file_path, stringsAsFactors = FALSE, header = TRUE,
                          comment.char = "#", check.names = FALSE)
        return(data)
    }, error = function(e) {
        return(NULL)
    })
}

# Function to safely read JSON files
safe_read_json <- function(file_path) {
    tryCatch({
        data <- fromJSON(file_path)
        return(data)
    }, error = function(e) {
        return(NULL)
    })
}

# Read MultiQC general stats (contains most metrics)
general_stats_file <- file.path(multiqc_data_dir, "multiqc_general_stats.txt")
if (file.exists(general_stats_file)) {
    cat("Reading MultiQC general stats...\n")
    general_stats <- safe_read_tsv(general_stats_file)
    if (!is.null(general_stats)) {
        cat("  Columns:", paste(colnames(general_stats), collapse = ", "), "\n")
        cat("  Rows:", nrow(general_stats), "\n\n")
    }
}

# Read MultiQC general stats JSON for more detailed metrics
general_stats_json_file <- file.path(multiqc_data_dir, "multiqc_general_stats.json")
general_stats_json <- NULL
if (file.exists(general_stats_json_file)) {
    cat("Reading MultiQC general stats JSON...\n")
    general_stats_json <- safe_read_json(general_stats_json_file)
    if (!is.null(general_stats_json)) {
        cat("  Loaded JSON data\n\n")
    }
}

# Process each sample
for (sample_name in sample_names) {
    cat("Processing sample:", sample_name, "\n")

    # Initialize metrics for this sample
    sample_metrics <- list(
        sample_name = sample_name,
        project_id = project
    )

    # Add metadata from samples.tsv
    meta_cols <- c("dataset_id", "dataset_type", "condition", "patient",
                   "batch", "seq_type", "raw_data")
    for (col in meta_cols) {
        if (col %in% colnames(samples_df)) {
            val <- samples_df[[col]][samples_df$sample_name == sample_name]
            sample_metrics[[col]] <- if (length(val) > 0) val[1] else NA
        }
    }

    # Add info columns dynamically
    info_cols <- grep("^info_", colnames(samples_df), value = TRUE)
    for (col in info_cols) {
        val <- samples_df[[col]][samples_df$sample_name == sample_name]
        sample_metrics[[col]] <- if (length(val) > 0) val[1] else NA
    }

    # === Read from MultiQC general stats ===
    if (!is.null(general_stats) && "Sample" %in% colnames(general_stats)) {
        # Find the row for this sample
        sample_row_idx <- which(general_stats$Sample == sample_name)

        if (length(sample_row_idx) > 0) {
            sample_row <- general_stats[sample_row_idx[1], ]

            # Add all available metrics from general stats
            # Exclude fastp plot data columns (containing "-plot-" in name)
            for (col in colnames(sample_row)) {
                if (col != "Sample") {
                    # Skip fastp plot data columns
                    if (grepl("plot", col, ignore.case = TRUE)) {
                        next
                    }
                    # Convert column name to a more readable format if needed
                    # MultiQC column names are already descriptive
                    metric_name <- gsub(" ", "_", col)
                    sample_metrics[[metric_name]] <- sample_row[[col]]
                }
            }
        }
    }

    # === Read additional metrics from individual MultiQC module files ===

    # fastp data
    fastp_file <- file.path(multiqc_data_dir, "multiqc_fastp.txt")
    # if (file.exists(fastp_file)) {
    #     fastp_data <- safe_read_tsv(fastp_file)
    #     if (!is.null(fastp_data) && "Sample" %in% colnames(fastp_data)) {
    #         sample_row_idx <- which(fastp_data$Sample == sample_name)
    #         if (length(sample_row_idx) > 0) {
    #             sample_row <- fastp_data[sample_row_idx[1], ]
    #             for (col in colnames(sample_row)) {
    #                 if (col != "Sample") {
    #                     sample_metrics[[paste0("fastp_", col)]] <- sample_row[[col]]
    #                 }
    #             }
    #         }
    #     }
    # }

    # STAR data
    star_file <- file.path(multiqc_data_dir, "multiqc_star.txt")
    if (file.exists(star_file)) {
        star_data <- safe_read_tsv(star_file)
        if (!is.null(star_data) && "Sample" %in% colnames(star_data)) {
            sample_row_idx <- which(star_data$Sample == sample_name)
            if (length(sample_row_idx) > 0) {
                sample_row <- star_data[sample_row_idx[1], ]
                for (col in colnames(sample_row)) {
                    if (col != "Sample") {
                        sample_metrics[[paste0("star_", col)]] <- sample_row[[col]]
                    }
                }
            }
        }
    }

    # Qualimap RNA-seq data
    qualimap_file <- file.path(multiqc_data_dir, "multiqc_qualimap_rnaseq.txt")
    if (file.exists(qualimap_file)) {
        qualimap_data <- safe_read_tsv(qualimap_file)
        if (!is.null(qualimap_data) && "Sample" %in% colnames(qualimap_data)) {
            sample_row_idx <- which(qualimap_data$Sample == sample_name)
            if (length(sample_row_idx) > 0) {
                sample_row <- qualimap_data[sample_row_idx[1], ]
                for (col in colnames(sample_row)) {
                    if (col != "Sample") {
                        sample_metrics[[paste0("qualimap_", col)]] <- sample_row[[col]]
                    }
                }
            }
        }
    }

    # Picard RNA-seq data
    picard_rnaseq_file <- file.path(multiqc_data_dir, "multiqc_picard_RNASeq.txt")
    if (file.exists(picard_rnaseq_file)) {
        picard_data <- safe_read_tsv(picard_rnaseq_file)
        if (!is.null(picard_data) && "Sample" %in% colnames(picard_data)) {
            sample_row_idx <- which(picard_data$Sample == sample_name)
            if (length(sample_row_idx) > 0) {
                sample_row <- picard_data[sample_row_idx[1], ]
                for (col in colnames(sample_row)) {
                    if (col != "Sample") {
                        sample_metrics[[paste0("picard_", col)]] <- sample_row[[col]]
                    }
                }
            }
        }
    }

    # Picard insert size data
    picard_insert_file <- file.path(multiqc_data_dir, "multiqc_picard_insertSize.txt")
    if (file.exists(picard_insert_file)) {
        picard_insert_data <- safe_read_tsv(picard_insert_file)
        if (!is.null(picard_insert_data) && "Sample" %in% colnames(picard_insert_data)) {
            sample_row_idx <- which(picard_insert_data$Sample == sample_name)
            if (length(sample_row_idx) > 0) {
                sample_row <- picard_insert_data[sample_row_idx[1], ]
                for (col in colnames(sample_row)) {
                    if (col != "Sample") {
                        sample_metrics[[paste0("picard_insert_", col)]] <- sample_row[[col]]
                    }
                }
            }
        }
    }

    # Picard GC bias data
    picard_gc_file <- file.path(multiqc_data_dir, "multiqc_picard_GC_bias.txt")
    if (file.exists(picard_gc_file)) {
        picard_gc_data <- safe_read_tsv(picard_gc_file)
        if (!is.null(picard_gc_data) && "Sample" %in% colnames(picard_gc_data)) {
            sample_row_idx <- which(picard_gc_data$Sample == sample_name)
            if (length(sample_row_idx) > 0) {
                sample_row <- picard_gc_data[sample_row_idx[1], ]
                for (col in colnames(sample_row)) {
                    if (col != "Sample") {
                        sample_metrics[[paste0("picard_gc_", col)]] <- sample_row[[col]]
                    }
                }
            }
        }
    }

    # RSeQC data
    rseqc_file <- file.path(multiqc_data_dir, "multiqc_rseqc.txt")
    if (file.exists(rseqc_file)) {
        rseqc_data <- safe_read_tsv(rseqc_file)
        if (!is.null(rseqc_data) && "Sample" %in% colnames(rseqc_data)) {
            sample_row_idx <- which(rseqc_data$Sample == sample_name)
            if (length(sample_row_idx) > 0) {
                sample_row <- rseqc_data[sample_row_idx[1], ]
                for (col in colnames(sample_row)) {
                    if (col != "Sample") {
                        sample_metrics[[paste0("rseqc_", col)]] <- sample_row[[col]]
                    }
                }
            }
        }
    }

    # Add to list
    qc_metrics_list[[sample_name]] <- sample_metrics
}

# Convert list to data frame
cat("\n")
cat(paste0(rep("=", 78), collapse = ""), "\n")
cat("Creating summary table...\n")

qc_summary_df <- data.table::rbindlist(
    lapply(qc_metrics_list, function(x) as.data.frame(x, stringsAsFactors = FALSE)),
    fill = TRUE
)

# Reorder columns
# Put sample metadata first
meta_cols_keep <- c("sample_name", "project_id", "dataset_id", "dataset_type",
                    "condition", "patient", "batch", "seq_type", "raw_data")
info_cols_keep <- grep("^info_", colnames(qc_summary_df), value = TRUE)
qc_cols_keep <- setdiff(colnames(qc_summary_df), c(meta_cols_keep, info_cols_keep))

col_order <- c(meta_cols_keep[meta_cols_keep %in% colnames(qc_summary_df)],
               info_cols_keep,
               qc_cols_keep)

qc_summary_df <- qc_summary_df[, col_order, with = FALSE]

# Write output
cat("Writing QC summary to:", output_summary, "\n")
write.table(qc_summary_df, output_summary, sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary statistics
cat("\n")
cat(paste0(rep("=", 78), collapse = ""), "\n")
cat("QC Summary Statistics\n")
cat(paste0(rep("=", 78), collapse = ""), "\n\n")

cat("Total samples processed:", nrow(qc_summary_df), "\n")
cat("Total metrics collected:", ncol(qc_summary_df), "\n\n")

# Show data completeness
cat("Metric completeness:\n")
metric_completeness <- sapply(qc_summary_df, function(x) sum(!is.na(x)) / length(x) * 100)
completeness_df <- data.frame(
    Metric = names(metric_completeness),
    Completeness = round(metric_completeness, 1),
    stringsAsFactors = FALSE
)
completeness_df <- completeness_df[order(-completeness_df$Completeness), ]
colnames(completeness_df) <- c("Metric", "Completeness_%")
print(completeness_df, row.names = FALSE)

# Condition distribution
if ("condition" %in% colnames(qc_summary_df)) {
    cat("\nCondition distribution:\n")
    print(table(qc_summary_df$condition))
}

cat("\n")
cat(paste0(rep("=", 78), collapse = ""), "\n")
cat("QC Summary Complete!\n")
cat(paste0(rep("=", 78), collapse = ""), "\n")
cat("\nOutput files:\n")
cat("  Summary table:", output_summary, "\n")

# Close sinks
sink()
sink(type = "message")
