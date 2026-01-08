#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
    library(jsonlite)
    library(data.table)
    library(dplyr)
    library(tidyr)
})

# Logging
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat("=" , rep("=", 78), "\n", sep = "")
cat("QC Summary Aggregation\n")
cat("=" , rep("=", 78), "\n\n", sep = "")

# Get input parameters
samples_file <- snakemake@params$samples
project <- snakemake@params$project
output_summary <- snakemake@output$summary

cat("Project:", project, "\n")
cat("Samples file:", samples_file, "\n")
cat("Output file:", output_summary, "\n\n")

# Load samples metadata
cat("Loading sample metadata...\n")
samples_df <- read.delim(samples_file, comment.char = "#", stringsAsFactors = FALSE)
samples_df <- samples_df[samples_df$project_id == project, ]
cat("  Loaded", nrow(samples_df), "samples\n\n")

# Initialize list to store QC metrics for all samples
qc_metrics_list <- list()

# Process each sample
for (i in 1:nrow(samples_df)) {
    sample_name <- samples_df$sample_name[i]
    sample_project <- samples_df$project_id[i]

    cat(paste0(rep("-", 78), collapse = ""), "\n")
    cat("Processing sample:", sample_name, "\n")

    # Initialize metrics dictionary for this sample
    sample_metrics <- list(
        sample_name = sample_name,
        project_id = sample_project
    )

    # Add metadata from samples.tsv
    meta_cols <- c("dataset_id", "dataset_type", "condition", "patient",
                   "batch", "seq_type", "raw_data")
    for (col in meta_cols) {
        if (col %in% colnames(samples_df)) {
            sample_metrics[[col]] <- samples_df[[col]][i]
        }
    }

    # Add info columns dynamically
    info_cols <- grep("^info_", colnames(samples_df), value = TRUE)
    for (col in info_cols) {
        sample_metrics[[col]] <- samples_df[[col]][i]
    }

    # === FASTP QC ===
    fastp_json <- file.path(sample_project, "qc", "fastp", sample_name,
                            paste0(sample_name, ".fastp.json"))
    if (file.exists(fastp_json)) {
        cat("  Processing fastp JSON...\n")
        tryCatch({
            fastp_data <- fromJSON(fastp_json)

            # Basic statistics
            sample_metrics$fastp_total_reads <- fastp_data$summary$before_filtering$total_reads
            sample_metrics$fastp_total_bases <- fastp_data$summary$before_filtering$total_bases
            sample_metrics$fastp_q20_rate <- fastp_data$summary$before_filtering$q20_rate
            sample_metrics$fastp_q30_rate <- fastp_data$summary$before_filtering$q30_rate
            sample_metrics$fastp_gc_content <- fastp_data$summary$before_filtering$gc_content

            # After filtering
            sample_metrics$fastp_total_reads_after <- fastp_data$summary$after_filtering$total_reads
            sample_metrics$fastp_total_bases_after <- fastp_data$summary$after_filtering$total_bases
            sample_metrics$fastp_q20_rate_after <- fastp_data$summary$after_filtering$q20_rate
            sample_metrics$fastp_q30_rate_after <- fastp_data$summary$after_filtering$q30_rate
            sample_metrics$fastp_gc_content_after <- fastp_data$summary$after_filtering$gc_content

            # Adapter content
            if (length(fastp_data$adapter_cutting) > 0) {
                sample_metrics$fastp_adapter_reads <- fastp_data$adapter_cutting$adapter_trimmed_reads
            }

            cat("    fastp QC loaded\n")
        }, error = function(e) {
            cat("    Error loading fastp:", e$message, "\n")
        })
    }

    # === STAR Alignment ===
    star_log <- file.path(sample_project, "qc", "STAR", sample_name,
                         paste0(sample_name, ".Log.final.out"))
    if (file.exists(star_log)) {
        cat("  Processing STAR alignment log...\n")
        tryCatch({
            star_lines <- readLines(star_log)

            # Parse key metrics
            for (line in star_lines) {
                if (grepl("Uniquely mapped reads number", line)) {
                    sample_metrics$star_uniquely_mapped <-
                        as.integer(gsub(".*\\|\\s*(\\d+).*", "\\1", line))
                }
                if (grepl("Uniquely mapped reads %", line)) {
                    sample_metrics$star_uniquely_mapped_pct <-
                        as.numeric(gsub(".*\\|\\s*([0-9.]+)%.*", "\\1", line))
                }
                if (grepl("Number of reads mapped to multiple loci", line)) {
                    sample_metrics$star_multi_mapped <-
                        as.integer(gsub(".*\\|\\s*(\\d+).*", "\\1", line))
                }
                if (grepl("% of reads mapped to too many loci", line)) {
                    sample_metrics$star_too_many_loci_pct <-
                        as.numeric(gsub(".*\\|\\s*([0-9.]+)%.*", "\\1", line))
                }
                if (grepl("Unmapped reads number", line)) {
                    sample_metrics$star_unmapped <-
                        as.integer(gsub(".*\\|\\s*(\\d+).*", "\\1", line))
                }
                if (grepl("% of reads unmapped: too many mismatches", line)) {
                    sample_metrics$star_unmapped_mismatch_pct <-
                        as.numeric(gsub(".*\\|\\s*([0-9.]+)%.*", "\\1", line))
                }
                if (grepl("% of reads unmapped: too short", line)) {
                    sample_metrics$star_unmapped_short_pct <-
                        as.numeric(gsub(".*\\|\\s*([0-9.]+)%.*", "\\1", line))
                }
                if (grepl("Number of splicing events:", line)) {
                    sample_metrics$star_splicing_events <-
                        as.integer(gsub(".*:\\s*(\\d+).*", "\\1", line))
                }
            }

            cat("    STAR alignment metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading STAR log:", e$message, "\n")
        })
    }

    # === Qualimap RNA-seq ===
    qualimap_file <- file.path(sample_project, "qc", "qualimap-rnaseq", sample_name,
                                "rnaseq_qc_results.txt")
    if (file.exists(qualimap_file)) {
        cat("  Processing Qualimap RNA-seq results...\n")
        tryCatch({
            qualimap_data <- read.delim(qualimap_file, comment.char = "#",
                                        stringsAsFactors = FALSE, header = FALSE)
            colnames(qualimap_data) <- c("metric", "value")

            # Extract key metrics
            for (j in 1:nrow(qualimap_data)) {
                metric <- qualimap_data$metric[j]
                value <- qualimap_data$value[j]

                if (grepl("Total number of reads", metric)) {
                    sample_metrics$qualimap_total_reads <- value
                }
                if (grepl("Mapped reads", metric)) {
                    sample_metrics$qualimap_mapped_reads <- value
                }
                if (grepl("Mapping rate", metric)) {
                    sample_metrics$qualimap_mapping_rate <- as.numeric(sub("%", "", value))
                }
                if (grepl("Mean coverage", metric)) {
                    sample_metrics$qualimap_mean_coverage <- as.numeric(value)
                }
                if (grepl("Mean insert size", metric)) {
                    sample_metrics$qualimap_mean_insert_size <- as.numeric(value)
                }
                if (grepl("Duplication rate", metric)) {
                    sample_metrics$qualimap_duplication_rate <- as.numeric(sub("%", "", value))
                }
            }

            cat("    Qualimap metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Qualimap:", e$message, "\n")
        })
    }

    # === Picard Alignment Summary ===
    picard_align <- file.path(sample_project, "qc", "picard", sample_name,
                              paste0(sample_name, ".alignment_summary_metrics.txt"))
    if (file.exists(picard_align)) {
        cat("  Processing Picard alignment metrics...\n")
        tryCatch({
            picard_align_data <- read.delim(picard_align, comment.char = "#",
                                            stringsAsFactors = FALSE)
            # Get first of data (after header)
            if (nrow(picard_align_data) > 0) {
                first_cat <- picard_align_data[1, ]
                sample_metrics$picard_total_reads <- first_cat$FIRST_OF_PAIR_READS if "FIRST_OF_PAIR_READS" %in% colnames(picard_align_data) else NA
                sample_metrics$picard_pct_aligned <- first_cat$PCT_PF_READS_ALIGNED if "PCT_PF_READS_ALIGNED" %in% colnames(picard_align_data) else NA
            }

            cat("    Picard alignment metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Picard alignment:", e$message, "\n")
        })
    }

    # === Picard RNA-seq Metrics ===
    picard_rna <- file.path(sample_project, "qc", "picard", sample_name,
                            paste0(sample_name, ".rnaseq_metrics.txt"))
    if (file.exists(picard_rna)) {
        cat("  Processing Picard RNA-seq metrics...\n")
        tryCatch({
            picard_rna_data <- read.delim(picard_rna, comment.char = "#",
                                          stringsAsFactors = FALSE)
            # Skip header rows
            if (nrow(picard_rna_data) > 1) {
                metrics_row <- picard_rna_data[2, ]
                sample_metrics$picard_pct_rRNA <- metrics_row$PCT_RIBOSOMAL_BASES if "PCT_RIBOSOMAL_BASES" %in% colnames(picard_rna_data) else NA
                sample_metrics$picard_pct_mRNA <- metrics_row$PCT_MRNA_BASES if "PCT_MRNA_BASES" %in% colnames(picard_rna_data) else NA
                sample_metrics$picard_pct_intronic <- metrics_row$PCT_INTRONIC_BASES if "PCT_INTRONIC_BASES" %in% colnames(picard_rna_data) else NA
                sample_metrics$picard_pct_intergenic <- metrics_row$PCT_INTERGENIC_BASES if "PCT_INTERGENIC_BASES" %in% colnames(picard_rna_data) else NA
                sample_metrics$picard_median_5prime <- metrics_row$MEDIAN_5PRIME_BIAS if "MEDIAN_5PRIME_BIAS" %in% colnames(picard_rna_data) else NA
                sample_metrics$picard_median_3prime <- metrics_row$MEDIAN_3PRIME_BIAS if "MEDIAN_3PRIME_BIAS" %in% colnames(picard_rna_data) else NA
            }

            cat("    Picard RNA-seq metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Picard RNA-seq:", e$message, "\n")
        })
    }

    # === Picard Insert Size ===
    picard_insert <- file.path(sample_project, "qc", "picard", sample_name,
                               paste0(sample_name, ".insert_size_metrics.txt"))
    if (file.exists(picard_insert)) {
        cat("  Processing Picard insert size metrics...\n")
        tryCatch({
            picard_insert_data <- read.delim(picard_insert, comment.char = "#",
                                             stringsAsFactors = FALSE)
            # Skip header rows and insert size histogram
            if (nrow(picard_insert_data) > 1) {
                metrics_row <- picard_insert_data[2, ]
                sample_metrics$picard_median_insert_size <- metrics_row$MEDIAN_INSERT_SIZE if "MEDIAN_INSERT_SIZE" %in% colnames(picard_insert_data) else NA
                sample_metrics$picard_mean_insert_size <- metrics_row$MEAN_INSERT_SIZE if "MEAN_INSERT_SIZE" %in% colnames(picard_insert_data) else NA
                sample_metrics$picard_min_insert_size <- metrics_row$MIN_INSERT_SIZE if "MIN_INSERT_SIZE" %in% colnames(picard_insert_data) else NA
                sample_metrics$picard_max_insert_size <- metrics_row$MAX_INSERT_SIZE if "MAX_INSERT_SIZE" %in% colnames(picard_insert_data) else NA
            }

            cat("    Picard insert size metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Picard insert size:", e$message, "\n")
        })
    }

    # === Picard GC Bias ===
    picard_gc <- file.path(sample_project, "qc", "picard", sample_name,
                          paste0(sample_name, ".gc_bias_summary_metrics.txt"))
    if (file.exists(picard_gc)) {
        cat("  Processing Picard GC bias metrics...\n")
        tryCatch({
            picard_gc_data <- read.delim(picard_gc, comment.char = "#",
                                        stringsAsFactors = FALSE)
            # Get summary metrics
            if (nrow(picard_gc_data) > 0) {
                sample_metrics$picard_gc_bias <- picard_gc_data$GC_BIAS_METRIC[1] if "GC_BIAS_METRIC" %in% colnames(picard_gc_data) else NA
                sample_metrics$picard_at_dropout <- picard_gc_data$AT_DROPOUT_METRIC[1] if "AT_DROPOUT_METRIC" %in% colnames(picard_gc_data) else NA
            }

            cat("    Picard GC bias metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Picard GC bias:", e$message, "\n")
        })
    }

    # === RNA-SeQC 2 ===
    rnaseqc2_dir <- file.path(sample_project, "qc", "rnaseqc2", sample_name)
    if (dir.exists(rnaseqc2_dir)) {
        cat("  Processing RNA-SeQC 2 results...\n")
        rnaseqc2_metrics <- file.path(rnaseqc2_dir, "metrics.tsv")
        if (file.exists(rnaseqc2_metrics)) {
            tryCatch({
                rnaseqc2_data <- read.delim(rnaseqc2_metrics, stringsAsFactors = FALSE)
                if (nrow(rnaseqc2_data) > 0) {
                    sample_metrics$rnaseqc2_genes_detected <- rnaseqc2_data$Genes.Detected[1] if "Genes.Detected" %in% colnames(rnaseqc2_data) else NA
                    sample_metrics$rnaseqc2_expression_profiling_efficiency <- rnaseqc2_data$Expression.Profiling.Efficiency[1] if "Expression.Profiling.Efficiency" %in% colnames(rnaseqc2_data) else NA
                    sample_metrics$rnaseqc2_intragenic_rate <- rnaseqc2_data$Intragenic.rate[1] if "Intragenic.rate" %in% colnames(rnaseqc2_data) else NA
                    sample_metrics$rnaseqc2_exonic_rate <- rnaseqc2_data$Exonic.Rate[1] if "Exonic.Rate" %in% colnames(rnaseqc2_data) else NA
                    sample_metrics$rnaseqc2_rRNA_rate <- rnaseqc2_data$rRNA.rate[1] if "rRNA.rate" %in% colnames(rnaseqc2_data) else NA
                    sample_metrics$rnaseqc2_5prime_3prime_bias <- rnaseqc2_data$`5'.3'bias`[1] if "`5'.3'bias`" %in% colnames(rnaseqc2_data) else NA
                }

                cat("    RNA-SeQC 2 metrics loaded\n")
            }, error = function(e) {
                cat("    Error loading RNA-SeQC 2:", e$message, "\n")
            })
        }
    }

    # Add to list
    qc_metrics_list[[i]] <- sample_metrics
    cat("  Sample", sample_name, "complete\n\n")
}

# Convert list to data frame
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

cat("\n")
cat(paste0(rep("=", 78), collapse = ""), "\n")
cat("QC Summary Complete!\n")
cat(paste0(rep("=", 78), collapse = ""), "\n")

# Close sinks
sink()
sink(type = "message")
