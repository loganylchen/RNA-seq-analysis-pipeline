#!/usr/bin/env Rscript

# Comprehensive QC Summary and Visualization
# Aggregates QC metrics from multiple tools and generates labeled visualizations

# Load required libraries
suppressPackageStartupMessages({
    library(jsonlite)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(gridExtra)
    library(readr)
    library(scales)
    library(ggsci)
})

# Logging
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat("==============================================================\n")
cat("Comprehensive QC Summary and Visualization\n")
cat("==============================================================\n\n")

# Get input parameters
samples_file <- snakemake@params$samples
project <- snakemake@params$project

# Get output paths
output_summary <- snakemake@output$summary
output_dir <- dirname(output_summary)
output_prefix <- tools::file_path_sans_ext(basename(output_summary))

# Create output directory for figures
figures_dir <- file.path(output_dir, paste0(output_prefix, "_figures"))
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

cat("Project:", project, "\n")
cat("Samples file:", samples_file, "\n")
cat("Output summary:", output_summary, "\n")
cat("Output figures directory:", figures_dir, "\n\n")

# Load sample metadata
cat("Loading sample metadata...\n")
samples_df <- read.delim(samples_file, comment.char = "#", stringsAsFactors = FALSE)
samples_df <- samples_df[samples_df$project_id == project, ]
cat("  Loaded", nrow(samples_df), "samples\n")
cat("  Columns:", paste(colnames(samples_df), collapse = ", "), "\n\n")

# Get list of QC files
qc_files_list <- snakemake@input$files
cat("Processing", length(qc_files_list), "QC files\n\n")

# Initialize list to store QC metrics for all samples
qc_metrics_list <- list()
sample_qc_files <- list()

# Parse sample names from QC file paths
for (qc_file in qc_files_list) {
    # Extract sample name from path
    # Path format: {project}/qc/{tool}/{sample}/{file}
    path_parts <- strsplit(qc_file, "/")[[1]]
    if (length(path_parts) >= 4) {
        sample_idx <- which(path_parts == "qc")
        if (length(sample_idx) > 0 && sample_idx[1] + 2 < length(path_parts)) {
            sample_name <- path_parts[sample_idx[1] + 2]
            tool_name <- path_parts[sample_idx[1] + 1]

            if (!sample_name %in% names(sample_qc_files)) {
                sample_qc_files[[sample_name]] <- list()
            }
            sample_qc_files[[sample_name]][[tool_name]] <- qc_file
        }
    }
}

cat("Found QC files for", length(sample_qc_files), "samples\n\n")

# Process each sample
for (i in 1:length(sample_qc_files)) {
    sample_name <- names(sample_qc_files)[i]
    sample_project <- samples_df$project_id[samples_df$sample_name == sample_name]

    if (length(sample_project) == 0) {
        sample_project <- project
    }

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

    # === FASTP QC ===
    if ("fastp" %in% names(sample_qc_files[[sample_name]])) {
        fastp_json <- sample_qc_files[[sample_name]]$fastp
        fastp_json <- paste0(dirname(fastp_json), "/", sample_name, ".fastp.json")

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
                cat("    Error loading fastp:", conditionMessage(e), "\n")
            })
        }
    }

    # === STAR Alignment ===
    if ("STAR" %in% names(sample_qc_files[[sample_name]])) {
        star_log <- sample_qc_files[[sample_name]]$STAR
        star_log <- paste0(dirname(star_log), "/", sample_name, ".Log.final.out")

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
                cat("    Error loading STAR log:", conditionMessage(e), "\n")
            })
        }
    }

    # === Qualimap RNA-seq ===
    if ("qualimap-rnaseq" %in% names(sample_qc_files[[sample_name]])) {
        qualimap_file <- sample_qc_files[[sample_name]][['qualimap-rnaseq']]
        qualimap_file <- file.path(qualimap_file, "rnaseq_qc_results.txt")

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
                cat("    Error loading Qualimap:", conditionMessage(e), "\n")
            })
        }
    }

    # === Picard Alignment Summary ===
    if ("picard" %in% names(sample_qc_files[[sample_name]])) {
        picard_align <- sample_qc_files[[sample_name]]$picard
        picard_align <- file.path(picard_align, paste0(sample_name, ".alignment_summary_metrics.txt"))

        if (file.exists(picard_align)) {
            cat("  Processing Picard alignment metrics...\n")
            tryCatch({
                picard_align_data <- read.delim(picard_align, comment.char = "#",
                                                stringsAsFactors = FALSE)
                # Get first of data (after header)
                if (nrow(picard_align_data) > 0) {
                    first_cat <- picard_align_data[1, ]
                    sample_metrics$picard_total_reads <- if ("FIRST_OF_PAIR_READS" %in% colnames(picard_align_data)) first_cat$FIRST_OF_PAIR_READS else NA
                    sample_metrics$picard_pct_aligned <- if ("PCT_PF_READS_ALIGNED" %in% colnames(picard_align_data)) first_cat$PCT_PF_READS_ALIGNED else NA
                }

                cat("    Picard alignment metrics loaded\n")
            }, error = function(e) {
                cat("    Error loading Picard alignment:", conditionMessage(e), "\n")
            })
        }
    }

    # === Picard RNA-seq Metrics ===
    picard_rna <- file.path(dirname(sample_qc_files[[sample_name]]$picard[1]),
                            paste0(sample_name, ".rnaseq_metrics.txt"))
    if (file.exists(picard_rna)) {
        cat("  Processing Picard RNA-seq metrics...\n")
        tryCatch({
            picard_rna_data <- read.delim(picard_rna, comment.char = "#",
                                          stringsAsFactors = FALSE)
            # Skip header rows
            if (nrow(picard_rna_data) > 1) {
                metrics_row <- picard_rna_data[2, ]
                sample_metrics$picard_pct_rRNA <- if ("PCT_RIBOSOMAL_BASES" %in% colnames(picard_rna_data)) metrics_row$PCT_RIBOSOMAL_BASES else NA
                sample_metrics$picard_pct_mRNA <- if ("PCT_MRNA_BASES" %in% colnames(picard_rna_data)) metrics_row$PCT_MRNA_BASES else NA
                sample_metrics$picard_pct_intronic <- if ("PCT_INTRONIC_BASES" %in% colnames(picard_rna_data)) metrics_row$PCT_INTRONIC_BASES else NA
                sample_metrics$picard_pct_intergenic <- if ("PCT_INTERGENIC_BASES" %in% colnames(picard_rna_data)) metrics_row$PCT_INTERGENIC_BASES else NA
                sample_metrics$picard_median_5prime <- if ("MEDIAN_5PRIME_BIAS" %in% colnames(picard_rna_data)) metrics_row$MEDIAN_5PRIME_BIAS else NA
                sample_metrics$picard_median_3prime <- if ("MEDIAN_3PRIME_BIAS" %in% colnames(picard_rna_data)) metrics_row$MEDIAN_3PRIME_BIAS else NA
            }

            cat("    Picard RNA-seq metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Picard RNA-seq:", conditionMessage(e), "\n")
        })
    }

    # === Picard Insert Size ===
    picard_insert <- file.path(dirname(sample_qc_files[[sample_name]]$picard[1]),
                               paste0(sample_name, ".insert_size_metrics.txt"))
    if (file.exists(picard_insert)) {
        cat("  Processing Picard insert size metrics...\n")
        tryCatch({
            picard_insert_data <- read.delim(picard_insert, comment.char = "#",
                                             stringsAsFactors = FALSE)
            # Skip header rows and insert size histogram
            if (nrow(picard_insert_data) > 1) {
                metrics_row <- picard_insert_data[2, ]
                sample_metrics$picard_median_insert_size <- if ("MEDIAN_INSERT_SIZE" %in% colnames(picard_insert_data)) metrics_row$MEDIAN_INSERT_SIZE else NA
                sample_metrics$picard_mean_insert_size <- if ("MEAN_INSERT_SIZE" %in% colnames(picard_insert_data)) metrics_row$MEAN_INSERT_SIZE else NA
                sample_metrics$picard_min_insert_size <- if ("MIN_INSERT_SIZE" %in% colnames(picard_insert_data)) metrics_row$MIN_INSERT_SIZE else NA
                sample_metrics$picard_max_insert_size <- if ("MAX_INSERT_SIZE" %in% colnames(picard_insert_data)) metrics_row$MAX_INSERT_SIZE else NA
            }

            cat("    Picard insert size metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Picard insert size:", conditionMessage(e), "\n")
        })
    }

    # === Picard GC Bias ===
    picard_gc <- file.path(dirname(sample_qc_files[[sample_name]]$picard[1]),
                          paste0(sample_name, ".gc_bias_summary_metrics.txt"))
    if (file.exists(picard_gc)) {
        cat("  Processing Picard GC bias metrics...\n")
        tryCatch({
            picard_gc_data <- read.delim(picard_gc, comment.char = "#",
                                        stringsAsFactors = FALSE)
            # Get summary metrics
            if (nrow(picard_gc_data) > 0) {
                sample_metrics$picard_gc_bias <- if ("GC_BIAS_METRIC" %in% colnames(picard_gc_data)) picard_gc_data$GC_BIAS_METRIC[1] else NA
                sample_metrics$picard_at_dropout <- if ("AT_DROPOUT_METRIC" %in% colnames(picard_gc_data)) picard_gc_data$AT_DROPOUT_METRIC[1] else NA
            }

            cat("    Picard GC bias metrics loaded\n")
        }, error = function(e) {
            cat("    Error loading Picard GC bias:", conditionMessage(e), "\n")
        })
    }

    # === RNA-SeQC 2 ===
    if ("rnaseqc2" %in% names(sample_qc_files[[sample_name]])) {
        rnaseqc2_dir <- sample_qc_files[[sample_name]]$rnaseqc2

        cat("  Processing RNA-SeQC 2 results...\n")
        rnaseqc2_metrics <- file.path(rnaseqc2_dir, "metrics.tsv")
        if (file.exists(rnaseqc2_metrics)) {
            tryCatch({
                rnaseqc2_data <- read.delim(rnaseqc2_metrics, stringsAsFactors = FALSE)
                if (nrow(rnaseqc2_data) > 0) {
                    sample_metrics$rnaseqc2_genes_detected <- if ("Genes.Detected" %in% colnames(rnaseqc2_data)) rnaseqc2_data$Genes.Detected[1] else NA
                    sample_metrics$rnaseqc2_expression_profiling_efficiency <- if ("Expression.Profiling.Efficiency" %in% colnames(rnaseqc2_data)) rnaseqc2_data$Expression.Profiling.Efficiency[1] else NA
                    sample_metrics$rnaseqc2_intragenic_rate <- if ("Intragenic.rate" %in% colnames(rnaseqc2_data)) rnaseqc2_data$Intragenic.rate[1] else NA
                    sample_metrics$rnaseqc2_exonic_rate <- if ("Exonic.Rate" %in% colnames(rnaseqc2_data)) rnaseqc2_data$Exonic.Rate[1] else NA
                    sample_metrics$rnaseqc2_rRNA_rate <- if ("rRNA.rate" %in% colnames(rnaseqc2_data)) rnaseqc2_data$rRNA.rate[1] else NA
                    sample_metrics$rnaseqc2_5prime_3prime_bias <- if ("`5'.3'bias`" %in% colnames(rnaseqc2_data)) rnaseqc2_data$`5'.3'bias`[1] else NA
                }

                cat("    RNA-SeQC 2 metrics loaded\n")
            }, error = function(e) {
                cat("    Error loading RNA-SeQC 2:", conditionMessage(e), "\n")
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
cat("Generating QC Visualizations\n")
cat(paste0(rep("=", 78), collapse = ""), "\n\n")

# Function to create labeled plots
create_qc_plot <- function(data, metric_col, plot_title, y_label,
                           color_by = "condition", shape_by = NULL,
                           plot_type = "boxplot") {

    # Check if metric exists and has data
    if (!metric_col %in% colnames(data)) {
        cat("  Skipping", metric_col, "- not found in data\n")
        return(NULL)
    }

    # Build columns to extract
    cols_to_get <- c("sample_name", color_by, metric_col)
    if (!is.null(shape_by) && shape_by %in% colnames(data)) {
        cols_to_get <- c(cols_to_get, shape_by)
    }

    metric_data <- data[, cols_to_get, drop = FALSE]

    # Remove NA values
    metric_data <- metric_data[!is.na(metric_data[[metric_col]]), ]

    if (nrow(metric_data) == 0) {
        cat("  Skipping", metric_col, "- no data available\n")
        return(NULL)
    }

    # Create plot
    p <- ggplot(metric_data, aes_string(x = color_by, y = metric_col, fill = color_by))

    if (plot_type == "boxplot") {
        p <- p + geom_boxplot(alpha = 0.7) +
            geom_point(position = position_jitter(width = 0.2), size = 2)
    } else if (plot_type == "bar") {
        p <- p + geom_bar(stat = "identity", alpha = 0.7)
    }

    # Add shape if specified
    if (!is.null(shape_by) && shape_by %in% colnames(metric_data)) {
        p <- p + aes_string(shape = shape_by)
    }

    # Add labels
    p <- p + labs(
        title = plot_title,
        x = "",
        y = y_label
    ) +
    theme_minimal(base_size = 12) +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    ) +
    scale_fill_pal_d() +
    scale_color_pal_d()

    return(p)
}

# Define color and shape columns
color_col <- "condition"
shape_col <- "dataset_type"

# Check what columns are available for coloring
if (!color_col %in% colnames(qc_summary_df)) {
    color_col <- colnames(qc_summary_df)[4]  # Use first available metadata column
}

cat("Coloring by:", color_col, "\n")
cat("Shape by:", shape_col, "\n\n")

# Define metrics to visualize
metrics_to_plot <- list(
    list(col = "fastp_q30_rate", title = "Q30 Rate", ylab = "Q30 Rate (%)", type = "boxplot"),
    list(col = "fastp_gc_content", title = "GC Content", ylab = "GC Content (%)", type = "boxplot"),
    list(col = "star_uniquely_mapped_pct", title = "Uniquely Mapped Reads", ylab = "Uniquely Mapped (%)", type = "boxplot"),
    list(col = "star_multi_mapped", title = "Multi-Mapped Reads", ylab = "Number of Reads", type = "boxplot"),
    list(col = "star_unmapped", title = "Unmapped Reads", ylab = "Number of Reads", type = "boxplot"),
    list(col = "qualimap_mapping_rate", title = "Mapping Rate", ylab = "Mapping Rate (%)", type = "boxplot"),
    list(col = "qualimap_mean_coverage", title = "Mean Coverage", ylab = "Mean Coverage", type = "boxplot"),
    list(col = "picard_pct_rRNA", title = "Ribosomal RNA Content", ylab = "rRNA (%)", type = "boxplot"),
    list(col = "picard_pct_mRNA", title = "mRNA Content", ylab = "mRNA (%)", type = "boxplot"),
    list(col = "picard_pct_intronic", title = "Intronic Content", ylab = "Intronic (%)", type = "boxplot"),
    list(col = "picard_pct_intergenic", title = "Intergenic Content", ylab = "Intergenic (%)", type = "boxplot"),
    list(col = "picard_median_insert_size", title = "Insert Size", ylab = "Median Insert Size (bp)", type = "boxplot"),
    list(col = "rnaseqc2_genes_detected", title = "Genes Detected", ylab = "Number of Genes", type = "boxplot"),
    list(col = "rnaseqc2_expression_profiling_efficiency", title = "Expression Profiling Efficiency", ylab = "EPE", type = "boxplot")
)

# Generate individual plots
cat("Generating individual QC metric plots...\n")
plots_list <- list()

for (i in 1:length(metrics_to_plot)) {
    metric_info <- metrics_to_plot[[i]]
    cat("  Creating plot for:", metric_info$title, "\n")

    p <- create_qc_plot(
        data = qc_summary_df,
        metric_col = metric_info$col,
        plot_title = metric_info$title,
        y_label = metric_info$ylab,
        color_by = color_col,
        shape_by = shape_col,
        plot_type = metric_info$type
    )

    if (!is.null(p)) {
        # Save individual plot
        safe_name <- gsub("[^A-Za-z0-9]", "_", metric_info$title)
        plot_file_png <- file.path(figures_dir, paste0(safe_name, ".png"))
        plot_file_pdf <- file.path(figures_dir, paste0(safe_name, ".pdf"))

        ggsave(plot_file_png, p, width = 8, height = 6, dpi = 300)
        ggsave(plot_file_pdf, p, width = 8, height = 6)

        cat("    Saved:", plot_file_png, "\n")

        plots_list[[metric_info$title]] <- p
    }
}

# Create summary figure with multiple panels
cat("\nCreating summary multi-panel figure...\n")
n_plots <- min(length(plots_list), 12)  # Limit to 12 plots for readability

if (n_plots > 0) {
    # Arrange plots in grid
    n_cols <- min(4, n_plots)
    summary_plot <- do.call(gridExtra::grid.arrange, c(plots_list[1:n_plots], ncol = n_cols))

    summary_file_png <- file.path(figures_dir, "QC_summary_overview.png")
    summary_file_pdf <- file.path(figures_dir, "QC_summary_overview.pdf")

    ggsave(summary_file_png, summary_plot, width = n_cols * 4, height = ceiling(n_plots/n_cols) * 4, dpi = 300)
    ggsave(summary_file_pdf, summary_plot, width = n_cols * 4, height = ceiling(n_plots/n_cols) * 4)

    cat("  Saved summary figure:", summary_file_png, "\n")
}

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
cat("  Individual plots:", figures_dir, "/\n")
cat("  Summary overview:", file.path(figures_dir, "QC_summary_overview.png"), "\n")

# Close sinks
sink()
sink(type = "message")
