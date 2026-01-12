#!/usr/bin/env Rscript

# QC Summary Visualization
# Generates bar plots, box plots, and pairwise comparisons for key QC metrics

# Load required libraries
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(ggsci)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(gridExtra)
})

# Logging
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con)
sink(log_con, type = "message")

cat("==============================================================\n")
cat("QC Summary Visualization\n")
cat("==============================================================\n\n")

# Get input parameters
summary_file <- snakemake@input$summary
samples_file <- snakemake@params$samples
project <- snakemake@params$project
figures_dir <- snakemake@output$figures_dir

cat("Project:", project, "\n")
cat("Summary file:", summary_file, "\n")
cat("Samples file:", samples_file, "\n")
cat("Output directory:", figures_dir, "\n\n")

# Create output directory
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Load QC summary data
cat("Loading QC summary data...\n")
qc_data <- fread(summary_file, sep = "\t", header = TRUE, na.strings = c("NA", "", " "))
cat("  Loaded", nrow(qc_data), "samples with", ncol(qc_data), "columns\n\n")

# Load sample metadata
cat("Loading sample metadata...\n")
samples_df <- read.delim(samples_file, comment.char = "#", stringsAsFactors = FALSE)
samples_df <- samples_df[samples_df$project_id == project, ]
cat("  Loaded", nrow(samples_df), "samples\n\n")

# Merge metadata with QC data
qc_data <- merge(qc_data, samples_df[, c("sample_name", "patient", "dataset_id", "condition")],
                 by = "sample_name", all.x = TRUE)

# Define metrics to visualize
metrics_to_plot <- c(
    "qualimap_rnaseq.5_3_bias",
    "rna_seqc.Genes_Detected",
    "rna_seqc.rRNA_rate",
    "rna_seqc.Exonic_Rate",
    "rna_seqc.Intronic_Rate",
    "rna_seqc.Intergenic_Rate",
    "rna_seqc.Intragenic_Rate",
    "picard_insertsizemetrics.summed_median",
    "picard_rnaseqmetrics.PCT_MRNA_BASES",
    "star.mapped_percent",
    "fastp.after_filtering_gc_content"
)

# Create metric labels for plotting
metric_labels <- c(
    "qualimap_rnaseq.5_3_bias" = "5'-3' Bias",
    "rna_seqc.Genes_Detected" = "Genes Detected",
    "rna_seqc.rRNA_rate" = "rRNA Rate",
    "rna_seqc.Exonic_Rate" = "Exonic Rate",
    "rna_seqc.Intronic_Rate" = "Intronic Rate",
    "rna_seqc.Intergenic_Rate" = "Intergenic Rate",
    "rna_seqc.Intragenic_Rate" = "Intragenic Rate",
    "picard_insertsizemetrics.summed_median" = "Insert Size (Median)",
    "picard_rnaseqmetrics.PCT_MRNA_BASES" = "mRNA Bases (%)",
    "star.mapped_percent" = "Mapped (%)",
    "fastp.after_filtering_gc_content" = "GC Content"
)

cat("Metrics to visualize:\n")
for (m in metrics_to_plot) {
    cat("  -", m, "\n")
}
cat("\n")

# Function to create bar plot with error bars
create_bar_plot <- function(data, metric, metric_label) {
    # Filter to available data
    plot_data <- data[, .(sample_name, dataset_id, condition, patient, value = get(metric))]
    plot_data <- plot_data[!is.na(value)]

    if (nrow(plot_data) == 0) {
        return(NULL)
    }

    # Create grouping variable
    plot_data[, group := paste(dataset_id, condition, sep = "_")]

    # Calculate summary statistics
    summary_stats <- plot_data[, .(
        mean = mean(value, na.rm = TRUE),
        sd = sd(value, na.rm = TRUE),
        n = .N,
        se = sd(value, na.rm = TRUE) / sqrt(.N)
    ), by = .(dataset_id, condition)]

    summary_stats[, ymin := pmax(mean - se, 0)]
    summary_stats[, ymax := mean + se]

    # Create bar plot
    p <- ggplot(summary_stats, aes(x = condition, y = mean, fill = dataset_id)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
        geom_errorbar(aes(ymin = ymin, ymax = ymax),
                      position = position_dodge(width = 0.9),
                      width = 0.25) +
        scale_fill_pal_d(palette = "default") +
        labs(
            title = metric_label,
            x = "Condition",
            y = "Mean ± SE"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text = element_text(size = 11),
            axis.title = element_text(size = 12),
            legend.position = "right",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        )

    return(p)
}

# Function to create box plot with statistical comparisons
create_box_plot <- function(data, metric, metric_label) {
    # Filter to available data
    plot_data <- data[, .(sample_name, dataset_id, condition, patient, value = get(metric))]
    plot_data <- plot_data[!is.na(value)]

    if (nrow(plot_data) == 0) {
        return(NULL)
    }

    # Create grouping variable for pairwise comparisons
    plot_data[, group := paste(dataset_id, condition, sep = "_")]

    # Check if we have multiple groups for comparison
    unique_groups <- unique(plot_data$group)

    p <- ggplot(plot_data, aes(x = condition, y = value, fill = dataset_id)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_point(position = position_jitterdodge(dodge.width = 0.9),
                   size = 2, shape = 21, alpha = 0.8) +
        scale_fill_pal_d(palette = "default") +
        labs(
            title = metric_label,
            x = "Condition",
            y = "Value"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text = element_text(size = 11),
            axis.title = element_text(size = 12),
            legend.position = "right",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        )

    # Add statistical comparisons if we have multiple groups
    if (length(unique_groups) > 1) {
        # Add pairwise comparisons
        p <- p + stat_compare_means(
            method = "wilcox.test",
            label = "p.signif",
            symnum.args = list(
                cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                symbols = c("***", "**", "*", "ns")
            ),
            step.increase = 0.08,
            tip.length = 0.01,
            size = 3
        )
    }

    return(p)
}

# Function to create patient-wise pairwise comparison plot
create_patient_plot <- function(data, metric, metric_label) {
    # Filter to available data
    plot_data <- data[, .(sample_name, dataset_id, condition, patient, value = get(metric))]
    plot_data <- plot_data[!is.na(value)]

    # Filter to patients with multiple samples
    patient_counts <- plot_data[, .N, by = patient]
    multi_sample_patients <- patient_counts[N > 1, patient]
    plot_data <- plot_data[patient %in% multi_sample_patients]

    if (nrow(plot_data) == 0) {
        return(NULL)
    }

    # Sort by patient and create color mapping
    plot_data <- plot_data[order(patient, dataset_id, condition)]
    plot_data[, group := paste(dataset_id, condition, sep = "_")]

    # Create grouped box plot by patient
    p <- ggplot(plot_data, aes(x = factor(patient), y = value, fill = group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_point(aes(group = group), position = position_dodge(width = 0.8),
                   size = 2, shape = 21, alpha = 0.8) +
        geom_line(aes(group = patient), position = position_dodge(width = 0.8),
                  alpha = 0.5, linetype = "dashed") +
        scale_fill_pal_d(palette = "default") +
        labs(
            title = paste(metric_label, "(by Patient)"),
            x = "Patient",
            y = "Value"
        ) +
        theme_bw() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
            axis.text.y = element_text(size = 11),
            axis.title = element_text(size = 12),
            legend.position = "right",
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 10)
        )

    # Add statistical test comparing across conditions
    unique_conditions <- unique(plot_data$condition)
    if (length(unique_conditions) > 1) {
        p <- p + stat_compare_means(
            aes(label = after_stat(p.signif)),
            method = "wilcox.test",
            symnum.args = list(
                cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                symbols = c("***", "**", "*", "ns")
            ),
            label.y = max(plot_data$value, na.rm = TRUE) * 1.1,
            size = 3
        )
    }

    return(p)
}

# Generate plots for each metric
cat("Generating plots...\n\n")

all_plots <- list()
metric_count <- 0

for (metric in metrics_to_plot) {
    metric_label <- metric_labels[metric]

    cat(sprintf("Processing: %s (%s)\n", metric_label, metric))

    # Check if metric exists in data
    if (!metric %in% colnames(qc_data)) {
        cat(sprintf("  WARNING: Metric '%s' not found in data. Skipping.\n\n", metric))
        next
    }

    # Check if we have any data for this metric
    non_na_count <- sum(!is.na(qc_data[[metric]]))
    if (non_na_count == 0) {
        cat(sprintf("  WARNING: No data available for metric '%s'. Skipping.\n\n", metric_label))
        next
    }

    cat(sprintf("  Available data points: %d\n", non_na_count))

    # Create bar plot
    bar_p <- create_bar_plot(qc_data, metric, metric_label)
    if (!is.null(bar_p)) {
        bar_file_pdf <- file.path(figures_dir, paste0(gsub(" ", "_", metric_label), "_barplot.pdf"))
        bar_file_png <- file.path(figures_dir, paste0(gsub(" ", "_", metric_label), "_barplot.png"))
        ggsave(bar_file_pdf, bar_p, width = 8, height = 6)
        ggsave(bar_file_png, bar_p, width = 8, height = 6, dpi = 300)
        cat(sprintf("  Saved bar plot: %s\n", basename(bar_file_pdf)))
        all_plots[[paste0(metric, "_bar")]] <- bar_p
    }

    # Create box plot with statistics
    box_p <- create_box_plot(qc_data, metric, metric_label)
    if (!is.null(box_p)) {
        box_file_pdf <- file.path(figures_dir, paste0(gsub(" ", "_", metric_label), "_boxplot.pdf"))
        box_file_png <- file.path(figures_dir, paste0(gsub(" ", "_", metric_label), "_boxplot.png"))
        ggsave(box_file_pdf, box_p, width = 8, height = 6)
        ggsave(box_file_png, box_p, width = 8, height = 6, dpi = 300)
        cat(sprintf("  Saved box plot: %s\n", basename(box_file_pdf)))
        all_plots[[paste0(metric, "_box")]] <- box_p
    }

    # Create patient-wise plot
    patient_p <- create_patient_plot(qc_data, metric, metric_label)
    if (!is.null(patient_p)) {
        patient_file_pdf <- file.path(figures_dir, paste0(gsub(" ", "_", metric_label), "_by_patient.pdf"))
        patient_file_png <- file.path(figures_dir, paste0(gsub(" ", "_", metric_label), "_by_patient.png"))
        ggsave(patient_file_pdf, patient_p, width = 10, height = 6)
        ggsave(patient_file_png, patient_p, width = 10, height = 6, dpi = 300)
        cat(sprintf("  Saved patient plot: %s\n", basename(patient_file_pdf)))
        all_plots[[paste0(metric, "_patient")]] <- patient_p
    }

    metric_count <- metric_count + 1
    cat("\n")
}

# Create combined summary plots
cat("Creating combined summary plots...\n")

# Combine all bar plots
bar_plots <- Filter(function(x) grepl("_bar$", names(x)), all_plots)
if (length(bar_plots) > 0) {
    n <- length(bar_plots)
    n_cols <- min(3, n)
    n_rows <- ceiling(n / n_cols)

    combined_bar_pdf <- file.path(figures_dir, "combined_barplots.pdf")
    combined_bar_png <- file.path(figures_dir, "combined_barplots.png")

    combined_bar <- do.call(gridExtra::grid.arrange, c(bar_plots, ncol = n_cols))
    ggsave(combined_bar_pdf, combined_bar, width = 8 * n_cols, height = 6 * n_rows)
    ggsave(combined_bar_png, combined_bar, width = 8 * n_cols, height = 6 * n_rows, dpi = 300)

    cat("  Saved combined bar plots: combined_barplots.pdf\n")
}

# Combine all box plots
box_plots <- Filter(function(x) grepl("_box$", names(x)), all_plots)
if (length(box_plots) > 0) {
    n <- length(box_plots)
    n_cols <- min(3, n)
    n_rows <- ceiling(n / n_cols)

    combined_box_pdf <- file.path(figures_dir, "combined_boxplots.pdf")
    combined_box_png <- file.path(figures_dir, "combined_boxplots.png")

    combined_box <- do.call(gridExtra::grid.arrange, c(box_plots, ncol = n_cols))
    ggsave(combined_box_pdf, combined_box, width = 8 * n_cols, height = 6 * n_rows)
    ggsave(combined_box_png, combined_box, width = 8 * n_cols, height = 6 * n_rows, dpi = 300)

    cat("  Saved combined box plots: combined_boxplots.pdf\n")
}

# Summary statistics
cat("\n")
cat("==============================================================\n")
cat("QC Summary Visualization Complete!\n")
cat("==============================================================\n\n")

cat("Total metrics processed:", metric_count, "\n")
cat("Total plots generated:", length(all_plots), "\n\n")

cat("Output directory:", figures_dir, "\n")
cat("Individual plots saved for each metric:\n")
cat("  - Bar plots (mean ± SE by dataset_id and condition)\n")
cat("  - Box plots with statistical comparisons\n")
cat("  - Patient-wise comparison plots\n\n")

cat("Combined summary plots:\n")
cat("  - combined_barplots.pdf/png\n")
cat("  - combined_boxplots.pdf/png\n\n")

# Close sinks
sink()
sink(type = "message")
