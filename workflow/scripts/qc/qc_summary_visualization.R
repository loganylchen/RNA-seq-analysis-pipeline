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
figures_dir <- snakemake@output$figures_dir


cat("Summary file:", summary_file, "\n")

cat("Output directory:", figures_dir, "\n\n")

# Create output directory
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Load QC summary data
cat("Loading QC summary data...\n")
qc_data <- fread(summary_file, sep = "\t", header = TRUE, na.strings = c("NA", "", " "))
cat("  Loaded", nrow(qc_data), "samples with", ncol(qc_data), "columns\n")
cat("  Column names:\n")
print(colnames(qc_data))
cat("\n")

# Load sample metadata



# Check if key columns exist after merge
cat("Checking key columns in merged data...\n")
for (col in c("sample_name", "dataset_id", "condition", "patient")) {
    if (col %in% colnames(qc_data)) {
        cat("  ✓", col, "- found\n")
        cat("    Unique values:", length(unique(qc_data[[col]])), "\n")
        if (length(unique(qc_data[[col]])) <= 10) {
            cat("    Values:", paste(unique(qc_data[[col]]), collapse = ", "), "\n")
        }
    } else {
        cat("  ✗", col, "- NOT FOUND!\n")
    }
}
cat("\n")

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
    "qualimap_rnaseq.5_3_bias" = "5prime_3prime_Bias",
    "rna_seqc.Genes_Detected" = "Genes_Detected",
    "rna_seqc.rRNA_rate" = "rRNA_Rate",
    "rna_seqc.Exonic_Rate" = "Exonic_Rate",
    "rna_seqc.Intronic_Rate" = "Intronic_Rate",
    "rna_seqc.Intergenic_Rate" = "Intergenic_Rate",
    "rna_seqc.Intragenic_Rate" = "Intragenic_Rate",
    "picard_insertsizemetrics.summed_median" = "Insert_Size_Median",
    "picard_rnaseqmetrics.PCT_MRNA_BASES" = "mRNA_Bases_Pct",
    "star.mapped_percent" = "Mapped_Percent",
    "fastp.after_filtering_gc_content" = "GC_Content"
)

cat("Metrics to visualize:\n")
for (m in metrics_to_plot) {
    status <- if (m %in% colnames(qc_data)) "✓" else "✗"
    cat(" ", status, m)
    if (m %in% colnames(qc_data)) {
        non_na <- sum(!is.na(qc_data[[m]]))
        cat(sprintf(" (%d non-NA values)", non_na))
    }
    cat("\n")
}
cat("\n")

# Function to create bar plot by patient
create_bar_plot <- function(data, metric, metric_label) {
    cat("  -> create_bar_plot called for metric:", metric, "\n")
    cat("     Data dimensions:", nrow(data), "x", ncol(data), "\n")
    cat("     Data columns:", paste(colnames(data), collapse = ", "), "\n")

    # Convert to data.frame for easier manipulation
    plot_df <- as.data.frame(data)

    # Check if metric exists
    if (!metric %in% colnames(plot_df)) {
        cat("     ERROR: Metric", metric, "not found in data!\n")
        return(NULL)
    }

    # Check if patient column exists
    if (!"patient" %in% colnames(plot_df)) {
        cat("     ERROR: 'patient' column not found in data!\n")
        return(NULL)
    }

    # Select required columns
    required_cols <- c("sample_name", "dataset_id", "condition", "patient", metric)
    missing_cols <- setdiff(required_cols, colnames(plot_df))
    if (length(missing_cols) > 0) {
        cat("     ERROR: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
        return(NULL)
    }

    plot_df <- plot_df[, required_cols]
    colnames(plot_df)[colnames(plot_df) == metric] <- "value"

    cat("     After column selection:", nrow(plot_df), "rows\n")

    # Remove NA values
    plot_df <- plot_df[!is.na(plot_df$value), ]
    cat("     After removing NA:", nrow(plot_df), "rows\n")

    if (nrow(plot_df) == 0) {
        cat("     WARNING: No data points after NA removal\n")
        return(NULL)
    }

    # Create grouping variable for color (dataset_id + condition)
    plot_df$group <- paste(plot_df$dataset_id, plot_df$condition, sep = "_")

    # Convert to factors
    plot_df$patient <- factor(plot_df$patient)
    if (!is.factor(plot_df$dataset_id)) {
        plot_df$dataset_id <- as.factor(plot_df$dataset_id)
    }
    if (!is.factor(plot_df$condition)) {
        plot_df$condition <- as.factor(plot_df$condition)
    }

    # Sort by patient for better visualization
    plot_df <- plot_df[order(plot_df$patient, plot_df$dataset_id, plot_df$condition), ]

    cat("     Number of patients:", length(unique(plot_df$patient)), "\n")
    cat("     Unique groups:", paste(unique(plot_df$group), collapse = ", "), "\n")

    # Create bar plot grouped by patient
    cat("     Creating ggplot bar plot by patient...\n")
    p <- ggplot(plot_df, aes(x = patient, y = value, fill = group)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.8) +
        scale_fill_npg() +
        labs(
            title = metric_label,
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

    cat("     Bar plot created successfully\n")
    return(p)
}

# Function to create box plot with statistical comparisons
create_box_plot <- function(data, metric, metric_label) {
    cat("  -> create_box_plot called for metric:", metric, "\n")

    # Convert to data.frame for easier manipulation
    plot_df <- as.data.frame(data)

    # Check if metric exists
    if (!metric %in% colnames(plot_df)) {
        cat("     ERROR: Metric", metric, "not found in data!\n")
        return(NULL)
    }

    # Select required columns
    required_cols <- c("sample_name", "dataset_id", "condition", "patient", metric)
    plot_df <- plot_df[, required_cols]
    colnames(plot_df)[colnames(plot_df) == metric] <- "value"

    # Remove NA values
    plot_df <- plot_df[!is.na(plot_df$value), ]

    if (nrow(plot_df) == 0) {
        cat("     WARNING: No data points after NA removal\n")
        return(NULL)
    }

    # Convert to factors if needed
    if (!is.factor(plot_df$dataset_id)) {
        plot_df$dataset_id <- as.factor(plot_df$dataset_id)
    }
    if (!is.factor(plot_df$condition)) {
        plot_df$condition <- as.factor(plot_df$condition)
    }

    # Create grouping variable for pairwise comparisons
    plot_df$group <- paste(plot_df$dataset_id, plot_df$condition, sep = "_")

    # Check if we have multiple groups for comparison
    unique_groups <- unique(plot_df$group)
    cat("     Unique groups:", length(unique_groups), "\n")
    print(unique_groups)

    # Create box plot
    cat("     Creating ggplot box plot...\n")
    p <- ggplot(plot_df, aes(x = condition, y = value, fill = dataset_id)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_point(position = position_jitterdodge(dodge.width = 0.9),
                   size = 2, shape = 21, alpha = 0.8) +
        scale_fill_npg() +
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
        cat("     Adding statistical comparisons...\n")
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

    cat("     Box plot created successfully\n")
    return(p)
}

# Function to create patient-wise pairwise comparison plot
create_patient_plot <- function(data, metric, metric_label) {
    cat("  -> create_patient_plot called for metric:", metric, "\n")

    # Convert to data.frame for easier manipulation
    plot_df <- as.data.frame(data)

    # Check if metric exists
    if (!metric %in% colnames(plot_df)) {
        cat("     ERROR: Metric", metric, "not found in data!\n")
        return(NULL)
    }

    # Select required columns
    required_cols <- c("sample_name", "dataset_id", "condition", "patient", metric)
    plot_df <- plot_df[, required_cols]
    colnames(plot_df)[colnames(plot_df) == metric] <- "value"

    # Remove NA values
    plot_df <- plot_df[!is.na(plot_df$value), ]

    if (nrow(plot_df) == 0) {
        cat("     WARNING: No data points after NA removal\n")
        return(NULL)
    }

    # Filter to patients with multiple samples
    patient_counts <- plot_df %>%
        group_by(patient) %>%
        summarise(n = n(), .groups = "drop")

    multi_sample_patients <- patient_counts$patient[patient_counts$n > 1]
    cat("     Patients with multiple samples:", length(multi_sample_patients), "\n")

    plot_df <- plot_df[plot_df$patient %in% multi_sample_patients, ]

    if (nrow(plot_df) == 0) {
        cat("     WARNING: No multi-sample patients found\n")
        return(NULL)
    }

    # Sort by patient
    plot_df <- plot_df[order(plot_df$patient, plot_df$dataset_id, plot_df$condition), ]
    plot_df$group <- paste(plot_df$dataset_id, plot_df$condition, sep = "_")

    # Convert to factors
    plot_df$patient <- factor(plot_df$patient)
    if (!is.factor(plot_df$dataset_id)) {
        plot_df$dataset_id <- as.factor(plot_df$dataset_id)
    }
    if (!is.factor(plot_df$condition)) {
        plot_df$condition <- as.factor(plot_df$condition)
    }

    # Create grouped box plot by patient
    cat("     Creating patient-wise plot...\n")
    p <- ggplot(plot_df, aes(x = patient, y = value, fill = group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_point(aes(group = group), position = position_dodge(width = 0.8),
                   size = 2, shape = 21, alpha = 0.8) +
        geom_line(aes(group = patient), position = position_dodge(width = 0.8),
                  alpha = 0.5, linetype = "dashed") +
        scale_fill_npg() +
        labs(
            title = paste(metric_label, "- by Patient"),
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
    unique_conditions <- unique(plot_df$condition)
    if (length(unique_conditions) > 1) {
        cat("     Adding overall statistical test...\n")
        max_val <- max(plot_df$value, na.rm = TRUE)
        p <- p + stat_compare_means(
            aes(label = after_stat(p.signif)),
            method = "wilcox.test",
            symnum.args = list(
                cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                symbols = c("***", "**", "*", "ns")
            ),
            label.y = max_val * 1.1,
            size = 3
        )
    }

    cat("     Patient plot created successfully\n")
    return(p)
}

# Generate plots for each metric
cat("==============================================================\n")
cat("Generating plots...\n")
cat("==============================================================\n\n")

all_plots <- list()
metric_count <- 0
skipped_count <- 0

for (metric in metrics_to_plot) {
    metric_label <- metric_labels[metric]

    cat("\n")
    cat("--------------------------------------------------------------\n")
    cat(sprintf("Processing: %s (%s)\n", metric_label, metric))
    cat("--------------------------------------------------------------\n")

    # Check if metric exists in data
    if (!metric %in% colnames(qc_data)) {
        cat(sprintf("WARNING: Metric '%s' not found in data. Skipping.\n", metric))
        skipped_count <- skipped_count + 1
        next
    }

    # Check if we have any data for this metric
    non_na_count <- sum(!is.na(qc_data[[metric]]))
    if (non_na_count == 0) {
        cat(sprintf("WARNING: No data available for metric '%s'. Skipping.\n", metric_label))
        skipped_count <- skipped_count + 1
        next
    }

    cat(sprintf("Available data points: %d\n", non_na_count))

    # Try to create bar plot
    cat("\nCreating bar plot...\n")
    tryCatch({
        bar_p <- create_bar_plot(qc_data, metric, metric_label)
        if (!is.null(bar_p)) {
            safe_label <- gsub("'", "_", gsub('"', "_", metric_label))
            bar_file_pdf <- file.path(figures_dir, paste0(safe_label, "_barplot.pdf"))
            bar_file_png <- file.path(figures_dir, paste0(safe_label, "_barplot.png"))
            ggsave(bar_file_pdf, bar_p, width = 8, height = 6)
            ggsave(bar_file_png, bar_p, width = 8, height = 6, dpi = 300)
            cat(sprintf("  ✓ Saved bar plot: %s\n", basename(bar_file_pdf)))
            all_plots[[paste0(metric, "_bar")]] <- bar_p
        } else {
            cat("  ✗ Bar plot returned NULL\n")
        }
    }, error = function(e) {
        cat("  ✗ ERROR creating bar plot:", conditionMessage(e), "\n")
    })

    # Try to create box plot
    cat("\nCreating box plot...\n")
    tryCatch({
        box_p <- create_box_plot(qc_data, metric, metric_label)
        if (!is.null(box_p)) {
            safe_label <- gsub("'", "_", gsub('"', "_", metric_label))
            box_file_pdf <- file.path(figures_dir, paste0(safe_label, "_boxplot.pdf"))
            box_file_png <- file.path(figures_dir, paste0(safe_label, "_boxplot.png"))
            ggsave(box_file_pdf, box_p, width = 8, height = 6)
            ggsave(box_file_png, box_p, width = 8, height = 6, dpi = 300)
            cat(sprintf("  ✓ Saved box plot: %s\n", basename(box_file_pdf)))
            all_plots[[paste0(metric, "_box")]] <- box_p
        } else {
            cat("  ✗ Box plot returned NULL\n")
        }
    }, error = function(e) {
        cat("  ✗ ERROR creating box plot:", conditionMessage(e), "\n")
    })

    # Try to create patient-wise plot
    cat("\nCreating patient-wise plot...\n")
    tryCatch({
        patient_p <- create_patient_plot(qc_data, metric, metric_label)
        if (!is.null(patient_p)) {
            safe_label <- gsub("'", "_", gsub('"', "_", metric_label))
            patient_file_pdf <- file.path(figures_dir, paste0(safe_label, "_by_patient.pdf"))
            patient_file_png <- file.path(figures_dir, paste0(safe_label, "_by_patient.png"))
            ggsave(patient_file_pdf, patient_p, width = 10, height = 6)
            ggsave(patient_file_png, patient_p, width = 10, height = 6, dpi = 300)
            cat(sprintf("  ✓ Saved patient plot: %s\n", basename(patient_file_pdf)))
            all_plots[[paste0(metric, "_patient")]] <- patient_p
        } else {
            cat("  ✗ Patient plot returned NULL\n")
        }
    }, error = function(e) {
        cat("  ✗ ERROR creating patient plot:", conditionMessage(e), "\n")
    })

    metric_count <- metric_count + 1
    cat(sprintf("Completed metric %d: %s\n\n", metric_count, metric_label))
}

# Create combined summary plots
cat("\n")
cat("==============================================================\n")
cat("Creating combined summary plots...\n")
cat("==============================================================\n\n")

# Combine all bar plots
bar_plots <- Filter(function(x) grepl("_bar$", names(x)), all_plots)
cat("Found", length(bar_plots), "bar plots\n")
if (length(bar_plots) > 0) {
    n <- length(bar_plots)
    n_cols <- min(3, n)
    n_rows <- ceiling(n / n_cols)

    combined_bar_pdf <- file.path(figures_dir, "combined_barplots.pdf")
    combined_bar_png <- file.path(figures_dir, "combined_barplots.png")

    cat("  Creating combined bar plot grid (", n_cols, "x", n_rows, ")...\n", sep = "")
    tryCatch({
        combined_bar <- do.call(gridExtra::grid.arrange, c(bar_plots, ncol = n_cols))
        ggsave(combined_bar_pdf, combined_bar, width = 8 * n_cols, height = 6 * n_rows)
        ggsave(combined_bar_png, combined_bar, width = 8 * n_cols, height = 6 * n_rows, dpi = 300)
        cat("  ✓ Saved combined bar plots\n")
    }, error = function(e) {
        cat("  ✗ ERROR creating combined bar plots:", conditionMessage(e), "\n")
    })
} else {
    cat("  No bar plots to combine\n")
}

# Combine all box plots
box_plots <- Filter(function(x) grepl("_box$", names(x)), all_plots)
cat("Found", length(box_plots), "box plots\n")
if (length(box_plots) > 0) {
    n <- length(box_plots)
    n_cols <- min(3, n)
    n_rows <- ceiling(n / n_cols)

    combined_box_pdf <- file.path(figures_dir, "combined_boxplots.pdf")
    combined_box_png <- file.path(figures_dir, "combined_boxplots.png")

    cat("  Creating combined box plot grid (", n_cols, "x", n_rows, ")...\n", sep = "")
    tryCatch({
        combined_box <- do.call(gridExtra::grid.arrange, c(box_plots, ncol = n_cols))
        ggsave(combined_box_pdf, combined_box, width = 8 * n_cols, height = 6 * n_rows)
        ggsave(combined_box_png, combined_box, width = 8 * n_cols, height = 6 * n_rows, dpi = 300)
        cat("  ✓ Saved combined box plots\n")
    }, error = function(e) {
        cat("  ✗ ERROR creating combined box plots:", conditionMessage(e), "\n")
    })
} else {
    cat("  No box plots to combine\n")
}

# Summary statistics
cat("\n")
cat("==============================================================\n")
cat("QC Summary Visualization Complete!\n")
cat("==============================================================\n\n")

cat("Total metrics requested:", length(metrics_to_plot), "\n")
cat("Metrics processed:", metric_count, "\n")
cat("Metrics skipped:", skipped_count, "\n")
cat("Total individual plots generated:", length(all_plots), "\n\n")

cat("Output directory:", figures_dir, "\n\n")

cat("Individual plots saved for each metric:\n")
cat("  - Bar plots (by patient, colored by dataset_id and condition)\n")
cat("  - Box plots with statistical comparisons\n")
cat("  - Patient-wise comparison plots\n\n")

cat("Combined summary plots:\n")
cat("  - combined_barplots.pdf/png\n")
cat("  - combined_boxplots.pdf/png\n\n")

cat("List of generated plots:\n")
if (length(all_plots) > 0) {
    for (plot_name in names(all_plots)) {
        cat("  -", plot_name, "\n")
    }
} else {
    cat("  (no plots generated)\n")
}
cat("\n")

# Close sinks
sink()
sink(type = "message")
