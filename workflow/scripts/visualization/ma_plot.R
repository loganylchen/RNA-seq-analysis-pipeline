#!/usr/bin/env Rscript
# MA Plot Visualization for DEG Results
# Creates MA plots for each DEG tool (DESeq2, edgeR, limma-trend, limma-voom)
# M: average expression (A = (log2(x) + log2(y))/2)
# A: log fold change (M = log2(x) - log2(y))

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("==============================================================\n")
cat("MA Plot Visualization\n")
cat("==============================================================\n\n")

suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(dplyr)
    library(readr)
    library(RColorBrewer)
    library(scales)
})

# Get parameters from Snakemake
deg_file <- snakemake@input[["deg_file"]]
output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]
project <- snakemake@params[["project"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
dataset <- snakemake@params[["dataset"]]
tool <- snakemake@params[["tool"]]

cat("Parameters:\n")
cat("  Project:", project, "\n")
cat("  Dataset:", dataset, "\n")
cat("  Tool:", tool, "\n")
cat("  Log2FC threshold:", log2fc_threshold, "\n")
cat("  Padj threshold:", padj_threshold, "\n")
cat("  DEG file:", deg_file, "\n")
cat("  Output PDF:", output_pdf, "\n")
cat("  Output PNG:", output_png, "\n\n")

# ============================================================================
# READ DEG DATA
# ============================================================================

cat("Reading DEG results...\n")
deg_data <- read_tsv(deg_file, show_col_types = FALSE)
cat("  Dimensions:", nrow(deg_data), "genes x", ncol(deg_data), "columns\n")

# Check required columns
required_cols <- c("baseMean", "log2FoldChange", "padj")
missing_cols <- setdiff(required_cols, colnames(deg_data))
if (length(missing_cols) > 0) {
    # Try alternative column names
    if ("baseMean" %in% missing_cols && "mean_counts" %in% colnames(deg_data)) {
        deg_data <- rename(deg_data, baseMean = mean_counts)
        missing_cols <- setdiff(missing_cols, "baseMean")
    }
    if (length(missing_cols) > 0) {
        stop("ERROR: Missing required columns: ", paste(missing_cols, collapse=", "))
    }
}

# ============================================================================
# PREPARE DATA FOR MA PLOT
# ============================================================================

cat("\nPreparing MA plot data...\n")

# Calculate A (average expression) = log2(baseMean)
# M (log fold change) = log2FoldChange
deg_data <- deg_data %>%
    mutate(
        A = log2(baseMean),
        M = log2FoldChange,
        significant = (padj < padj_threshold) &
                      (abs(log2FoldChange) >= log2fc_threshold),
        direction = case_when(
            log2FoldChange >= log2fc_threshold & padj < padj_threshold ~ "Up",
            log2FoldChange <= -log2fc_threshold & padj < padj_threshold ~ "Down",
            TRUE ~ "NS"
        )
    ) %>%
    filter(!is.infinite(A) & !is.infinite(M))

# Remove NA values for plotting
deg_data <- deg_data %>% filter(!is.na(A) & !is.na(M))

cat("  Valid genes for plotting:", nrow(deg_data), "\n")

# Count genes by category
n_up <- sum(deg_data$direction == "Up")
n_down <- sum(deg_data$direction == "Down")
n_ns <- sum(deg_data$direction == "NS")

cat("  Up-regulated:", n_up, "\n")
cat("  Down-regulated:", n_down, "\n")
cat("  Not significant:", n_ns, "\n")

# ============================================================================
# CREATE MA PLOT
# ============================================================================

cat("\nCreating MA plot...\n")

# Define colors
colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

# Create MA plot
p <- ggplot(deg_data, aes(x = A, y = M)) +
    geom_point(aes(color = direction), alpha = 0.4, size = 0.8) +
    scale_color_manual(
        values = colors,
        name = "Regulation",
        labels = c("Up" = paste0("Up (", n_up, ")"),
                   "Down" = paste0("Down (", n_down, ")"),
                   "NS" = paste0("NS (", n_ns, ")"))
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
    geom_hline(yintercept = c(-log2fc_threshold, log2fc_threshold),
               linetype = "dashed", color = "red", size = 0.3, alpha = 0.5) +
    labs(
        title = paste("MA Plot:", tools::toTitleCase(gsub("_", " ", tool)), "-", dataset),
        subtitle = paste0(project, "\n",
                       "Red lines: |log2FC| >= ", log2fc_threshold, ", padj < ", padj_threshold),
        x = expression(log[2]~(Mean Expression)),
        y = expression(log[2]~(Fold Change))
    ) +
    theme_bw(base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey90", size = 0.3),
        panel.grid.minor = element_line(color = "grey95", size = 0.3),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11)
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 3))
    )

# ============================================================================
# SAVE PLOTS
# ============================================================================

cat("\nSaving plots...\n")

# Save as PDF
ggsave(output_pdf, plot = p, width = 10, height = 8, dpi = 300)
cat("  Saved PDF:", output_pdf, "\n")

# Save as PNG
ggsave(output_png, plot = p, width = 10, height = 8, dpi = 300)
cat("  Saved PNG:", output_png, "\n")

# ============================================================================
# SUMMARY STATISTICS
# ============================================================================

cat("\n==============================================================\n")
cat("Summary Statistics\n")
cat("==============================================================\n\n")

cat("Tool:", tool, "\n")
cat("Dataset:", dataset, "\n")
cat("Project:", project, "\n\n")

cat("DEG Summary:\n")
cat("  Total genes tested:", nrow(deg_data) + sum(is.na(deg_data$A)) + sum(is.na(deg_data$M)), "\n")
cat("  Up-regulated (padj <", padj_threshold, ", log2FC >=", log2fc_threshold, "):", n_up, "\n")
cat("  Down-regulated (padj <", padj_threshold, ", log2FC <=-", log2fc_threshold, "):", n_down, "\n")
cat("  Not significant:", n_ns, "\n\n")

# Expression statistics
cat("Expression Statistics:\n")
cat("  Mean expression range:", round(min(deg_data$A, na.rm=TRUE), 2), "to",
    round(max(deg_data$A, na.rm=TRUE), 2), "log2(mean)\n")
cat("  Median log2FC:", round(median(deg_data$M), 3), "\n")
cat("  Log2FC range:", round(min(deg_data$M), 2), "to",
    round(max(deg_data$M), 2), "\n\n")

cat("==============================================================\n")
cat("MA Plot Generation Complete!\n")
cat("==============================================================\n\n")

sink()
sink(type="message")
close(log)
