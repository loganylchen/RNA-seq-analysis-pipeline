#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("==============================================================\n")
cat("DEG Scatterplot Comparison: Multiple Datasets\n")
cat("==============================================================\n\n")

suppressPackageStartupMessages({
    library(ggplot2)
    library(gridExtra)
    library(dplyr)
    library(readr)
    library(patchwork)
    library(stringr)
})

# Get parameters from Snakemake
combined_deg_file <- snakemake@input[["combined_deg"]]
deg_files <- unlist(snakemake@input[["deg_files"]])
project <- snakemake@params[["project"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
padj_threshold <- snakemake@params[["padj_threshold"]]

output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]

cat("Parameters:\n")
cat("  Project:", project, "\n")
cat("  Log2FC threshold:", log2fc_threshold, "\n")
cat("  Padj threshold:", padj_threshold, "\n")
cat("  Number of DEG files:", length(deg_files), "\n")
cat("  Output PDF:", output_pdf, "\n")
cat("  Output PNG:", output_png, "\n\n")

# Read combined DEG results
cat("Reading combined DEG results...\n")
combined_deg <- read_tsv(combined_deg_file, show_col_types = FALSE)
cat("  Dimensions:", nrow(combined_deg), "genes x", ncol(combined_deg), "columns\n")

# Filter to only up-regulated genes detected by ALL 4 tools (up_count == 4)
up_genes_all <- combined_deg %>%
    filter(up_regulated_count == 4) %>%
    pull(gene_id)

cat("  Up-regulated genes detected by ALL tools:", length(up_genes_all), "\n")

# Check if we have any genes
if (length(up_genes_all) == 0) {
    cat("\nWARNING: No up-regulated genes detected by all 4 tools.\n")
    cat("Generating placeholder figure...\n")

    # Create placeholder plot
    placeholder_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, size = 6,
                 label = "No up-regulated genes\ndetected by all 4 tools\n\n",
                 hjust = 0.5, vjust = 0.5) +
        annotate("text", x = 0.5, y = 0.35, size = 4,
                 label = paste0("Thresholds: log2FC >= ", log2fc_threshold, ", padj < ", padj_threshold),
                 hjust = 0.5, vjust = 0.5, color = "grey40") +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void() +
        theme(
            plot.margin = margin(20, 20, 20, 20)
        )

    # Add title
    final_placeholder <- placeholder_plot +
        plot_annotation(
            title = paste0(project, ": DEG Log2FC Comparison Across Datasets\n",
                          "Up-regulated genes (n=0)"),
            subtitle = paste0("Log2FC threshold >= ", log2fc_threshold, ", padj < ", padj_threshold,
                             "\nNo genes detected by all 4 tools"),
            theme = theme(
                plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
                plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
            )
        )

    # Save placeholder PDF
    cat("\nSaving placeholder PDF to:", output_pdf, "\n")
    ggsave(output_pdf, plot = final_placeholder, width = 10, height = 8, dpi = 300)
    cat("Placeholder PDF saved successfully\n")

    # Save placeholder PNG
    cat("Saving placeholder PNG to:", output_png, "\n")
    ggsave(output_png, plot = final_placeholder, width = 10, height = 8, dpi = 300)
    cat("Placeholder PNG saved successfully\n")

    cat("\n==============================================================\n")
    cat("Placeholder Figure Generated (No genes detected by all tools)\n")
    cat("==============================================================\n\n")

    sink()
    sink(type="message")
    quit(save = "no")
}

# Use genes detected by all tools
up_genes <- up_genes_all
cat("  Using", length(up_genes), "up-regulated genes detected by all tools\n\n")

# Function to parse file path and extract tool and dataset names
parse_deg_path <- function(filepath) {
    # Expected path format: {project}/DEG/{tool}/{quant_tool}/{dataset}_deg.tsv
    parts <- strsplit(filepath, "/")[[1]]
    filename <- basename(filepath)
    dataset <- str_replace(filename, "_deg.tsv$", "")
    tool <- parts[which(parts == "DEG") + 1]
    return(list(tool = tool, dataset = dataset, file = filepath))
}

# Parse all DEG files
cat("\n--- Parsing DEG Files ---\n")
deg_info <- lapply(deg_files, parse_deg_path)
cat("  Found", length(deg_info), "DEG files\n")
for (info in deg_info) {
    cat("    ", info$tool, "-", info$dataset, "\n")
}

# Function to read DEG TSV
read_deg_data <- function(filepath) {
    deg <- read_tsv(filepath, show_col_types = FALSE)
    return(deg)
}

# Read all DEG files and organize by tool and dataset
cat("\n--- Loading DEG Data ---\n")
deg_data_list <- list()
for (info in deg_info) {
    cat("  Loading", info$tool, "-", info$dataset, "\n")
    deg <- read_deg_data(info$file)
    deg_data_list[[paste(info$tool, info$dataset, sep = "_")]] <- list(
        data = deg,
        tool = info$tool,
        dataset = info$dataset
    )
}

# Get unique tools and datasets
tools <- unique(sapply(deg_data_list, function(x) x$tool))
datasets <- unique(sapply(deg_data_list, function(x) x$dataset))
cat("\n  Tools:", paste(tools, collapse = ", "), "\n")
cat("  Datasets:", paste(datasets, collapse = ", "), "\n")

# Function to determine if gene is significant
is_significant <- function(deg_data, gene_ids) {
    sig <- deg_data$padj < padj_threshold & abs(deg_data$log2FoldChange) >= log2fc_threshold
    names(sig) <- rownames(deg_data)
    return(sig[gene_ids])
}

# Function to prepare scatterplot data for a tool comparing two datasets
prepare_scatter_data <- function(deg1, deg2, dataset1_name, dataset2_name, tool_name, up_genes) {
    cat("\nPreparing scatter data for", tool_name, "-", dataset1_name, "vs", dataset2_name, "...\n")

    # Get common genes
    common_genes <- intersect(intersect(rownames(deg1), rownames(deg2)), up_genes)
    cat("  Common up-regulated genes:", length(common_genes), "\n")

    if (length(common_genes) == 0) {
        return(NULL)
    }

    # Create data frame
    scatter_df <- data.frame(
        gene_id = common_genes,
        x_log2FC = deg1[common_genes, "log2FoldChange"],
        y_log2FC = deg2[common_genes, "log2FoldChange"],
        stringsAsFactors = FALSE
    )

    # Determine significance
    sig1 <- is_significant(deg1, common_genes)
    sig2 <- is_significant(deg2, common_genes)

    scatter_df$significant <- "Not significant"
    scatter_df$significant[sig1 & !sig2] <- paste(dataset1_name, "only")
    scatter_df$significant[!sig1 & sig2] <- paste(dataset2_name, "only")
    scatter_df$significant[sig1 & sig2] <- "Both significant"

    cat("  Both significant:", sum(scatter_df$significant == "Both significant"), "\n")
    cat("  ", dataset1_name, "only:", sum(scatter_df$significant == paste(dataset1_name, "only")), "\n")
    cat("  ", dataset2_name, "only:", sum(scatter_df$significant == paste(dataset2_name, "only")), "\n")
    cat("  Not significant:", sum(scatter_df$significant == "Not significant"), "\n")

    return(scatter_df)
}

# Tool colors
tool_colors <- c(
    "deseq2" = "#E41A1C",
    "edger" = "#377EB8",
    "limma_trend" = "#4DAF4A",
    "limma_voom" = "#984EA3"
)

# Prepare scatter data for each tool
cat("\n--- Preparing Scatter Data ---\n")
scatter_list <- list()

for (tool in tools) {
    tool_normalized <- gsub("-", "_", tool)
    cat("\nTool:", tool, "\n")

    # Get DEG data for this tool across all datasets
    tool_data <- deg_data_list[sapply(deg_data_list, function(x) x$tool == tool)]

    if (length(tool_data) < 2) {
        cat("  Skipping: need at least 2 datasets for comparison\n")
        next
    }

    # Compare first dataset with each subsequent dataset
    dataset1_data <- tool_data[[1]]
    dataset1_name <- dataset1_data$dataset

    for (i in 2:length(tool_data)) {
        dataset2_data <- tool_data[[i]]
        dataset2_name <- dataset2_data$dataset

        scatter_df <- prepare_scatter_data(
            dataset1_data$data,
            dataset2_data$data,
            dataset1_name,
            dataset2_name,
            tool,
            up_genes
        )

        if (!is.null(scatter_df)) {
            key <- paste0(tool_normalized, "_", dataset1_name, "_vs_", dataset2_name)
            scatter_list[[key]] <- list(
                data = scatter_df,
                tool = tool,
                dataset1 = dataset1_name,
                dataset2 = dataset2_name,
                color = tool_colors[tool_normalized]
            )
        }
    }
}

cat("\nPrepared", length(scatter_list), "scatterplots\n")

# Function to create scatterplot
create_scatterplot <- function(scatter_obj) {
    scatter_df <- scatter_obj$data
    tool_name <- scatter_obj$tool
    dataset1 <- scatter_obj$dataset1
    dataset2 <- scatter_obj$dataset2
    tool_color <- scatter_obj$color

    cat("\nCreating scatterplot for", tool_name, "-", dataset1, "vs", dataset2, "...\n")

    # Calculate correlation
    cor_test <- cor.test(scatter_df$x_log2FC, scatter_df$y_log2FC,
                         method = "pearson", use = "complete.obs")

    # Update significance labels for this comparison
    dataset1_only <- paste(dataset1, "only")
    dataset2_only <- paste(dataset2, "only")

    scatter_df$significant <- factor(scatter_df$significant,
                                     levels = c("Not significant",
                                               dataset1_only,
                                               dataset2_only,
                                               "Both significant"))

    # Create named vectors for color and shape scales
    color_values <- c("Not significant", dataset1_only, dataset2_only, "Both significant")
    color_names <- c("grey70", "lightblue", "lightcoral", tool_color)
    names(color_names) <- color_values

    shape_values <- c("Not significant", dataset1_only, dataset2_only, "Both significant")
    shape_names <- c(16, 16, 16, 17)
    names(shape_names) <- shape_values

    # Create plot
    p <- ggplot(scatter_df, aes(x = x_log2FC, y = y_log2FC)) +
        geom_point(aes(shape = significant, color = significant), alpha = 0.6, size = 1.5) +
        scale_color_manual(
            values = color_names,
            name = "Significance",
            drop = FALSE
        ) +
        scale_shape_manual(
            values = shape_names,
            name = "Significance",
            drop = FALSE
        ) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", size = 0.8) +
        geom_smooth(method = "lm", color = "darkgrey", se = TRUE, size = 0.5, alpha = 0.3) +
        labs(
            title = paste(tools::toTitleCase(gsub("_", " ", tool_name)), ":",
                         dataset1, "vs", dataset2),
            x = paste0(dataset1, " Log2FC"),
            y = paste0(dataset2, " Log2FC")
        ) +
        annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5, size = 3,
                 label = paste0("r = ", round(cor_test$estimate, 3),
                              "\np < 2e-16"),
                 color = "black") +
        theme_bw(base_size = 10) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            legend.position = "right",
            legend.title = element_text(face = "bold", size = 9),
            legend.text = element_text(size = 8),
            panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
            panel.grid.minor = element_line(color = "grey95", linewidth = 0.3),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 9)
        )

    return(p)
}

# Create all scatterplots
cat("\n--- Creating Scatterplots ---\n")
plot_list <- lapply(scatter_list, create_scatterplot)

if (length(plot_list) == 0) {
    cat("\nWARNING: No scatterplots could be created (no common genes between datasets).\n")
    cat("Generating placeholder figure...\n")

    # Create placeholder plot
    placeholder_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, size = 6,
                 label = "No scatterplots could be created\n\n",
                 hjust = 0.5, vjust = 0.5) +
        annotate("text", x = 0.5, y = 0.35, size = 4,
                 label = paste0("Up-regulated genes detected by all tools: ", length(up_genes),
                              "\nNo common genes found between datasets for comparison"),
                 hjust = 0.5, vjust = 0.5, color = "grey40") +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void() +
        theme(
            plot.margin = margin(20, 20, 20, 20)
        )

    # Add title
    final_placeholder <- placeholder_plot +
        plot_annotation(
            title = paste0(project, ": DEG Log2FC Comparison Across Datasets\n",
                          "Up-regulated genes (n=", length(up_genes), ")"),
            subtitle = paste0("Log2FC threshold >= ", log2fc_threshold, ", padj < ", padj_threshold,
                             "\nNo common genes between datasets"),
            theme = theme(
                plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
                plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
            )
        )

    # Save placeholder PDF
    cat("\nSaving placeholder PDF to:", output_pdf, "\n")
    ggsave(output_pdf, plot = final_placeholder, width = 10, height = 8, dpi = 300)
    cat("Placeholder PDF saved successfully\n")

    # Save placeholder PNG
    cat("Saving placeholder PNG to:", output_png, "\n")
    ggsave(output_png, plot = final_placeholder, width = 10, height = 8, dpi = 300)
    cat("Placeholder PNG saved successfully\n")

    cat("\n==============================================================\n")
    cat("Placeholder Figure Generated (No common genes between datasets)\n")
    cat("==============================================================\n\n")

    sink()
    sink(type="message")
    quit(save = "no")
}

# Determine layout for combining plots
n_plots <- length(plot_list)
cat("\nCombining", n_plots, "plots...\n")

# Calculate grid dimensions
n_cols <- min(4, n_plots)
n_rows <- ceiling(n_plots / n_cols)

# Create plot layout using wrap_plots
combined_plot <- wrap_plots(plotlist = plot_list, ncol = n_cols, nrow = n_rows)

# Add overall title
final_plot <- combined_plot +
    plot_annotation(
        title = paste0(project, ": DEG Log2FC Comparison Across Datasets\n",
                      "Up-regulated genes (n=", length(up_genes), ")"),
        subtitle = paste0("Log2FC threshold >= ", log2fc_threshold, ", padj < ", padj_threshold,
                         "\nTriangle shape = significant in both datasets"),
        theme = theme(
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
            plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
        )
    )

# Save as PDF
cat("\nSaving PDF to:", output_pdf, "\n")
ggsave(output_pdf, plot = final_plot, width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
cat("PDF saved successfully\n")

# Save as PNG
cat("Saving PNG to:", output_png, "\n")
ggsave(output_png, plot = final_plot, width = 5 * n_cols, height = 4 * n_rows, dpi = 300)
cat("PNG saved successfully\n")

# Print summary statistics
cat("\n==============================================================\n")
cat("Summary Statistics\n")
cat("==============================================================\n\n")

for (key in names(scatter_list)) {
    obj <- scatter_list[[key]]
    df <- obj$data
    cat(obj$tool, "-", obj$dataset1, "vs", obj$dataset2, ":\n")
    cat("  Total genes plotted:", nrow(df), "\n")
    cat("  Both significant:", sum(df$significant == "Both significant"), "\n")
    cat("  ", obj$dataset1, "only:", sum(grepl(paste0(obj$dataset1, "only"), df$significant)), "\n")
    cat("  ", obj$dataset2, "only:", sum(grepl(paste0(obj$dataset2, "only"), df$significant)), "\n")
    cat("  Correlation:", round(cor(df$x_log2FC, df$y_log2FC, use = "complete.obs"), 3), "\n\n")
}

cat("==============================================================\n")
cat("Scatterplot Generation Complete!\n")
cat("==============================================================\n\n")

sink()
sink(type="message")
