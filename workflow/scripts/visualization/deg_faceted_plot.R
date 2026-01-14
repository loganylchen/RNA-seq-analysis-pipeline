#!/usr/bin/env Rscript
# Faceted DEG visualization with top genes labeled
# Creates a multi-panel figure comparing DEG results from all tools
# Labels the top 3 genes by log2FC in each tool

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(readr)
    library(stringr)
    library(ggrepel)
})

cat("==============================================================\n")
cat("Faceted DEG Visualization with Top Gene Labels\n")
cat("==============================================================\n\n")

# Get parameters from Snakemake
combined_deg_file <- snakemake@input[["combined_deg"]]
deg_files <- snakemake@input[["deg_files"]]
output_pdf <- snakemake@output[["pdf"]]
output_png <- snakemake@output[["png"]]
project <- snakemake@params[["project"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
padj_threshold <- snakemake@params[["padj_threshold"]]

cat("Parameters:\n")
cat("  Project:", project, "\n")
cat("  Combined DEG file:", combined_deg_file, "\n")
cat("  Number of DEG files:", length(deg_files), "\n")
cat("  log2FC threshold:", log2fc_threshold, "\n")
cat("  padj threshold:", padj_threshold, "\n")
cat("  Output PDF:", output_pdf, "\n")
cat("  Output PNG:", output_png, "\n\n")

# ============================================================================
# LOAD COMBINED DEG DATA
# ============================================================================

cat("Loading combined DEG data...\n")
combined_deg <- read.delim(combined_deg_file, stringsAsFactors = FALSE)
cat("  Loaded", nrow(combined_deg), "genes\n\n")

# ============================================================================
# LOAD INDIVIDUAL DEG FILES AND PREPARE DATA
# ============================================================================

cat("Loading individual DEG files...\n")

# Function to parse file path and extract tool name
parse_tool_from_path <- function(filepath) {
    parts <- strsplit(filepath, "/")[[1]]
    # Tool is between "DEG" and the dataset name
    deg_idx <- which(parts == "DEG")
    if (length(deg_idx) > 0) {
        return(parts[deg_idx + 1])
    }
    return(NA)
}

# Read and combine all DEG files
all_deg_data <- list()

for (i in seq_along(deg_files)) {
    deg_file <- deg_files[i]
    tool <- parse_tool_from_path(deg_file)

    cat("  Loading:", tool, "-", basename(deg_file), "\n")

    deg_data <- read.delim(deg_file, stringsAsFactors = FALSE)

    # Standardize column names across tools
    if ("log2FoldChange" %in% colnames(deg_data)) {
        deg_data$log2FC <- deg_data$log2FoldChange
    }
    if ("logFC" %in% colnames(deg_data)) {
        deg_data$log2FC <- deg_data$logFC
    }
    if ("FDR" %in% colnames(deg_data)) {
        deg_data$padj <- deg_data$FDR
    }
    if ("PValue" %in% colnames(deg_data)) {
        deg_data$pvalue <- deg_data$PValue
    }

    # Add tool identifier
    deg_data$tool <- tool

    # Determine significance
    deg_data <- deg_data %>%
        mutate(
            significant = (padj < padj_threshold) & (abs(log2FC) >= log2fc_threshold),
            direction = case_when(
                log2FC >= log2fc_threshold & padj < padj_threshold ~ "Up",
                log2FC <= -log2fc_threshold & padj < padj_threshold ~ "Down",
                TRUE ~ "NS"
            )
        )

    # Gene ID column might vary
    if ("gene_id" %in% colnames(deg_data)) {
        deg_data$gene <- deg_data$gene_id
    } else if ("rowname" %in% colnames(deg_data)) {
        deg_data$gene <- deg_data$rowname
    } else if ("X" %in% colnames(deg_data)) {
        deg_data$gene <- deg_data$X
    } else {
        deg_data$gene <- rownames(deg_data)
    }

    all_deg_data[[tool]] <- deg_data
}

cat("\n")

# ============================================================================
# IDENTIFY TOP 3 GENES BY LOG2FC FOR EACH TOOL
# ============================================================================

cat("Identifying top 3 genes by log2FC for each tool...\n")

top_genes_list <- list()

for (tool in names(all_deg_data)) {
    deg_data <- all_deg_data[[tool]]

    # Get top 3 up-regulated genes
    top_up <- deg_data %>%
        filter(significant & direction == "Up") %>%
        arrange(desc(log2FC)) %>%
        head(3) %>%
        pull(gene)

    # Get top 3 down-regulated genes
    top_down <- deg_data %>%
        filter(significant & direction == "Down") %>%
        arrange(log2FC) %>%
        head(3) %>%
        pull(gene)

    # Get top 3 by absolute log2FC (for general labeling)
    top_abs <- deg_data %>%
        arrange(desc(abs(log2FC))) %>%
        head(3) %>%
        pull(gene)

    top_genes_list[[tool]] <- unique(c(top_up, top_down, top_abs))

    cat("  ", tools::toTitleCase(tool), ":", paste(head(top_genes_list[[tool]], 3), collapse = ", "), "\n")
}

cat("\n")

# ============================================================================
# COMBINE DATA FOR PLOTTING
# ============================================================================

cat("Preparing data for plotting...\n")

# Combine all DEG data
combined_plot_data <- bind_rows(all_deg_data)

# Calculate baseMean if not present (average of all samples if available)
if (!"baseMean" %in% colnames(combined_plot_data)) {
    # Use absolute log2FC as size aesthetic
    combined_plot_data$baseMean <- abs(combined_plot_data$log2FC)
}

# Add label column for top genes
combined_plot_data <- combined_plot_data %>%
    mutate(
        is_top_gene = gene %in% unlist(top_genes_list),
        tool_label = tools::toTitleCase(gsub("_", " ", tool))
    )

# Create label text (only for top genes)
combined_plot_data <- combined_plot_data %>%
    group_by(tool) %>%
    mutate(
        gene_label = ifelse(
            is_top_gene,
            gene,
            NA
        )
    ) %>%
    ungroup()

cat("  Total genes in plot:", nrow(combined_plot_data), "\n")
cat("  Genes to label:", sum(!is.na(combined_plot_data$gene_label)), "\n\n")

# ============================================================================
# CREATE FACETED PLOT
# ============================================================================

cat("Creating faceted plot...\n")

# Define colors for tools
tool_colors <- c(
    "deseq2" = "#E41A1C",
    "edger" = "#377EB8",
    "limma_trend" = "#4DAF4A",
    "limma_voom" = "#984EA3"
)

# Define colors for direction
direction_colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

# Create the plot
p <- ggplot(combined_plot_data, aes(x = baseMean, y = log2FC)) +
    # Color by significance/direction
    geom_point(
        aes(color = direction, alpha = !is_top_gene),
        size = 0.8
    ) +
    # Highlight top genes with larger points
    geom_point(
        data = combined_plot_data %>% filter(is_top_gene),
        aes(color = direction),
        size = 2.5,
        shape = 1,
        stroke = 1
    ) +
    # Add text labels for top genes
    geom_text_repel(
        data = combined_plot_data %>% filter(!is.na(gene_label)),
        aes(label = gene_label),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.3,
        force = 2
    ) +
    # Facet by tool
    facet_wrap(~ tool_label, ncol = 2, scales = "free") +
    # Color scale
    scale_color_manual(
        values = direction_colors,
        name = "Regulation",
        drop = FALSE
    ) +
    scale_alpha_manual(
        values = c("TRUE" = 0.4, "FALSE" = 1),
        guide = "none"
    ) +
    # Threshold lines
    geom_hline(
        yintercept = c(-log2fc_threshold, log2fc_threshold),
        linetype = "dashed",
        color = "grey50",
        size = 0.5
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30", size = 0.3) +
    # Theme
    theme_bw(base_size = 12) +
    theme(
        legend.position = "right",
        strip.background = element_rect(fill = "grey90", color = "grey70"),
        strip.text = element_text(face = "bold", size = 11),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
        title = paste("Differential Expression Results -", project),
        x = expression("Mean Expression (log"[2]*" baseMean)"),
        y = expression("Log"[2]*" Fold Change"),
        color = "Regulation"
    )

# ============================================================================
# SAVE PLOT
# ============================================================================

cat("Saving plot...\n")

# Save as PDF
ggsave(
    filename = output_pdf,
    plot = p,
    width = 14,
    height = 12,
    units = "in",
    dpi = 300
)
cat("  Saved PDF:", output_pdf, "\n")

# Save as PNG
ggsave(
    filename = output_png,
    plot = p,
    width = 14,
    height = 12,
    units = "in",
    dpi = 300
)
cat("  Saved PNG:", output_png, "\n")

# ============================================================================
# GENERATE SUMMARY TABLE
# ============================================================================

cat("\nGenerating summary table...\n")

summary_table <- data.frame(
    Tool = character(),
    Total_DEGs = integer(),
    Up_regulated = integer(),
    Down_regulated = integer(),
    Top_up_gene = character(),
    Top_down_gene = character(),
    stringsAsFactors = FALSE
)

for (tool in names(all_deg_data)) {
    deg_data <- all_deg_data[[tool]]

    summary_table <- rbind(summary_table, data.frame(
        Tool = tools::toTitleCase(gsub("_", " ", tool)),
        Total_DEGs = sum(deg_data$significant),
        Up_regulated = sum(deg_data$direction == "Up"),
        Down_regulated = sum(deg_data$direction == "Down"),
        Top_up_gene = deg_data %>%
            filter(direction == "Up") %>%
            arrange(desc(log2FC)) %>%
            head(1) %>%
            pull(gene),
        Top_down_gene = deg_data %>%
            filter(direction == "Down") %>%
            arrange(log2FC) %>%
            head(1) %>%
            pull(gene),
        stringsAsFactors = FALSE
    ))
}

# Save summary table
summary_file <- sub("\\.(pdf|png)$", "_summary.tsv", output_pdf)
write.table(summary_table, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat("  Saved summary:", summary_file, "\n")

# Print summary
cat("\nSummary:\n")
print(summary_table)
cat("\n")

cat("==============================================================\n")
cat("Faceted DEG Visualization Complete!\n")
cat("==============================================================\n\n")

sink()
sink(type="message")
close(log)
