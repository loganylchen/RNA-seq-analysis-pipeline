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
    library(tibble)
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
# LOAD GENE ID TO GENE NAME MAPPING
# ============================================================================

cat("Loading gene ID to gene name mapping...\n")
gene_name_file <- "resources/gene_id_to_gene_name.tsv"

if (file.exists(gene_name_file)) {
    gene_name_map <- read.delim(gene_name_file, stringsAsFactors = FALSE)
    cat("  Loaded", nrow(gene_name_map), "gene mappings\n")
    # Create named vector for lookup
    gene_name_lookup <- setNames(gene_name_map$gene_name, gene_name_map$gene_id)
} else {
    cat("  Warning: Gene name file not found, using gene IDs\n")
    gene_name_lookup <- NULL
}
cat("\n")

# ============================================================================
# LOAD INDIVIDUAL DEG FILES AND PREPARE DATA
# ============================================================================

cat("Loading individual DEG files...\n")

# Get the target dataset from output path
target_dataset <- basename(sub("_faceted_plot\\.(pdf|png)$", "", output_pdf))
cat("  Target dataset:", target_dataset, "\n")

# Function to parse file path and extract tool name
parse_tool_from_path <- function(filepath) {
    parts <- strsplit(filepath, "/")[[1]]
    # Tool is between "DEG" and the quantification tool/dataset
    deg_idx <- which(parts == "DEG")
    if (length(deg_idx) > 0) {
        return(parts[deg_idx + 1])
    }
    return(NA)
}

# Function to check if file matches target dataset
matches_dataset <- function(filepath, target) {
    filename <- basename(filepath)
    # Check if filename starts with the quant tool and ends with target dataset
    # Format: {quant_tool}_{target}_deg.tsv
    return(grepl(paste0(target, "_deg\\.tsv$"), filename) ||
           grepl(paste0("_", target, "_deg\\.tsv$"), filename))
}

# Read and combine all DEG files
all_deg_data <- list()

for (i in seq_along(deg_files)) {
    deg_file <- deg_files[i]

    # Skip files that don't match target dataset
    if (!matches_dataset(deg_file, target_dataset)) {
        cat("  Skipping:", basename(deg_file), "(not target dataset)\n")
        next
    }

    tool <- parse_tool_from_path(deg_file)

    cat("  Loading:", tool, "-", basename(deg_file), "\n")

    deg_data <- read.delim(deg_file, stringsAsFactors = FALSE, row.names = 1)

    # Standardize column names across tools
    if ("log2FoldChange" %in% colnames(deg_data)) {
        deg_data$log2FC <- deg_data$log2FoldChange
    } else if ("logFC" %in% colnames(deg_data)) {
        deg_data$log2FC <- deg_data$logFC
    }

    if ("FDR" %in% colnames(deg_data)) {
        deg_data$padj <- deg_data$FDR
    } else if ("PValue" %in% colnames(deg_data)) {
        deg_data$padj <- deg_data$PValue
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

    # Add gene_id column from rownames
    deg_data$gene_id <- rownames(deg_data)

    all_deg_data[[tool]] <- deg_data

    cat("    Genes loaded:", nrow(deg_data), "\n")
    cat("    Significant:", sum(deg_data$significant), "\n")
}

cat("\n  Tools loaded:", paste(names(all_deg_data), collapse = ", "), "\n\n")

# ============================================================================
# IDENTIFY TOP 3 GENES BY ABSOLUTE LOG2FC FOR EACH TOOL
# ============================================================================

cat("Identifying top 3 genes by absolute log2FC for each tool...\n")

top_genes_list <- list()

for (tool in names(all_deg_data)) {
    deg_data <- all_deg_data[[tool]]

    # Get top 3 up-regulated genes
    top_up <- deg_data %>%
        filter(significant & direction == "Up") %>%
        arrange(desc(log2FC)) %>%
        head(3) %>%
        pull(gene_id)

    # Get top 3 down-regulated genes
    top_down <- deg_data %>%
        filter(significant & direction == "Down") %>%
        arrange(log2FC) %>%
        head(3) %>%
        pull(gene_id)

    # Get top 3 by absolute log2FC (for general labeling)
    top_abs <- deg_data %>%
        arrange(desc(abs(log2FC))) %>%
        head(3) %>%
        pull(gene_id)

    top_genes_list[[tool]] <- unique(c(top_up, top_down, top_abs))

    cat("  ", tools::toTitleCase(gsub("_", " ", tool)), ":", length(top_genes_list[[tool]]), "genes\n")
}

cat("\n")

# ============================================================================
# COMBINE DATA FOR PLOTTING
# ============================================================================

cat("Preparing data for plotting...\n")

# Combine all DEG data
combined_plot_data <- bind_rows(all_deg_data, .id = "source")

# Calculate expression if baseMean not available
if (!"baseMean" %in% colnames(combined_plot_data)) {
    # Use a default value for visualization
    combined_plot_data$baseMean <- 10
    cat("  Warning: baseMean not found, using default value\n")
}

# Log10 transform expression (add small value to avoid log(0))
combined_plot_data <- combined_plot_data %>%
    mutate(
        log10_expr = log10(baseMean + 1)
    )

# Add gene names
if (!is.null(gene_name_lookup)) {
    combined_plot_data <- combined_plot_data %>%
        mutate(
            gene_display = ifelse(
                gene_id %in% names(gene_name_lookup),
                gene_name_lookup[gene_id],
                gene_id
            )
        )
} else {
    combined_plot_data$gene_display <- combined_plot_data$gene_id
}

# Add label column for top genes
combined_plot_data <- combined_plot_data %>%
    mutate(
        is_top_gene = gene_id %in% unlist(top_genes_list),
        tool_label = case_when(
            tool == "deseq2" ~ "DESeq2",
            tool == "edger" ~ "edgeR",
            tool == "limma_trend" ~ "limma-trend",
            tool == "limma_voom" ~ "limma-voom",
            TRUE ~ tools::toTitleCase(gsub("_", " ", tool))
        )
    )

# Create label text (only for top genes)
combined_plot_data <- combined_plot_data %>%
    group_by(tool) %>%
    mutate(
        gene_label = ifelse(
            is_top_gene,
            gene_display,
            NA
        )
    ) %>%
    ungroup()

cat("  Total genes in plot:", nrow(combined_plot_data), "\n")
cat("  Genes to label:", sum(!is.na(combined_plot_data$gene_label)), "\n")
cat("  Tools in plot:", paste(unique(combined_plot_data$tool_label), collapse = ", "), "\n\n")

# ============================================================================
# CREATE FACETED PLOT
# ============================================================================

cat("Creating faceted plot...\n")

# Define colors for direction
direction_colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

# Create the plot
p <- ggplot(combined_plot_data, aes(x = log10_expr, y = log2FC)) +
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
    # Facet by tool in one row
    facet_wrap(~ tool_label, nrow = 1, scales = "free_x") +
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
        x = expression("Mean Expression (log"[10]*" baseMean)"),
        y = expression("Log"[2]*" Fold Change"),
        color = "Regulation"
    )

# ============================================================================
# SAVE PLOT
# ============================================================================

cat("Saving plot...\n")

# Calculate dimensions based on number of tools
n_tools <- length(unique(combined_plot_data$tool_label))
plot_width <- 4 * n_tools
plot_height <- 5

# Save as PDF
ggsave(
    filename = output_pdf,
    plot = p,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300
)
cat("  Saved PDF:", output_pdf, "\n")

# Save as PNG
ggsave(
    filename = output_png,
    plot = p,
    width = plot_width,
    height = plot_height,
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

    # Get gene names for top genes (handle empty results)
    top_up_gene <- tryCatch({
        gene <- deg_data %>%
            filter(direction == "Up") %>%
            arrange(desc(log2FC)) %>%
            head(1) %>%
            pull(gene_id)
        if (length(gene) == 0) "NA" else gene
    }, error = function(e) "NA")

    top_down_gene <- tryCatch({
        gene <- deg_data %>%
            filter(direction == "Down") %>%
            arrange(log2FC) %>%
            head(1) %>%
            pull(gene_id)
        if (length(gene) == 0) "NA" else gene
    }, error = function(e) "NA")

    # Convert to gene names if mapping available
    if (!is.null(gene_name_lookup)) {
        top_up_gene <- ifelse(top_up_gene %in% names(gene_name_lookup),
                             gene_name_lookup[top_up_gene], top_up_gene)
        top_down_gene <- ifelse(top_down_gene %in% names(gene_name_lookup),
                               gene_name_lookup[top_down_gene], top_down_gene)
    }

    summary_table <- rbind(summary_table, data.frame(
        Tool = case_when(
            tool == "deseq2" ~ "DESeq2",
            tool == "edger" ~ "edgeR",
            tool == "limma_trend" ~ "limma-trend",
            tool == "limma_voom" ~ "limma-voom",
            TRUE ~ tools::toTitleCase(gsub("_", " ", tool))
        ),
        Total_DEGs = sum(deg_data$significant),
        Up_regulated = sum(deg_data$direction == "Up"),
        Down_regulated = sum(deg_data$direction == "Down"),
        Top_up_gene = top_up_gene,
        Top_down_gene = top_down_gene,
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
