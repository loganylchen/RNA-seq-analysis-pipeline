#!/usr/bin/env Rscript

# DEG Visualization Script
# Creates:
# - Heatmap of DEGs using ComplexHeatmap
# - Scatter plot comparing discovery vs validation log2FC
# - Boxplots for genes significant in both cohorts

suppressPackageStartupMessages({
  library(dplyr)
  library(ComplexHeatmap)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(ggrepel)
  library(tidyr)
  library(readr)
  library(tibble)
})

# Get parameters from Snakemake
discovery_deg_rds <- snakemake@input[["discovery_deg_rds"]]
validation_deg_rds <- snakemake@input[["validation_deg_rds"]]
discovery_vst_rds <- snakemake@input[["discovery_vst_rds"]]
validation_vst_rds <- snakemake@input[["validation_vst_rds"]]
samples_file <- snakemake@input[["samples"]]

project <- snakemake@params[["project"]]
case_condition <- snakemake@params[["case_condition"]]
control_condition <- snakemake@params[["control_condition"]]
discovery_sample_type <- snakemake@params[["discovery_sample_type"]]
log2fc_threshold <- as.numeric(snakemake@params[["log2fc_threshold"]])
padj_threshold <- as.numeric(snakemake@params[["padj_threshold"]])
heatmap_top_n <- as.integer(snakemake@params[["heatmap_top_n"]])

# Output files
heatmap_output <- snakemake@output[["heatmap"]]
scatter_output <- snakemake@output[["scatter"]]
boxplot_output <- snakemake@output[["boxplot"]]
discovery_deg_list <- snakemake@output[["discovery_deg_list"]]
validation_deg_list <- snakemake@output[["validation_deg_list"]]
comparison_table <- snakemake@output[["comparison_table"]]

cat("=== DEG Visualization ===\n")
cat("Parameters:\n")
cat(sprintf("  log2FC threshold: %.2f\n", log2fc_threshold))
cat(sprintf("  padj threshold: %.5f\n", padj_threshold))
cat(sprintf("  Heatmap top genes: %d\n", heatmap_top_n))

# Load DEG results
discovery_deg <- readRDS(discovery_deg_rds) %>% as.data.frame()
validation_deg <- readRDS(validation_deg_rds) %>% as.data.frame()

# Load VST expression matrices
discovery_expression <- readRDS(discovery_vst_rds) %>%
  assay() %>%
  as.matrix()
validation_expression <- readRDS(validation_vst_rds) %>%
  assay() %>%
  as.matrix()

# Load sample information
sample_info <- read_tsv(samples_file, show_col_types = FALSE) %>%
  column_to_rownames(var = "sample_name")

# Create gene ID to name mapping (if available)
# For now, use gene_id as gene_name
gene_id_to_name <- data.frame(
  gene_id = rownames(discovery_deg),
  gene_name = rownames(discovery_deg),
  stringsAsFactors = FALSE
)

# ============================================================================
# 1. Write DEG lists with gene names
# ============================================================================
cat("\n--- Generating DEG lists ---\n")

discovery_deg_list_df <- discovery_deg %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(gene_id_to_name %>% distinct(gene_id, gene_name), by = "gene_id")

validation_deg_list_df <- validation_deg %>%
  tibble::rownames_to_column("gene_id") %>%
  left_join(gene_id_to_name %>% distinct(gene_id, gene_name), by = "gene_id")

write_tsv(discovery_deg_list_df, discovery_deg_list)
write_tsv(validation_deg_list_df, validation_deg_list)

# ============================================================================
# 2. Create scatter plot comparing discovery vs validation
# ============================================================================
cat("\n--- Creating discovery vs validation scatter plot ---\n")

comparison_data <- discovery_deg %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::select(gene_id, log2FoldChange_discovery = log2FoldChange, padj_discovery = padj) %>%
  inner_join(
    validation_deg %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::select(gene_id, log2FoldChange_validation = log2FoldChange, padj_validation = padj),
    by = "gene_id"
  ) %>%
  mutate(
    same_direction = sign(log2FoldChange_discovery) == sign(log2FoldChange_validation),
    significance = case_when(
      padj_discovery < padj_threshold & padj_validation < padj_threshold & same_direction ~ "Both significant (same direction)",
      padj_discovery < padj_threshold & padj_validation < padj_threshold & !same_direction ~ "Both significant (opposite direction)",
      padj_discovery < padj_threshold ~ "Discovery only",
      padj_validation < padj_threshold ~ "Validation only",
      TRUE ~ "Not significant"
    )
  ) %>%
  left_join(gene_id_to_name %>% distinct(gene_id, gene_name), by = "gene_id")

# Identify genes to label: both significant and abs(log2FC) > threshold
genes_to_label <- comparison_data %>%
  filter(
    padj_discovery < padj_threshold,
    padj_validation < padj_threshold,
    abs(log2FoldChange_discovery) > log2fc_threshold,
    abs(log2FoldChange_validation) > log2fc_threshold
  )

# Calculate axis limits
max_abs_fc <- max(abs(c(comparison_data$log2FoldChange_discovery,
                       comparison_data$log2FoldChange_validation)), na.rm = TRUE)

# Create scatter plot
p_scatter <- ggplot(comparison_data, aes(x = log2FoldChange_discovery,
                                          y = log2FoldChange_validation,
                                          color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Both significant (same direction)" = "#9370DB",
    "Both significant (opposite direction)" = "#E41A1C",
    "Discovery only" = "#377EB8",
    "Validation only" = "#4DAF4A",
    "Not significant" = "#999999"
  )) +
  geom_text_repel(
    data = genes_to_label,
    aes(label = gene_name),
    size = 3,
    fontface = "bold",
    max.overlaps = 30,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
  geom_hline(yintercept = c(-log2fc_threshold, log2fc_threshold),
             linetype = "dashed", color = "gray30", alpha = 0.7) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
             linetype = "dashed", color = "gray30", alpha = 0.7) +
  coord_fixed(ratio = 1,
              xlim = c(-max_abs_fc * 1.1, max_abs_fc * 1.1),
              ylim = c(-max_abs_fc * 1.1, max_abs_fc * 1.1)) +
  labs(
    x = sprintf("log2 Fold Change (Discovery [%s])", discovery_sample_type),
    y = "log2 Fold Change (Validation [cfRNA])",
    color = sprintf("Significance (padj < %.2g)", padj_threshold),
    title = sprintf("Comparison of log2 Fold Changes: %s", project)
  ) +
  theme_pubr() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(scatter_output, p_scatter, width = 10, height = 8, dpi = 300)

cat(sprintf("Scatter plot saved to: %s\n", scatter_output))
cat("Summary of gene categories:\n")
print(table(comparison_data$significance))

# ============================================================================
# 3. Create heatmap of top DEGs
# ============================================================================
cat("\n--- Creating DEG heatmap ---\n")

# Get discovery DEGs
discovery_deg_names <- discovery_deg %>%
  filter(abs(log2FoldChange) > log2fc_threshold, padj < padj_threshold) %>%
  rownames()

# Take top N genes by absolute log2FC for heatmap
if (length(discovery_deg_names) > 0) {
  top_n <- min(heatmap_top_n, length(discovery_deg_names))
  top_deg_names <- discovery_deg[discovery_deg_names, ] %>%
    arrange(desc(abs(log2FoldChange))) %>%
    head(top_n) %>%
    rownames()

  # Prepare heatmap matrix
  mat <- discovery_expression[top_deg_names, ]
  mat_scaled <- t(scale(t(mat)))

  # Create annotation
  ha <- HeatmapAnnotation(
    Condition = sample_info[colnames(discovery_expression), "condition"],
    Patient = sample_info[colnames(discovery_expression), "patient"],
    col = list(
      Condition = c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C")
    )
  )

  # Create heatmap
  ht <- Heatmap(
    mat_scaled,
    top_annotation = ha,
    name = "Z-score",
    show_row_names = nrow(mat_scaled) <= 50,
    show_column_names = FALSE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    use_raster = TRUE,
    raster_quality = 2,
    column_title = sprintf("Discovery Cohort (%s)", discovery_sample_type),
    row_title = sprintf("Top %d DEGs", top_n),
    heatmap_legend_param = list(title = "Z-score"),
    border = TRUE
  )

  # Save heatmap
  png(heatmap_output, width = 12, height = 10, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  cat(sprintf("Heatmap saved to: %s\n", heatmap_output))
} else {
  cat("No DEGs found for heatmap\n")
  # Create empty placeholder
  plot.new()
  text(0.5, 0.5, "No DEGs found")
  dev.copy(png, heatmap_output, width = 6, height = 4, units = "in", res = 150)
  dev.off()
}

# ============================================================================
# 4. Create boxplots for genes significant in both datasets
# ============================================================================
cat("\n--- Creating boxplots for significant genes ---\n")

both_sig_genes <- comparison_data %>%
  filter(
    padj_discovery < padj_threshold,
    padj_validation < padj_threshold,
    abs(log2FoldChange_discovery) >= log2fc_threshold,
    abs(log2FoldChange_validation) >= log2fc_threshold
  )

cat(sprintf("Found %d genes significant in both datasets\n", nrow(both_sig_genes)))

if (nrow(both_sig_genes) > 0) {
  # Prepare discovery boxplot data
  discovery_boxplot_data <- discovery_expression[both_sig_genes$gene_id, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
    left_join(gene_id_to_name %>% distinct(gene_id, gene_name), by = "gene_id") %>%
    left_join(sample_info %>% tibble::rownames_to_column("sample"), by = "sample") %>%
    mutate(dataset = "Discovery")

  # Prepare validation boxplot data
  validation_boxplot_data <- validation_expression[both_sig_genes$gene_id, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
    left_join(gene_id_to_name %>% distinct(gene_id, gene_name), by = "gene_id") %>%
    left_join(sample_info %>% tibble::rownames_to_column("sample"), by = "sample") %>%
    mutate(dataset = "Validation")

  # Combine datasets
  combined_data <- bind_rows(discovery_boxplot_data, validation_boxplot_data) %>%
    left_join(both_sig_genes %>% dplyr::select(gene_name, significance), by = "gene_name")

  # Create faceted boxplot
  # Order genes by significance and log2FC
  gene_order <- both_sig_genes %>%
    arrange(significance, desc(abs(log2FoldChange_discovery))) %>%
    pull(gene_name)

  combined_data$gene_name <- factor(combined_data$gene_name, levels = gene_order)

  # Calculate number of facets
  n_facets <- length(unique(combined_data$significance))

  p_boxplot <- ggplot(combined_data,
                      aes(x = gene_name, y = expression, fill = condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~ significance, scales = "free_x", drop = FALSE) +
    scale_fill_manual(values = c("Normal" = "#4DAF4A", "Tumor" = "#E41A1C")) +
    coord_flip() +
    labs(
      x = "",
      y = "VST Expression (Z-score)",
      fill = "Condition",
      title = sprintf("Genes Significant in Both Discovery and Validation (%s)", project)
    ) +
    theme_pubr() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 8),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # Save boxplot
  ggsave(boxplot_output, p_boxplot, width = 12, height = max(6, nrow(both_sig_genes) * 0.3), dpi = 300)

  cat(sprintf("Boxplot saved to: %s\n", boxplot_output))
} else {
  cat("No genes significant in both datasets for boxplot\n")
  # Create empty placeholder
  plot.new()
  text(0.5, 0.5, "No genes significant in both datasets")
  dev.copy(png, boxplot_output, width = 6, height = 4, units = "in", res = 150)
  dev.off()
}

# ============================================================================
# 5. Write comparison table
# ============================================================================
cat("\n--- Writing comparison table ---\n")

comparison_table_df <- comparison_data %>%
  dplyr::select(gene_name, log2FoldChange_discovery, padj_discovery,
                log2FoldChange_validation, padj_validation, significance)

write_tsv(comparison_table_df, comparison_table)
cat(sprintf("Comparison table saved to: %s\n", comparison_table))

cat("\n=== DEG Visualization Complete ===\n")
