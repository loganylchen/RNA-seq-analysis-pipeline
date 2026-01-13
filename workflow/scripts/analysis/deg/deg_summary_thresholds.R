#!/usr/bin/env Rscript

# DEG Summary Across Multiple Thresholds
# Generates comprehensive summaries of DEG results across different
# padj and log2FC thresholds for all DEG methods (DESeq2, edgeR, limma-trend, limma-voom)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(gridExtra)
})

# Get parameters from Snakemake
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("=== DEG Summary Across Multiple Thresholds ===\n")

# Input files
deseq2_file <- snakemake@input[["deseq2"]]
edger_file <- snakemake@input[["edger"]]
limma_trend_file <- snakemake@input[["limma_trend"]]
limma_voom_file <- snakemake@input[["limma_voom"]]

# Output files
summary_table <- snakemake@output[["summary_table"]]
upset_data <- snakemake@output[["upset_data"]]
heatmap_plot <- snakemake@output[["heatmap"]]
comparison_plot <- snakemake@output[["comparison_plot"]]

# Parameters
project <- snakemake@params[["project"]]
dataset <- snakemake@params[["dataset"]]
tool <- snakemake@params[["tool"]]

cat("\nProject:", project)
cat("\nDataset:", dataset)
cat("\nQuantification tool:", tool)
cat("\n")

# ============================================================================
# Define thresholds to test
# ============================================================================
padj_thresholds <- c(0.05, 0.01, 0.001, 0.0001)
log2fc_thresholds <- c(1, 1.2, 1.5, 2)

cat("Thresholds to test:\n")
cat("  padj:", paste(padj_thresholds, collapse = ", "), "\n")
cat("  log2FC:", paste(log2fc_thresholds, collapse = ", "), "\n\n")

# ============================================================================
# Function to read DEG file
# ============================================================================
read_deg_file <- function(file_path, method_name) {
  cat("Reading", method_name, "DEGs from:", file_path, "\n")

  deg_data <- read.csv(file_path, sep='\t', header=TRUE, row.names=1,
                       check.names=FALSE, comment.char="#")

  # Remove rows with NA in critical columns
  na_before <- nrow(deg_data)
  deg_data <- deg_data[!is.na(deg_data$padj) & !is.na(deg_data$log2FoldChange), ]
  na_removed <- na_before - nrow(deg_data)

  cat("  Total genes:", nrow(deg_data), "\n")
  if (na_removed > 0) {
    cat("  Removed", na_removed, "rows with NA values\n")
  }

  return(deg_data)
}

# ============================================================================
# Function to count DEGs at different thresholds
# ============================================================================
count_degs <- function(deg_data, padj_thresh, log2fc_thresh) {
  # Count upregulated
  up <- sum(deg_data$padj < padj_thresh & deg_data$log2FoldChange >= log2fc_thresh, na.rm=TRUE)

  # Count downregulated
  down <- sum(deg_data$padj < padj_thresh & deg_data$log2FoldChange <= -log2fc_thresh, na.rm=TRUE)

  # Count total
  total <- up + down

  return(c(up = up, down = down, total = total))
}

# ============================================================================
# Function to get DEG gene lists at different thresholds
# ============================================================================
get_deg_genes <- function(deg_data, padj_thresh, log2fc_thresh) {
  # Upregulated genes
  up_genes <- rownames(deg_data)[deg_data$padj < padj_thresh &
                                deg_data$log2FoldChange >= log2fc_thresh]

  # Downregulated genes
  down_genes <- rownames(deg_data)[deg_data$padj < padj_thresh &
                                  deg_data$log2FoldChange <= -log2fc_thresh]

  return(list(up = up_genes, down = down_genes,
              up_count = length(up_genes), down_count = length(down_genes)))
}

# ============================================================================
# Load all DEG data
# ============================================================================
cat("\n=== Loading DEG Data ===\n")
deseq2_deg <- read_deg_file(deseq2_file, "DESeq2")
edger_deg <- read_deg_file(edger_file, "edgeR")
limma_trend_deg <- read_deg_file(limma_trend_file, "limma-trend")
limma_voom_deg <- read_deg_file(limma_voom_file, "limma-voom")

# ============================================================================
# Generate summary table across all thresholds
# ============================================================================
cat("\n=== Generating Summary Table ===\n")

summary_list <- list()

for (padj_thresh in padj_thresholds) {
  for (log2fc_thresh in log2fc_thresholds) {
    cat("\nThreshold: padj <", padj_thresh, ", |log2FC| >=", log2fc_thresh, "\n")

    # Count DEGs for each method
    deseq2_counts <- count_degs(deseq2_deg, padj_thresh, log2fc_thresh)
    edger_counts <- count_degs(edger_deg, padj_thresh, log2fc_thresh)
    limma_trend_counts <- count_degs(limma_trend_deg, padj_thresh, log2fc_thresh)
    limma_voom_counts <- count_degs(limma_voom_deg, padj_thresh, log2fc_thresh)

    summary_list[[length(summary_list) + 1]] <- data.frame(
      padj_threshold = padj_thresh,
      log2fc_threshold = log2fc_thresh,
      deseq2_up = deseq2_counts["up"],
      deseq2_down = deseq2_counts["down"],
      deseq2_total = deseq2_counts["total"],
      edger_up = edger_counts["up"],
      edger_down = edger_counts["down"],
      edger_total = edger_counts["total"],
      limma_trend_up = limma_trend_counts["up"],
      limma_trend_down = limma_trend_counts["down"],
      limma_trend_total = limma_trend_counts["total"],
      limma_voom_up = limma_voom_counts["up"],
      limma_voom_down = limma_voom_counts["down"],
      limma_voom_total = limma_voom_counts["total"],
      stringsAsFactors = FALSE
    )

    # Print to console
    cat("  DESeq2:  ", deseq2_counts["total"], " (up:", deseq2_counts["up"],
        ", down:", deseq2_counts["down"], ")\n")
    cat("  edgeR:   ", edger_counts["total"], " (up:", edger_counts["up"],
        ", down:", edger_counts["down"], ")\n")
    cat("  limma-trend:", limma_trend_counts["total"], " (up:", limma_trend_counts["up"],
        ", down:", limma_trend_counts["down"], ")\n")
    cat("  limma-voom:", limma_voom_counts["total"], " (up:", limma_voom_counts["up"],
        ", down:", limma_voom_counts["down"], ")\n")
  }
}

summary_df <- bind_rows(summary_list)

# Save summary table
write_tsv(summary_df, summary_table)
cat("\nSummary table saved to:", summary_table, "\n")

# ============================================================================
# Generate upset-style data for overlap analysis
# ============================================================================
cat("\n=== Generating Overlap Analysis ===\n")

# Use default thresholds (padj < 0.05, |log2FC| >= 1)
default_padj <- 0.05
default_log2fc <- 1

cat("Using default thresholds: padj <", default_padj, ", |log2FC| >=", default_log2fc, "\n")

# Get DEG genes for each method
deseq2_genes <- get_deg_genes(deseq2_deg, default_padj, default_log2fc)
edger_genes <- get_deg_genes(edger_deg, default_padj, default_log2fc)
limma_trend_genes <- get_deg_genes(limma_trend_deg, default_padj, default_log2fc)
limma_voom_genes <- get_deg_genes(limma_voom_deg, default_padj, default_log2fc)

# Calculate overlaps
all_up <- Reduce(union, list(
  deseq2_genes$up,
  edger_genes$up,
  limma_trend_genes$up,
  limma_voom_genes$up
))

all_down <- Reduce(union, list(
  deseq2_genes$down,
  edger_genes$down,
  limma_trend_genes$down,
  limma_voom_genes$down
))

# Count how many tools identify each gene
up_tool_count <- sapply(all_up, function(gene) {
  count <- 0
  if (gene %in% deseq2_genes$up) count <- count + 1
  if (gene %in% edger_genes$up) count <- count + 1
  if (gene %in% limma_trend_genes$up) count <- count + 1
  if (gene %in% limma_voom_genes$up) count <- count + 1
  return(count)
})

down_tool_count <- sapply(all_down, function(gene) {
  count <- 0
  if (gene %in% deseq2_genes$down) count <- count + 1
  if (gene %in% edger_genes$down) count <- count + 1
  if (gene %in% limma_trend_genes$down) count <- count + 1
  if (gene %in% limma_voom_genes$down) count <- count + 1
  return(count)
})

# Create upset data
upset_df <- data.frame(
  gene = c(all_up, all_down),
  direction = c(rep("up", length(all_up)), rep("down", length(all_down))),
  n_tools = c(up_tool_count, down_tool_count),
  deseq2 = c(all_up %in% deseq2_genes$up, all_down %in% deseq2_genes$down),
  edger = c(all_up %in% edger_genes$up, all_down %in% edger_genes$down),
  limma_trend = c(all_up %in% limma_trend_genes$up, all_down %in% limma_trend_genes$down),
  limma_voom = c(all_up %in% limma_voom_genes$up, all_down %in% limma_voom_genes$down),
  stringsAsFactors = FALSE
)

# Save upset data
write_tsv(upset_df, upset_data)
cat("Upset data saved to:", upset_data, "\n")

# Print overlap statistics
cat("\nOverlap statistics (padj <", default_padj, ", |log2FC| >=", default_log2fc, "):\n")
cat("  Upregulated genes identified by:\n")
cat("    All 4 tools:", sum(up_tool_count == 4), "\n")
cat("    3 tools:     ", sum(up_tool_count == 3), "\n")
cat("    2 tools:     ", sum(up_tool_count == 2), "\n")
cat("    1 tool:      ", sum(up_tool_count == 1), "\n")
cat("  Downregulated genes identified by:\n")
cat("    All 4 tools:", sum(down_tool_count == 4), "\n")
cat("    3 tools:     ", sum(down_tool_count == 3), "\n")
cat("    2 tools:     ", sum(down_tool_count == 2), "\n")
cat("    1 tool:      ", sum(down_tool_count == 1), "\n")

# ============================================================================
# Create heatmap of DEG counts across thresholds
# ============================================================================
cat("\n=== Creating Heatmap ===\n")

# Prepare matrix for heatmap
# Rows: log2FC thresholds, Columns: padj thresholds
# We'll create separate matrices for each method

heatmap_data_list <- list()

for (method in c("deseq2", "edger", "limma_trend", "limma_voom")) {
  method_data <- summary_df[, c("padj_threshold", "log2fc_threshold",
                                  paste0(method, "_total"))]
  colnames(method_data)[3] <- "n_degs"

  # Pivot to matrix format
  mat <- method_data %>%
    select(padj_threshold, log2fc_threshold, n_degs) %>%
    tidyr::pivot_wider(names_from = padj_threshold,
                       values_from = n_degs) %>%
    column_to_rownames(var = "log2fc_threshold") %>%
    as.matrix()

  # Convert log2fc_threshold to character for row names
  rownames(mat) <- paste0("FC", rownames(mat))

  heatmap_data_list[[method]] <- mat
}

# Combine all methods
combined_heatmap_data <- do.call(cbind, heatmap_data_list)
colnames(combined_heatmap_data) <- c(
  rep("DESeq2", length(padj_thresholds)),
  rep("edgeR", length(padj_thresholds)),
  rep("limma-trend", length(padj_thresholds)),
  rep("limma-voom", length(padj_thresholds))
)

# Create annotation for columns
col_anno <- data.frame(
  Method = rep(c("DESeq2", "edgeR", "limma-trend", "limma-voom"),
               each = length(padj_thresholds)),
  padj = rep(padj_thresholds, 4)
)
rownames(col_anno) <- colnames(combined_heatmap_data)

# Create heatmap
pdf(heatmap_plot, width = 12, height = 8)
pheatmap(combined_heatmap_data,
         annotation_col = col_anno,
         annotation_names_col = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         display_numbers = TRUE,
         number_format = "%.0f",
         main = paste0("DEG Counts Across Thresholds\n", project, " - ", dataset, " (", tool, ")"),
         fontsize = 10,
         fontsize_number = 8,
         angle_col = "45",
         color = colorRampPalette(c("white", "yellow", "orange", "red"))(100))
dev.off()

cat("Heatmap saved to:", heatmap_plot, "\n")

# ============================================================================
# Create comparison plot
# ============================================================================
cat("\n=== Creating Comparison Plot ===\n")

# Prepare data for plotting
plot_data <- summary_df %>%
  mutate(
    padj_threshold = factor(padj_threshold),
    log2fc_threshold = factor(log2fc_threshold)
  ) %>%
  select(padj_threshold, log2fc_threshold,
         deseq2_total, edger_total, limma_trend_total, limma_voom_total) %>%
  pivot_longer(
    cols = c(deseq2_total, edger_total, limma_trend_total, limma_voom_total),
    names_to = "method",
    values_to = "n_degs"
  ) %>%
  mutate(method = factor(method,
                         levels = c("deseq2_total", "edger_total",
                                   "limma_trend_total", "limma_voom_total"),
                         labels = c("DESeq2", "edgeR", "limma-trend", "limma-voom")))

# Create grouped bar plot
p <- ggplot(plot_data, aes(x = log2fc_threshold, y = n_degs, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
  facet_wrap(~ padj_threshold, scales = "free_y",
             labeller = labeller(padj_threshold = function(x) paste("padj <", x))) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = paste("DEG Counts Across Thresholds", project, "-", dataset, "(", tool, ")"),
    x = "log2FC threshold",
    y = "Number of DEGs",
    fill = "Method"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    strip.text = element_text(face = "bold")
  )

ggsave(comparison_plot, p, width = 14, height = 10, dpi = 300)
cat("Comparison plot saved to:", comparison_plot, "\n")

# ============================================================================
# Print final summary
# ============================================================================
cat("\n=== Summary Statistics ===\n")
cat("Total DEGs per method (padj < 0.05, |log2FC| >= 1):\n")
cat("  DESeq2:     ", deseq2_genes$up_count + deseq2_genes$down_count, "\n")
cat("  edgeR:      ", edger_genes$up_count + edger_genes$down_count, "\n")
cat("  limma-trend:", limma_trend_genes$up_count + limma_trend_genes$down_count, "\n")
cat("  limma-voom: ", limma_voom_genes$up_count + limma_voom_genes$down_count, "\n")

cat("\nIntersection (all 4 tools):\n")
cat("  Up:   ", sum(up_tool_count == 4), "\n")
cat("  Down: ", sum(down_tool_count == 4), "\n")
cat("  Total:", sum(up_tool_count == 4) + sum(down_tool_count == 4), "\n")

cat("\n=== DEG Summary Complete ===\n")

sink()
sink(type="message")
