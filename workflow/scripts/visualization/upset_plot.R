#!/usr/bin/env Rscript

# DEG Upset Plot Script
# Compares DEGs from DESeq2, edgeR, limma-trend, and limma-voom
# Creates 8 sets: up-regulated and down-regulated for each tool

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(UpSetR)
  library(ggplot2)
  library(gridExtra)
  library(RColorBrewer)
})

# Get parameters
project <- snakemake@params[["project"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
padj_threshold <- snakemake@params[["padj_threshold"]]

# Input files
deseq2_file <- snakemake@input[["deseq2"]]
edger_file <- snakemake@input[["edger"]]
limma_trend_file <- snakemake@input[["limma_trend"]]
limma_voom_file <- snakemake@input[["limma_voom"]]

# Output files
upset_plot <- snakemake@output[["upset_plot"]]
upset_data_file <- snakemake@output[["upset_data"]]
summary_file <- snakemake@output[["summary"]]

cat("=== DEG Upset Plot Analysis (8 sets) ===\n")
cat("Project:", project, "\n")
cat("Log2FC threshold:", log2fc_threshold, "\n")
cat("Padj threshold:", padj_threshold, "\n")

# Function to read DEG TSV and extract significant genes by direction
read_deg_genes_by_direction <- function(file, tool_name) {
  cat("\n--- Reading", tool_name, "DEGs ---\n")

  # Read TSV file (gene_id is row name)
  deg_data <- read.csv(file, sep='\t', check.names=FALSE, row.names=1, comment.char='#')

  cat("Total genes in", tool_name, ":", nrow(deg_data), "\n")

  # Get significant DEGs and separate by direction
  sig_deg <- deg_data %>%
    tibble::rownames_to_column("gene_id") %>%
    filter(padj < padj_threshold)

  # Up-regulated (log2FC >= threshold)
  up_genes <- sig_deg %>%
    filter(log2FoldChange >= log2fc_threshold) %>%
    pull(gene_id)

  # Down-regulated (log2FC <= -threshold)
  down_genes <- sig_deg %>%
    filter(log2FoldChange <= -log2fc_threshold) %>%
    pull(gene_id)

  cat("Up-regulated in", tool_name, ":", length(up_genes), "\n")
  cat("Down-regulated in", tool_name, ":", length(down_genes), "\n")

  return(list(up = up_genes, down = down_genes))
}

# Read DEGs from each tool
cat("\n--- Loading DEG data ---\n")
deseq2_genes <- read_deg_genes_by_direction(deseq2_file, "DESeq2")
edger_genes <- read_deg_genes_by_direction(edger_file, "edgeR")
limma_trend_genes <- read_deg_genes_by_direction(limma_trend_file, "limma-trend")
limma_voom_genes <- read_deg_genes_by_direction(limma_voom_file, "limma-voom")

# Create 8 gene sets (4 tools x 2 directions)
# Use naming like "DESeq2_up", "DESeq2_down", etc.
gene_sets_8 <- list(
  `DESeq2_up` = deseq2_genes$up,
  `edgeR_up` = edger_genes$up,
  `limma_trend_up` = limma_trend_genes$up,
  `limma_voom_up` = limma_voom_genes$up,
  `DESeq2_down` = deseq2_genes$down,
  `edgeR_down` = edger_genes$down,
  `limma_trend_down` = limma_trend_genes$down,
  `limma_voom_down` = limma_voom_genes$down
)

# Write upset data with all 8 sets
all_genes <- unique(c(deseq2_genes$up, edger_genes$up, limma_trend_genes$up, limma_voom_genes$up,
                    deseq2_genes$down, edger_genes$down, limma_trend_genes$down, limma_voom_genes$down))

upset_df <- data.frame(
  gene_id = all_genes,
  DESeq2_up = all_genes %in% deseq2_genes$up,
  edgeR_up = all_genes %in% edger_genes$up,
  limma_trend_up = all_genes %in% limma_trend_genes$up,
  limma_voom_up = all_genes %in% limma_voom_genes$up,
  DESeq2_down = all_genes %in% deseq2_genes$down,
  edgeR_down = all_genes %in% edger_genes$down,
  limma_trend_down = all_genes %in% limma_trend_genes$down,
  limma_voom_down = all_genes %in% limma_voom_genes$down,
  stringsAsFactors = FALSE
)

write_tsv(upset_df, upset_data_file)
cat("\nUpset data written to:", upset_data_file, "\n")

# Calculate summary statistics
n_up_deseq2 <- length(deseq2_genes$up)
n_up_edger <- length(edger_genes$up)
n_up_limma_trend <- length(limma_trend_genes$up)
n_up_limma_voom <- length(limma_voom_genes$up)

n_down_deseq2 <- length(deseq2_genes$down)
n_down_edger <- length(edger_genes$down)
n_down_limma_trend <- length(limma_trend_genes$down)
n_down_limma_voom <- length(limma_voom_genes$down)

# Intersection of all 8 sets
all_eight <- Reduce(intersect, list(deseq2_genes$up, edger_genes$up, limma_trend_genes$up, limma_voom_genes$up,
                                  deseq2_genes$down, edger_genes$down, limma_trend_genes$down, limma_voom_genes$down))

# Intersection of all up-regulated (4 sets)
all_up_four <- Reduce(intersect, list(deseq2_genes$up, edger_genes$up,
                                     limma_trend_genes$up, limma_voom_genes$up))

# Intersection of all down-regulated (4 sets)
all_down_four <- Reduce(intersect, list(deseq2_genes$down, edger_genes$down,
                                       limma_trend_genes$down, limma_voom_genes$down))

# Create summary text
summary_text <- c(
  "=== DEG Upset Plot Summary (8 sets) ===",
  "",
  paste("Project:", project),
  paste("Log2FC threshold:", log2fc_threshold),
  paste("Padj threshold:", padj_threshold),
  paste("Analysis Date:", Sys.time()),
  "",
  "--- Up-regulated Genes ---",
  paste("DESeq2:", n_up_deseq2, "genes"),
  paste("edgeR:", n_up_edger, "genes"),
  paste("limma-trend:", n_up_limma_trend, "genes"),
  paste("limma-voom:", n_up_limma_voom, "genes"),
  "",
  paste("All 4 up-regulated tools:", length(all_up_four), "genes"),
  if (length(all_up_four) > 0) {
    paste("  Top 20:", paste(head(all_up_four, 20), collapse = ", "))
  } else {
    NULL
  },
  "",
  "--- Down-regulated Genes ---",
  paste("DESeq2:", n_down_deseq2, "genes"),
  paste("edgeR:", n_down_edger, "genes"),
  paste("limma-trend:", n_down_limma_trend, "genes"),
  paste("limma-voom:", n_down_limma_voom, "genes"),
  "",
  paste("All 4 down-regulated tools:", length(all_down_four), "genes"),
  if (length(all_down_four) > 0) {
    paste("  Top 20:", paste(head(all_down_four, 20), collapse = ", "))
  } else {
    NULL
  },
  "",
  "--- Note ---",
  "Genes cannot be both up and down regulated,",
  "so intersection of up and down sets will always be empty."
)

summary_text <- Filter(Negate(is.null), summary_text)
writeLines(summary_text, summary_file)
cat("\nSummary written to:", summary_file, "\n")

# Generate combined upset plot with 8 sets
cat("\n--- Generating 8-set upset plot ---\n")

png(upset_plot, width = 14, height = 10, units = "in", res = 300)

# Color palette: reds for up, blues for down
up_colors <- brewer.pal(4, "Reds")[2:4]
down_colors <- brewer.pal(4, "Blues")[2:4]
set_colors <- c(up_colors, down_colors)
names(set_colors) <- names(gene_sets_8)

upset(
  gene_sets_8,
  nsets = 8,
  nintersects = NA,
  order.by = "freq",
  keep.order = TRUE,
  point.size = 3,
  line.size = 1.2,
  mainbar.color = "#E41A1C",
  sets.bar.color = set_colors,
  matrix.color = "black",
  sets.x.label = "Number of DEGs",
  main.bar.y.label = "Number of Genes",
  text.scale = c(1.8, 1.3, 1.2, 1.2, 1.3, 1.2),
  query.legend = "top",
  queries = list(
    list(query = intersects, params = list("DESeq2_up", "edgeR_up", "limma_trend_up", "limma_voom_up"),
         color = "darkred", active = T, query.name = "All 4 up-regulated"),
    list(query = intersects, params = list("DESeq2_down", "edgeR_down", "limma_trend_down", "limma_voom_down"),
         color = "darkblue", active = T, query.name = "All 4 down-regulated"),
    list(query = intersects, params = list("DESeq2_up", "DESeq2_down"),
         color = "gray", active = F, query.name = "DESeq2_up ∩ DESeq2_down (empty)"),
    list(query = intersects, params = list("DESeq2_up", "edgeR_up"),
         color = "steelblue", active = T, query.name = "DESeq2_up ∩ edgeR_up"),
    list(query = intersects, params = list("DESeq2_down", "edgeR_down"),
         color = "seagreen", active = T, query.name = "DESeq2_down ∩ edgeR_down")
  )
)

# Add title
title(main = paste0("DEG Overlap: 8 Sets (", project, ")"),
      sub = paste("Up (red) & Down (blue): log2FC |>", log2fc_threshold, ", FDR <", padj_threshold),
      cex.main = 2.2, cex.sub = 1.5)

dev.off()
cat("\n8-set upset plot saved to:", upset_plot, "\n")

cat("\n=== DEG Upset Plot Analysis Complete ===\n")

sink()
