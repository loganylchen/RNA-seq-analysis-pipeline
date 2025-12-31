#!/usr/bin/env Rscript

# DEG Upset Plot Script
# Compares DEGs from DESeq2, edgeR, limma-trend, and limma-voom

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(ggplot2)
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

cat("=== DEG Upset Plot Analysis ===\n")
cat("Project:", project, "\n")
cat("Log2FC threshold:", log2fc_threshold, "\n")
cat("Padj threshold:", padj_threshold, "\n")

# Function to read DEG TSV and extract significant genes
read_deg_genes <- function(file, tool_name) {
  cat("\n--- Reading", tool_name, "DEGs ---\n")

  # Read TSV file (gene_id is row name)
  deg_data <- read.csv(file, sep='\t', check.names=FALSE, row.names=1, comment.char='#')

  cat("Total genes in", tool_name, ":", nrow(deg_data), "\n")

  # Get significant DEGs
  sig_deg <- deg_data %>%
    tibble::rownames_to_column("gene_id") %>%
    filter(padj < padj_threshold, abs(log2FoldChange) >= log2fc_threshold)

  cat("Significant DEGs in", tool_name, ":", nrow(sig_deg), "\n")

  # Return gene list
  return(sig_deg$gene_id)
}

# Read DEGs from each tool
cat("\n--- Loading DEG data ---\n")
deseq2_genes <- read_deg_genes(deseq2_file, "DESeq2")
edger_genes <- read_deg_genes(edger_file, "edgeR")
limma_trend_genes <- read_deg_genes(limma_trend_file, "limma-trend")
limma_voom_genes <- read_deg_genes(limma_voom_file, "limma-voom")

# Create list of gene sets
gene_sets <- list(
  DESeq2 = deseq2_genes,
  edgeR = edger_genes,
  limma_trend = limma_trend_genes,
  limma_voom = limma_voom_genes
)

# Create a data frame for upset plot
# Each column represents a tool, each row a gene
all_genes <- unique(c(deseq2_genes, edger_genes, limma_trend_genes, limma_voom_genes))

upset_df <- data.frame(
  gene_id = all_genes,
  DESeq2 = all_genes %in% deseq2_genes,
  edgeR = all_genes %in% edger_genes,
  limma_trend = all_genes %in% limma_trend_genes,
  limma_voom = all_genes %in% limma_voom_genes,
  stringsAsFactors = FALSE
)

# Calculate counts for each combination
upset_df$n_tools <- rowSums(upset_df[, c("DESeq2", "edgeR", "limma_trend", "limma_voom")])

# Write upset data
write_tsv(upset_df, upset_data_file)
cat("\nUpset data written to:", upset_data_file, "\n")

# Calculate summary statistics
n_deseq2 <- length(deseq2_genes)
n_edger <- length(edger_genes)
n_limma_trend <- length(limma_trend_genes)
n_limma_voom <- length(limma_voom_genes)

# Intersection statistics
all_four <- Reduce(intersect, list(deseq2_genes, edger_genes, limma_trend_genes, limma_voom_genes))
three_tools <- sum(upset_df$n_tools == 3)
two_tools <- sum(upset_df$n_tools == 2)
one_tool <- sum(upset_df$n_tools == 1)

# Pairwise intersections
deseq2_edger <- intersect(deseq2_genes, edger_genes)
deseq2_limma_trend <- intersect(deseq2_genes, limma_trend_genes)
deseq2_limma_voom <- intersect(deseq2_genes, limma_voom_genes)
edger_limma_trend <- intersect(edger_genes, limma_trend_genes)
edger_limma_voom <- intersect(edger_genes, limma_voom_genes)
limma_trend_limma_voom <- intersect(limma_trend_genes, limma_voom_genes)

# Create summary text
summary_text <- c(
  "=== DEG Upset Plot Summary ===",
  "",
  paste("Project:", project),
  paste("Quantification tool:", tool),
  paste("Log2FC threshold:", log2fc_threshold),
  paste("Padj threshold:", padj_threshold),
  paste("Analysis Date:", Sys.time()),
  "",
  "--- Individual Tool Results ---",
  paste("DESeq2:", n_deseq2, "genes"),
  paste("edgeR:", n_edger, "genes"),
  paste("limma-trend:", n_limma_trend, "genes"),
  paste("limma-voom:", n_limma_voom, "genes"),
  "",
  "--- Overlap Statistics ---",
  paste("All 4 tools:", length(all_four), "genes"),
  paste("Exactly 3 tools:", three_tools, "genes"),
  paste("Exactly 2 tools:", two_tools, "genes"),
  paste("Exactly 1 tool:", one_tool, "genes"),
  "",
  "--- Pairwise Intersections ---",
  paste("DESeq2 ∩ edgeR:", length(deseq2_edger), "genes"),
  paste("DESeq2 ∩ limma-trend:", length(deseq2_limma_trend), "genes"),
  paste("DESeq2 ∩ limma-voom:", length(deseq2_limma_voom), "genes"),
  paste("edgeR ∩ limma-trend:", length(edger_limma_trend), "genes"),
  paste("edgeR ∩ limma-voom:", length(edger_limma_voom), "genes"),
  paste("limma-trend ∩ limma-voom:", length(limma_trend_limma_voom), "genes")
)

writeLines(summary_text, summary_file)
cat("\nSummary written to:", summary_file, "\n")

# Generate upset plot using UpSetR
cat("\n--- Generating upset plot ---\n")

png(upset_plot, width = 10, height = 8, units = "in", res = 300)

# Create the upset plot
upset(
  gene_sets,
  nsets = 4,
  nintersects = NA,
  order.by = "freq",
  keep.order = TRUE,
  point.size = 4,
  line.size = 1.5,
  mainbar.color = "#377EB8",
  sets.bar.color = "#4DAF4A",
  matrix.color = "black",
  sets.x.label = "Number of DEGs",
  main.bar.y.label = "Number of Genes",
  text.scale = c(2, 1.5, 1.5, 1.5, 1.5, 1.5),
  query.legend = "top",
  queries = list(
    list(query = intersects, params = list("DESeq2", "edgeR", "limma_trend", "limma_voom"),
         color = "firebrick", active = T, query.name = "All 4 tools"),
    list(query = intersects, params = list("DESeq2", "edgeR"),
         color = "steelblue", active = T, query.name = "DESeq2 ∩ edgeR")
  )
)

# Add title
title(main = paste0("DEG Overlap:  (", project, ")"),
      sub = paste("Log2FC >", log2fc_threshold, ", FDR <", padj_threshold),
      cex.main = 2, cex.sub = 1.5)

dev.off()
cat("\nUpset plot saved to:", upset_plot, "\n")

cat("\n=== DEG Upset Plot Analysis Complete ===\n")

sink()
