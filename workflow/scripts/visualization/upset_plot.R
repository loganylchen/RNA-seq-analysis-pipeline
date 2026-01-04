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

# Check for empty sets and report
cat("\n--- Checking gene sets ---\n")
for (name in names(gene_sets_8)) {
  cat(sprintf("%s: %d genes\n", name, length(gene_sets_8[[name]])))
}

# Remove empty sets and warn
empty_sets <- sapply(gene_sets_8, function(x) length(x) == 0)
if (any(empty_sets)) {
  cat("\nWARNING: The following sets are empty and will be removed:\n")
  cat(paste(names(gene_sets_8)[empty_sets], collapse="\n"), "\n")
  gene_sets_8 <- gene_sets_8[!empty_sets]
  cat(sprintf("\nContinuing with %d sets\n", length(gene_sets_8)))
}

# Check if we have any sets left
if (length(gene_sets_8) == 0) {
  stop("ERROR: All gene sets are empty. Cannot generate upset plot.")
}

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

# Generate combined upset plot with available sets
n_sets <- length(gene_sets_8)
cat("\n--- Generating upset plot with", n_sets, "sets ---\n")

png(upset_plot, width = 14, height = 10, units = "in", res = 300)

# Color palette: assign colors based on set names
set_colors <- c()
for (name in names(gene_sets_8)) {
  if (grepl("_up$", name)) {
    # Use red shades for up-regulated
    set_colors[name] <- "#D73027"
  } else if (grepl("_down$", name)) {
    # Use blue shades for down-regulated
    set_colors[name] <- "#4575B4"
  } else {
    set_colors[name] <- "gray50"
  }
}

# Build queries dynamically based on available sets
queries <- list()

# Add query for all up-regulated (if all 4 up sets exist)
up_sets <- grep("_up$", names(gene_sets_8), value = TRUE)
if (length(up_sets) == 4) {
  queries[[length(queries) + 1]] <- list(
    query = intersects,
    params = as.list(up_sets),
    color = "darkred",
    active = TRUE,
    query.name = "All 4 up-regulated"
  )
}

# Add query for all down-regulated (if all 4 down sets exist)
down_sets <- grep("_down$", names(gene_sets_8), value = TRUE)
if (length(down_sets) == 4) {
  queries[[length(queries) + 1]] <- list(
    query = intersects,
    params = as.list(down_sets),
    color = "darkblue",
    active = TRUE,
    query.name = "All 4 down-regulated"
  )
}

# Add pairwise intersection queries
if ("DESeq2_up" %in% names(gene_sets_8) && "edgeR_up" %in% names(gene_sets_8)) {
  queries[[length(queries) + 1]] <- list(
    query = intersects,
    params = list("DESeq2_up", "edgeR_up"),
    color = "steelblue",
    active = TRUE,
    query.name = "DESeq2_up ∩ edgeR_up"
  )
}

if ("DESeq2_down" %in% names(gene_sets_8) && "edgeR_down" %in% names(gene_sets_8)) {
  queries[[length(queries) + 1]] <- list(
    query = intersects,
    params = list("DESeq2_down", "edgeR_down"),
    color = "seagreen",
    active = TRUE,
    query.name = "DESeq2_down ∩ edgeR_down"
  )
}

# Call upset with dynamic parameters
upset_call <- list(
  data = quote(gene_sets_8),
  nsets = n_sets,
  nintersects = NA,
  order.by = "freq",
  keep.order = TRUE,
  point.size = 3,
  line.size = 1.2,
  sets.bar.color = set_colors,
  matrix.color = "black",
  sets.x.label = "Number of DEGs",
  text.scale = c(1.8, 1.3, 1.2, 1.2, 1.3, 1.2),
  query.legend = "top"
)

# Add queries if any exist
if (length(queries) > 0) {
  upset_call$queries <- queries
}

do.call(upset, upset_call)

# Add title
title(main = paste0("DEG Overlap: 8 Sets (", project, ")"),
      sub = paste("Up (red) & Down (blue): log2FC |>", log2fc_threshold, ", FDR <", padj_threshold),
      cex.main = 2.2, cex.sub = 1.5)

dev.off()
cat("\n8-set upset plot saved to:", upset_plot, "\n")

cat("\n=== DEG Upset Plot Analysis Complete ===\n")

sink()
