#!/usr/bin/env Rscript

# Common DEGs ComplexHeatmap Visualization
# Visualizes expression of genes identified as DEGs by all 4 tools
# Annotates samples with clinical information and genes with binned statistics

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

start_time <- Sys.time()
cat("Analysis started at:", start_time, "\n")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
})

# Get parameters
project <- snakemake@params[["project"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
top_n <- snakemake@params[["top_n"]]

# Input files
deseq2_file <- snakemake@input[["deseq2"]]
edger_file <- snakemake@input[["edger"]]
limma_trend_file <- snakemake@input[["limma_trend"]]
limma_voom_file <- snakemake@input[["limma_voom"]]
tpm_file <- snakemake@input[["tpm"]]
samples_file <- snakemake@input[["samples"]]
gene_name_mapping <- snakemake@input[["gene_name_map"]]

# Output files
heatmap_pdf <- snakemake@output[["heatmap"]]
gene_list_tsv <- snakemake@output[["gene_list"]]
annotation_data <- snakemake@output[["annotation_data"]]

cat("=== Common DEGs ComplexHeatmap Visualization ===\n")
cat("Project:", project, "\n")
cat("Log2FC threshold:", log2fc_threshold, "\n")
cat("Padj threshold:", padj_threshold, "\n")
cat("Top N genes:", top_n, "\n")

# Check input files exist
cat("\n--- Checking input files ---\n")
input_files <- list(
  DESeq2 = deseq2_file,
  edgeR = edger_file,
  limma_trend = limma_trend_file,
  limma_voom = limma_voom_file,
  TPM = tpm_file,
  samples = samples_file,
  gene_name_map = gene_name_mapping
)

all_files_exist <- TRUE
for (name in names(input_files)) {
  file_path <- input_files[[name]]
  exists <- file.exists(file_path)
  status <- ifelse(exists, "OK", "MISSING")
  cat(sprintf("  %s: %s (%s)\n", name, file_path, status))
  if (!exists) {
    all_files_exist <- FALSE
  }
}

if (!all_files_exist) {
  stop("One or more input files are missing!")
}

# Check output directories exist or can be created
cat("\n--- Checking output directories ---\n")
output_files <- list(
  heatmap = heatmap_pdf,
  gene_list = gene_list_tsv,
  annotation_data = annotation_data
)

for (name in names(output_files)) {
  file_path <- output_files[[name]]
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    cat(sprintf("  Creating directory: %s\n", dir_path))
    dir.create(dir_path, recursive=TRUE, showWarnings=FALSE)
  }
}
cat("  All output directories OK\n")

# Function to read DEG TSV
read_deg <- function(file, tool_name) {
  cat("\nReading", tool_name, "DEGs...\n")
  cat("  File:", file, "\n")
  deg_data <- read.csv(file, sep='\t', check.names=FALSE, row.names=1, comment.char='#')

  # Remove rows with NA padj or log2FC
  na_count <- sum(is.na(deg_data$padj) | is.na(deg_data$log2FoldChange))
  if (na_count > 0) {
    cat("  Removing", na_count, "rows with NA values\n")
  }
  deg_data <- deg_data[!is.na(deg_data$padj) & !is.na(deg_data$log2FoldChange), ]

  cat("  Total genes:", nrow(deg_data), "\n")
  cat("  Columns:", paste(colnames(deg_data), collapse=", "), "\n")

  # Check required columns
  required_cols <- c("padj", "log2FoldChange")
  missing_cols <- setdiff(required_cols, colnames(deg_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse=", "))
  }

  return(deg_data)
}

# Read DEGs from all tools
cat("\n--- Loading DEG data ---\n")
deseq2_deg <- read_deg(deseq2_file, "DESeq2")
edger_deg <- read_deg(edger_file, "edgeR")
limma_trend_deg <- read_deg(limma_trend_file, "limma-trend")
limma_voom_deg <- read_deg(limma_voom_file, "limma-voom")

# Find common DEGs (significant in all 4 tools)
cat("\n--- Finding common DEGs ---\n")
deseq2_sig <- rownames(deseq2_deg)[deseq2_deg$padj < padj_threshold &
                                   abs(deseq2_deg$log2FoldChange) >= log2fc_threshold]
edger_sig <- rownames(edger_deg)[edger_deg$padj < padj_threshold &
                                 abs(edger_deg$log2FoldChange) >= log2fc_threshold]
limma_trend_sig <- rownames(limma_trend_deg)[limma_trend_deg$padj < padj_threshold &
                                             abs(limma_trend_deg$log2FoldChange) >= log2fc_threshold]
limma_voom_sig <- rownames(limma_voom_deg)[limma_voom_deg$padj < padj_threshold &
                                           abs(limma_voom_deg$log2FoldChange) >= log2fc_threshold]

cat("DESeq2 significant:", length(deseq2_sig), "\n")
cat("edgeR significant:", length(edger_sig), "\n")
cat("limma-trend significant:", length(limma_trend_sig), "\n")
cat("limma-voom significant:", length(limma_voom_sig), "\n")

# Find intersection
common_genes <- Reduce(intersect, list(deseq2_sig, edger_sig, limma_trend_sig, limma_voom_sig))
cat("\nCommon DEGs (all 4 tools):", length(common_genes), "\n")

if (length(common_genes) == 0) {
  cat("No common DEGs found! Using relaxed criteria...\n")
  # Try genes significant in at least 3 tools
  genes_list <- list(
    DESeq2 = deseq2_sig,
    edgeR = edger_sig,
    limma_trend = limma_trend_sig,
    limma_voom = limma_voom_sig
  )
  gene_counts <- table(unlist(genes_list))
  common_genes <- names(gene_counts[gene_counts >= 3])
  cat("DEGs in at least 3 tools:", length(common_genes), "\n")
}

if (length(common_genes) == 0) {
  stop("No DEGs found with the given thresholds!")
}

# Read TPM matrix
cat("\n--- Reading TPM matrix ---\n")
cat("  File:", tpm_file, "\n")
tpm_data <- read.csv(tpm_file, sep='\t', row.names=1, check.names=FALSE, comment.char="#")
cat("  TPM matrix dimensions:", nrow(tpm_data), "genes x", ncol(tpm_data), "samples\n")
cat("  Samples:", paste(colnames(tpm_data), collapse=", "), "\n")

# Check for negative or zero values
neg_values <- sum(tpm_data < 0, na.rm=TRUE)
zero_values <- sum(tpm_data == 0, na.rm=TRUE)
cat("  Negative values:", neg_values, "\n")
cat("  Zero values:", zero_values, "\n")

# Subset to common genes
common_genes_in_tpm <- intersect(common_genes, rownames(tpm_data))
cat("Common genes in TPM matrix:", length(common_genes_in_tpm), "\n")

if (length(common_genes_in_tpm) == 0) {
  stop("No common DEGs found in TPM matrix!")
}

# Select top N genes by variance (using log10(TPM+1) for variance calculation)
cat("\n--- Selecting top genes by variance ---\n")
cat("  Available common genes:", length(common_genes_in_tpm), "\n")
cat("  Requested top N:", top_n, "\n")

expr_matrix <- as.matrix(tpm_data[common_genes_in_tpm, , drop=FALSE])
cat("  Extracted expression matrix:", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")

log10_tpm <- log10(expr_matrix + 1)
cat("  Applied log10(TPM+1) transformation\n")

gene_vars <- apply(log10_tpm, 1, var)
cat("  Calculated variance for", length(gene_vars), "genes\n")
cat("  Variance range:", round(min(gene_vars), 4), "-", round(max(gene_vars), 4), "\n")

top_genes <- names(sort(gene_vars, decreasing=TRUE)[1:min(top_n, length(gene_vars))])
expr_matrix <- log10_tpm[top_genes, ]

cat("  Selected top", length(top_genes), "genes by variance\n")
cat("  Top gene variance:", round(gene_vars[top_genes[1]], 4), "\n")
cat("  Bottom gene variance:", round(gene_vars[top_genes[length(top_genes)]], 4), "\n")

# Normalize rows (genes) to Z-score
cat("\n--- Normalizing expression matrix ---\n")
expr_matrix <- t(scale(t(expr_matrix)))

cat("  Applied Z-score normalization per gene\n")
cat("  Expression matrix dimensions:", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")
cat("  Z-score range:", round(min(expr_matrix), 2), "-", round(max(expr_matrix), 2), "\n")
cat("  Mean Z-score:", round(mean(expr_matrix), 4), "(should be ~0)\n")
cat("  SD Z-score:", round(sd(expr_matrix), 4), "(should be ~1)\n")

# Map Ensembl IDs to gene names using gene_id_to_gene_name.tsv
cat("\n--- Mapping Ensembl IDs to gene names ---\n")
gene_ids <- rownames(expr_matrix)  # Keep version suffix if present

# Read gene name mapping file
gene_map_df <- read.csv(gene_name_mapping, sep='\t', header=TRUE, comment.char="#")
cat("Gene mapping file dimensions:", nrow(gene_map_df), "x", ncol(gene_map_df), "\n")

# Remove duplicates in gene_id and gene_name pair
gene_map_df <- unique(gene_map_df[, c("gene_id", "gene_name")])
cat("Unique gene_id-gene_name pairs:", nrow(gene_map_df), "\n")

# Create mapping vector
gene_name_map <- setNames(gene_map_df$gene_name, gene_map_df$gene_id)

# Map gene IDs to gene names
mapped_names <- gene_ids
mapped_names[gene_ids %in% names(gene_name_map)] <- gene_name_map[gene_ids[gene_ids %in% names(gene_name_map)]]

rownames(expr_matrix) <- mapped_names
cat("Mapped", sum(gene_ids %in% names(gene_name_map)), "gene names\n")
cat("Unmapped genes:", sum(!(gene_ids %in% names(gene_name_map))), "\n")

# Read sample information
cat("\n--- Reading sample information ---\n")
cat("  File:", samples_file, "\n")
samples_df <- read.csv(samples_file, sep='\t', comment.char="#")
cat("  Total samples in file:", nrow(samples_df), "\n")
cat("  Columns:", paste(colnames(samples_df), collapse=", "), "\n")

# Filter to project samples
samples_df <- samples_df[samples_df$project == project, ]
cat("  Samples in project:", nrow(samples_df), "\n")

# Check condition distribution
condition_counts <- table(samples_df$condition)
cat("  Condition distribution:\n")
for (cond in names(condition_counts)) {
  cat("    ", cond, ":", as.character(condition_counts[cond]), "samples\n")
}

# Reorder columns to match sample order in expression matrix
sample_order <- colnames(expr_matrix)
cat("\n  Matching samples to expression matrix...\n")
cat("  Expression matrix samples:", length(sample_order), "\n")

# Check if all samples match
missing_in_samples <- setdiff(sample_order, samples_df$sample_name)
missing_in_expr <- setdiff(samples_df$sample_name, sample_order)

if (length(missing_in_samples) > 0) {
  cat("  WARNING: Samples in expression matrix but not in samples file:\n")
  cat("    ", paste(missing_in_samples, collapse=", "), "\n")
}

if (length(missing_in_expr) > 0) {
  cat("  WARNING: Samples in samples file but not in expression matrix:\n")
  cat("    ", paste(missing_in_expr, collapse=", "), "\n")
}

samples_df <- samples_df[match(sample_order, samples_df$sample_name), ]
cat("  Successfully matched", sum(!is.na(samples_df$sample_name)), "samples\n")

# Prepare sample annotations
cat("\n--- Preparing sample annotations ---\n")
sample_annotations <- data.frame(
  Condition = samples_df$condition,
  stringsAsFactors = FALSE
)
rownames(sample_annotations) <- samples_df$sample_name

# Add batch if available
if ("batch" %in% colnames(samples_df)) {
  sample_annotations$Batch <- samples_df$batch
  cat("  Added Batch annotation\n")
}

# Add sample_type if available
if ("sample_type" %in% colnames(samples_df)) {
  sample_annotations$SampleType <- samples_df$sample_type
  cat("  Added SampleType annotation\n")
}

cat("  Final sample annotations:", ncol(sample_annotations), "columns\n")
cat("  Annotation columns:", paste(colnames(sample_annotations), collapse=", "), "\n")

# Prepare gene annotations from all 4 tools
cat("\n--- Preparing gene annotations ---\n")

# Bin log2FC values
bin_log2fc <- function(log2fc) {
  cut(abs(log2fc),
      breaks=c(0, 1, 1.2, 1.5, 2, Inf),
      labels=c("<1", "1-1.2", "1.2-1.5", "1.5-2", ">=2"),
      include.lowest=TRUE,
      right=FALSE)
}

# Bin padj values
bin_padj <- function(padj) {
  cut(-log10(padj + 1e-300),  # Use -log10 scale for better visualization
      breaks=c(0, -log10(0.05), -log10(0.01), -log10(0.001), -log10(0.0001), Inf),
      labels=c(">=0.05", "0.05-0.01", "0.01-0.001", "0.001-0.0001", "<0.0001"),
      include.lowest=TRUE,
      right=FALSE)
}

# Use original gene_ids for matching (these match the DEG file rownames)
gene_ids_for_matching <- gene_ids

# Extract annotations for each gene from each tool
get_tool_annotations <- function(deg_data, gene_ids) {
  # Create mapping from gene_id to rowname in deg_data
  matching_genes <- intersect(gene_ids, rownames(deg_data))

  log2fc_vec <- rep(NA, length(gene_ids))
  names(log2fc_vec) <- gene_ids
  log2fc_vec[matching_genes] <- deg_data[matching_genes, "log2FoldChange"]

  padj_vec <- rep(NA, length(gene_ids))
  names(padj_vec) <- gene_ids
  padj_vec[matching_genes] <- deg_data[matching_genes, "padj"]

  list(log2fc=log2fc_vec, padj=padj_vec)
}

deseq2_ann <- get_tool_annotations(deseq2_deg, gene_ids_for_matching)
edger_ann <- get_tool_annotations(edger_deg, gene_ids_for_matching)
limma_trend_ann <- get_tool_annotations(limma_trend_deg, gene_ids_for_matching)
limma_voom_ann <- get_tool_annotations(limma_voom_deg, gene_ids_for_matching)

cat("  Extracted annotations for all tools\n")
cat("  Genes with DESeq2 annotations:", sum(!is.na(deseq2_ann$log2fc)), "\n")
cat("  Genes with edgeR annotations:", sum(!is.na(edger_ann$log2fc)), "\n")
cat("  Genes with limma-trend annotations:", sum(!is.na(limma_trend_ann$log2fc)), "\n")
cat("  Genes with limma-voom annotations:", sum(!is.na(limma_voom_ann$log2fc)), "\n")

# Create gene annotation data frame
gene_annotations <- data.frame(
  DESeq2_log2FC = bin_log2fc(deseq2_ann$log2fc),
  DESeq2_padj = bin_padj(deseq2_ann$padj),
  edgeR_log2FC = bin_log2fc(edger_ann$log2fc),
  edgeR_padj = bin_padj(edger_ann$padj),
  limma_trend_log2FC = bin_log2fc(limma_trend_ann$log2fc),
  limma_trend_padj = bin_padj(limma_trend_ann$padj),
  limma_voom_log2FC = bin_log2fc(limma_voom_ann$log2fc),
  limma_voom_padj = bin_padj(limma_voom_ann$padj),
  stringsAsFactors = FALSE
)
rownames(gene_annotations) <- rownames(expr_matrix)

# Save annotation data
write_tsv(gene_annotations, annotation_data)
cat("\nAnnotation data saved to:", annotation_data, "\n")

# Save gene list (gene_id and gene_name)
gene_list_df <- data.frame(
  gene_name = rownames(expr_matrix),
  gene_id = gene_ids,
  stringsAsFactors = FALSE
)
write_tsv(gene_list_df, gene_list_tsv)
cat("Gene list saved to:", gene_list_tsv, "\n")

# Define colors for annotations
condition_colors <- structure(
  c("#F8766D", "#00BFC4"),  # Red for case, Cyan for control
  names = unique(sample_annotations$Condition)
)

# Define colors for log2FC bins (increasing intensity)
log2fc_colors <- c(
  "<1" = "#FFF7FB",
  "1-1.2" = "#FEE0D2",
  "1.2-1.5" = "#FC9272",
  "1.5-2" = "#DE2D26",
  ">=2" = "#67000D"
)

# Define colors for padj bins (increasing intensity)
padj_colors <- c(
  ">=0.05" = "#F7FBFF",
  "0.05-0.01" = "#C6DBEF",
  "0.01-0.001" = "#6BAED6",
  "0.001-0.0001" = "#2171B5",
  "<0.0001" = "#08306B"
)

cat("\n--- Generating ComplexHeatmap ---\n")
cat("  Heatmap dimensions:", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")
cat("  Number of row annotations:", 8, "(log2FC + padj for 4 tools)\n")
cat("  Number of column annotations:", ncol(sample_annotations), "\n")
cat("  Row clustering: k-means with", min(3, nrow(expr_matrix)), "clusters\n")
cat("  Column splitting: by condition\n")
cat("  PDF output:", heatmap_pdf, "\n\n")

# Create column annotation for samples
col_ha <- HeatmapAnnotation(
  df = sample_annotations,
  col = list(Condition = condition_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  simple_anno_size = unit(0.5, "cm")
)

# Create row annotation for genes (log2FC and padj from each tool)
row_ha <- rowAnnotation(
  # DESeq2
  DESeq2_log2FC = gene_annotations$DESeq2_log2FC,
  DESeq2_padj = gene_annotations$DESeq2_padj,
  # edgeR
  edgeR_log2FC = gene_annotations$edgeR_log2FC,
  edgeR_padj = gene_annotations$edgeR_padj,
  # limma-trend
  limma_trend_log2FC = gene_annotations$limma_trend_log2FC,
  limma_trend_padj = gene_annotations$limma_trend_padj,
  # limma-voom
  limma_voom_log2FC = gene_annotations$limma_voom_log2FC,
  limma_voom_padj = gene_annotations$limma_voom_padj,
  col = list(
    DESeq2_log2FC = log2fc_colors,
    DESeq2_padj = padj_colors,
    edgeR_log2FC = log2fc_colors,
    edgeR_padj = padj_colors,
    limma_trend_log2FC = log2fc_colors,
    limma_trend_padj = padj_colors,
    limma_voom_log2FC = log2fc_colors,
    limma_voom_padj = padj_colors
  ),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  simple_anno_size = unit(0.3, "cm"),
  gp = gpar(col = "gray")
)

# Create the main heatmap
ht <- Heatmap(
  expr_matrix,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
  top_annotation = col_ha,
  left_annotation = row_ha,
  show_row_names = ifelse(nrow(expr_matrix) <= 50, TRUE, FALSE),
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 45,
  column_title = "Samples",
  row_title = "Genes",
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  ),
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  column_split = sample_annotations$Condition,
  row_km = min(3, nrow(expr_matrix))
)

# Draw the heatmap
pdf(heatmap_pdf, width = max(12, ncol(expr_matrix) * 0.4),
    height = max(10, nrow(expr_matrix) * 0.15))
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     padding = unit(c(2, 2, 2, 2), "mm"))

# Add title
grid::grid.text(
  paste0("Common DEGs Expression Heatmap (", project, ")"),
  x = unit(0.5, "npc"),
  y = unit(0.98, "npc"),
  gp = grid::gpar(fontsize = 16, fontface = "bold")
)

# Add subtitle
grid::grid.text(
  paste0("Top ", nrow(expr_matrix), " genes | log10(TPM+1) | log2FC>=", log2fc_threshold,
         ", FDR<", padj_threshold),
  x = unit(0.5, "npc"),
  y = unit(0.95, "npc"),
  gp = grid::gpar(fontsize = 12)
)

dev.off()

cat("\n--- Heatmap generation complete ---\n")
cat("  Heatmap saved to:", heatmap_pdf, "\n")
cat("  Dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
cat("  File size:", round(file.info(heatmap_pdf)$size / 1024, 2), "KB\n")

cat("\n--- Output files ---\n")
cat("  Heatmap PDF:", heatmap_pdf, "\n")
cat("  Gene list TSV:", gene_list_tsv, "\n")
cat("  Annotation data TSV:", annotation_data, "\n")

cat("\n=== Common DEGs ComplexHeatmap Visualization Complete ===\n")
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units="secs")
cat("Total runtime:", round(as.numeric(elapsed_time), 2), "seconds\n")
cat("Analysis completed at:", end_time, "\n")

sink()
