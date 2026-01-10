#!/usr/bin/env Rscript

# DEG Intersection and Union Heatmap Visualization
# Creates comprehensive heatmaps for DEGs identified by:
# - Intersection: DEGs significant in ALL tools (DESeq2, edgeR, limma-trend, limma-voom)
# - Union: DEGs significant in AT LEAST ONE tool
#
# Features:
# - Uses corrected TPM matrices for expression values
# - Annotates genes with log2FC and -log10(padj) for each tool
# - NA values for non-significant genes in each tool
# - Annotation bar showing number of tools that identified each DEG
# - Sample clinical information as column annotations

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

start_time <- Sys.time()
cat("=== DEG Intersection/Union Heatmap Analysis ===\n")
cat("Started at:", start_time, "\n")

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(UpSetR)
})

# Get parameters
project <- snakemake@params[["project"]]
dataset <- snakemake@params[["dataset"]]
tool <- snakemake@params[["tool"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
padj_threshold <- snakemake@params[["padj_threshold"]]
top_n <- snakemake@params[["top_n"]]
discovery_sample_type <- snakemake@params[["discovery_sample_type"]]

# Input files
deseq2_file <- snakemake@input[["deseq2"]]
edger_file <- snakemake@input[["edger"]]
limma_trend_file <- snakemake@input[["limma_trend"]]
limma_voom_file <- snakemake@input[["limma_voom"]]
tpm_file <- snakemake@input[["tpm"]]
samples_file <- snakemake@input[["samples"]]
gene_name_mapping <- snakemake@input[["gene_name_map"]]

# Output files
intersection_heatmap_pdf <- snakemake@output[["intersection_heatmap"]]
union_heatmap_pdf <- snakemake@output[["union_heatmap"]]
intersection_gene_list <- snakemake@output[["intersection_gene_list"]]
union_gene_list <- snakemake@output[["union_gene_list"]]
intersection_annotation <- snakemake@output[["intersection_annotation"]]
union_annotation <- snakemake@output[["union_annotation"]]
summary_stats <- snakemake@output[["summary_stats"]]

cat("\n--- Configuration ---\n")
cat("Project:", project, "\n")
cat("Dataset:", dataset, "\n")
cat("Quantification tool:", tool, "\n")
cat("Discovery sample type:", discovery_sample_type, "\n")
cat("Log2FC threshold:", log2fc_threshold, "\n")
cat("Padj threshold:", padj_threshold, "\n")
cat("Top N genes:", top_n, "\n")

# ============================================================================
# Function to read DEG TSV files
# ============================================================================
read_deg <- function(file, tool_name) {
  cat("\nReading", tool_name, "DEGs from:", file, "\n")
  deg_data <- read.csv(file, sep='\t', check.names=FALSE, row.names=1, comment.char='#')

  # Remove rows with NA padj or log2FC
  na_count <- sum(is.na(deg_data$padj) | is.na(deg_data$log2FoldChange))
  if (na_count > 0) {
    cat("  Removed", na_count, "rows with NA values\n")
  }
  deg_data <- deg_data[!is.na(deg_data$padj) & !is.na(deg_data$log2FoldChange), ]

  cat("  Total genes:", nrow(deg_data), "\n")

  # Check required columns
  required_cols <- c("padj", "log2FoldChange")
  missing_cols <- setdiff(required_cols, colnames(deg_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", tool_name, ": ", paste(missing_cols, collapse=", "))
  }

  return(deg_data)
}

# ============================================================================
# Load DEG data from all tools
# ============================================================================
cat("\n--- Loading DEG data from all tools ---\n")
deseq2_deg <- read_deg(deseq2_file, "DESeq2")
edger_deg <- read_deg(edger_file, "edgeR")
limma_trend_deg <- read_deg(limma_trend_file, "limma-trend")
limma_voom_deg <- read_deg(limma_voom_file, "limma-voom")

# ============================================================================
# Find significant DEGs for each tool
# ============================================================================
cat("\n--- Identifying significant DEGs ---\n")
is_significant <- function(deg_data) {
  deg_data$padj < padj_threshold & abs(deg_data$log2FoldChange) >= log2fc_threshold
}

deseq2_sig <- is_significant(deseq2_deg)
edger_sig <- is_significant(edger_deg)
limma_trend_sig <- is_significant(limma_trend_deg)
limma_voom_sig <- is_significant(limma_voom_deg)

deseq2_sig_genes <- rownames(deseq2_deg)[deseq2_sig]
edger_sig_genes <- rownames(edger_deg)[edger_sig]
limma_trend_sig_genes <- rownames(limma_trend_deg)[limma_trend_sig]
limma_voom_sig_genes <- rownames(limma_voom_deg)[limma_voom_sig]

cat("DESeq2 significant:", length(deseq2_sig_genes), "genes\n")
cat("edgeR significant:", length(edger_sig_genes), "genes\n")
cat("limma-trend significant:", length(limma_trend_sig_genes), "genes\n")
cat("limma-voom significant:", length(limma_voom_sig_genes), "genes\n")

# ============================================================================
# Define INTERSECTION and UNION gene sets
# ============================================================================
cat("\n--- Computing intersection and union ---\n")

# Intersection: genes significant in ALL 4 tools
intersection_genes <- Reduce(intersect, list(
  deseq2_sig_genes,
  edger_sig_genes,
  limma_trend_sig_genes,
  limma_voom_sig_genes
))
cat("Intersection (all 4 tools):", length(intersection_genes), "genes\n")

# Union: genes significant in AT LEAST 1 tool
union_genes <- Reduce(union, list(
  deseq2_sig_genes,
  edger_sig_genes,
  limma_trend_sig_genes,
  limma_voom_sig_genes
))
cat("Union (at least 1 tool):", length(union_genes), "genes\n")

# ============================================================================
# Read TPM matrix (corrected)
# ============================================================================
cat("\n--- Reading TPM matrix ---\n")
cat("  File:", tpm_file, "\n")
tpm_data <- read.csv(tpm_file, sep='\t', row.names=1, check.names=FALSE, comment.char="#")
cat("  TPM dimensions:", nrow(tpm_data), "genes x", ncol(tpm_data), "samples\n")

# ============================================================================
# Read sample information
# ============================================================================
cat("\n--- Reading sample information ---\n")
samples_df <- read.csv(samples_file, sep='\t', comment.char="#")
cat("  Total samples in file:", nrow(samples_df), "\n")

# Filter to project samples
samples_df <- samples_df[samples_df$project_id == project, ]
cat("  Samples in project:", nrow(samples_df), "\n")

# Filter to discovery samples only
cat("  Filtering to discovery samples (", discovery_sample_type, ")...\n")
discovery_samples_df <- samples_df[samples_df$sample_type == discovery_sample_type, ]
cat("  Discovery samples:", nrow(discovery_samples_df), "\n")

# Update samples_df to only discovery samples
samples_df <- discovery_samples_df

# Check condition distribution
condition_counts <- table(samples_df$condition)
cat("  Condition distribution:\n")
for (cond in names(condition_counts)) {
  cat("    ", cond, ":", as.character(condition_counts[cond]), "samples\n")
}

# Get valid discovery samples that are in TPM matrix
discovery_sample_names <- samples_df$sample_name
valid_samples <- intersect(discovery_sample_names, colnames(tpm_data))
cat("  Valid samples in TPM matrix:", length(valid_samples), "\n")

if (length(valid_samples) == 0) {
  stop("No discovery samples found in TPM matrix!")
}

# Subset TPM to discovery samples
tpm_subset <- tpm_data[, valid_samples, drop=FALSE]
cat("  Filtered TPM matrix:", nrow(tpm_subset), "x", ncol(tpm_subset), "\n")

# Update samples_df to only include samples in expression matrix
samples_df <- samples_df[samples_df$sample_name %in% valid_samples, ]

# Reorder samples_df to match TPM column order
sample_order <- colnames(tpm_subset)
samples_df <- samples_df[match(sample_order, samples_df$sample_name), ]

# ============================================================================
# Read gene name mapping
# ============================================================================
cat("\n--- Reading gene name mapping ---\n")
gene_map_df <- read.csv(gene_name_mapping, sep='\t', header=TRUE, comment.char="#")
gene_map_df <- unique(gene_map_df[, c("gene_id", "gene_name")])
gene_name_map <- setNames(gene_map_df$gene_name, gene_map_df$gene_id)

# ============================================================================
# Function to extract annotations for a gene set
# ============================================================================
get_annotations_for_genes <- function(gene_ids, all_deg_list) {
  n_genes <- length(gene_ids)
  n_tools <- length(all_deg_list)
  tool_names <- names(all_deg_list)

  # Initialize matrices
  log2fc_matrix <- matrix(NA, nrow=n_genes, ncol=n_tools)
  neglog10padj_matrix <- matrix(NA, nrow=n_genes, ncol=n_tools)
  colnames(log2fc_matrix) <- paste0(tool_names, "_log2fc")
  colnames(neglog10padj_matrix) <- paste0(tool_names, "_neglog10padj")
  rownames(log2fc_matrix) <- gene_ids
  rownames(neglog10padj_matrix) <- gene_ids

  # Count how many tools identified each gene as significant
  tool_count <- integer(n_genes)
  names(tool_count) <- gene_ids

  # Fill in values for each tool
  for (i in seq_along(all_deg_list)) {
    tool_name <- tool_names[i]
    deg_data <- all_deg_list[[i]]

    for (j in seq_along(gene_ids)) {
      gene_id <- gene_ids[j]

      if (gene_id %in% rownames(deg_data)) {
        log2fc <- deg_data[gene_id, "log2FoldChange"]
        padj <- deg_data[gene_id, "padj"]

        # Check if significant
        is_sig <- padj < padj_threshold && abs(log2fc) >= log2fc_threshold
        if (is_sig) {
          tool_count[j] <- tool_count[j] + 1
        }

        # Always store log2FC and -log10(padj)
        log2fc_matrix[j, i] <- log2fc
        neglog10padj_matrix[j, i] <- -log10(padj)
      }
    }
  }

  return(list(
    log2fc = log2fc_matrix,
    neglog10padj = neglog10padj_matrix,
    tool_count = tool_count
  ))
}

# ============================================================================
# Function to create heatmap for a gene set
# ============================================================================
create_deg_heatmap <- function(gene_ids, gene_set_name, output_pdf,
                                gene_list_tsv, annotation_tsv,
                                tpm_matrix, samples_info, gene_name_map_vec) {

  cat("\n=== Creating heatmap for", gene_set_name, "===\n")
  cat("  Number of genes:", length(gene_ids), "\n")

  if (length(gene_ids) == 0) {
    cat("  No genes found! Skipping...\n")
    return(NULL)
  }

  # Check which genes are in TPM matrix
  genes_in_tpm <- intersect(gene_ids, rownames(tpm_matrix))
  cat("  Genes in TPM matrix:", length(genes_in_tpm), "\n")

  if (length(genes_in_tpm) == 0) {
    cat("  No genes found in TPM matrix! Skipping...\n")
    return(NULL)
  }

  # Select top N genes by variance
  top_n_genes <- min(top_n, length(genes_in_tpm))
  expr_matrix <- as.matrix(tpm_matrix[genes_in_tpm, , drop=FALSE])
  log10_tpm <- log10(expr_matrix + 1)

  gene_vars <- apply(log10_tpm, 1, var)
  top_genes <- names(sort(gene_vars, decreasing=TRUE)[1:top_n_genes])
  expr_matrix <- log10_tpm[top_genes, ]

  cat("  Selected top", length(top_genes), "genes by variance\n")

  # Z-score normalize
  expr_matrix <- t(scale(t(expr_matrix)))
  cat("  Applied Z-score normalization\n")

  # Get all DEG data for annotation extraction
  all_deg_list <- list(
    DESeq2 = deseq2_deg,
    edgeR = edger_deg,
    limma_trend = limma_trend_deg,
    limma_voom = limma_voom_deg
  )

  # Get annotations
  ann_data <- get_annotations_for_genes(top_genes, all_deg_list)

  # Create annotation data frame
  gene_annotations <- data.frame(
    n_tools = ann_data$tool_count,
    stringsAsFactors = FALSE
  )
  rownames(gene_annotations) <- top_genes

  # Add log2FC and -log10(padj) columns
  log2fc_df <- as.data.frame(ann_data$log2fc)
  neglog10padj_df <- as.data.frame(ann_data$neglog10padj)

  gene_annotations <- cbind(gene_annotations, log2fc_df, neglog10padj_df)

  # Map gene IDs to names
  gene_ids_for_names <- rownames(expr_matrix)
  gene_names <- gene_ids_for_names
  mapped <- gene_ids_for_names %in% names(gene_name_map_vec)
  gene_names[mapped] <- gene_name_map_vec[gene_ids_for_names[mapped]]
  rownames(expr_matrix) <- gene_names
  rownames(gene_annotations) <- gene_names

  # Save annotation data
  write_tsv(gene_annotations, annotation_tsv)
  cat("  Annotation data saved to:", annotation_tsv, "\n")

  # Save gene list
  gene_list_df <- data.frame(
    gene_name = gene_names,
    gene_id = gene_ids_for_names,
    n_tools = ann_data$tool_count,
    stringsAsFactors = FALSE
  )
  write_tsv(gene_list_df, gene_list_tsv)
  cat("  Gene list saved to:", gene_list_tsv, "\n")

  # Prepare sample annotations
  sample_annotations <- data.frame(
    Condition = samples_info$condition,
    stringsAsFactors = FALSE
  )
  rownames(sample_annotations) <- samples_info$sample_name

  # Add batch if available
  if ("batch" %in% colnames(samples_info)) {
    sample_annotations$Batch <- samples_info$batch
  }

  # Add sample_type if available
  if ("sample_type" %in% colnames(samples_info)) {
    sample_annotations$SampleType <- samples_info$sample_type
  }

  # Define colors
  condition_colors <- structure(
    c("#F8766D", "#00BFC4"),
    names = unique(sample_annotations$Condition)
  )

  # Create color function for n_tools (0-4)
  n_tools_colors <- colorRamp2(c(0, 1, 2, 3, 4), c("grey90", "lightblue", "blue", "darkblue", "red"))

  # Create color function for log2FC (diverging)
  log2fc_range <- range(gene_annotations[, grep("_log2fc$", colnames(gene_annotations))], na.rm=TRUE)
  log2fc_colors <- colorRamp2(c(-max(abs(log2fc_range)), 0, max(abs(log2fc_range))), c("blue", "white", "red"))

  # Create color function for -log10(padj) (sequential)
  neglog10padj_values <- gene_annotations[, grep("_neglog10padj$", colnames(gene_annotations))]
  neglog10padj_max <- max(neglog10padj_values, na.rm=TRUE)
  neglog10padj_colors <- colorRamp2(c(0, neglog10padj_max), c("white", "darkred"))

  # Create column annotation
  col_ha <- HeatmapAnnotation(
    df = sample_annotations,
    col = list(Condition = condition_colors),
    show_annotation_name = TRUE,
    annotation_name_side = "left",
    simple_anno_size = unit(0.5, "cm")
  )

  # Create row annotations
  # First: n_tools bar
  row_anno_list <- list(
    `n_tools` = gene_annotations$n_tools
  )

  # Add log2FC and -log10padj for each tool
  tool_names <- c("DESeq2", "edgeR", "limma_trend", "limma_voom")
  for (tool in tool_names) {
    log2fc_col <- paste0(tool, "_log2fc")
    neglog10padj_col <- paste0(tool, "_neglog10padj")
    row_anno_list[[log2fc_col]] <- gene_annotations[[log2fc_col]]
    row_anno_list[[neglog10padj_col]] <- gene_annotations[[neglog10padj_col]]
  }

  # Create color list for row annotations
  row_col_list <- list(`n_tools` = n_tools_colors)
  for (tool in tool_names) {
    row_col_list[[paste0(tool, "_log2fc")]] <- log2fc_colors
    row_col_list[[paste0(tool, "_neglog10padj")]] <- neglog10padj_colors
  }

  row_ha <- rowAnnotation(
    df = as.data.frame(row_anno_list),
    col = row_col_list,
    show_annotation_name = TRUE,
    annotation_name_side = "top",
    simple_anno_size = unit(0.2, "cm"),
    gp = gpar(col = "gray")
  )

  # Create main heatmap
  ht <- Heatmap(
    expr_matrix,
    name = "Z-score",
    col = colorRamp2(c(-2, 0, 2), c("#313695", "white", "#A50026")),
    top_annotation = col_ha,
    left_annotation = row_ha,
    show_row_names = nrow(expr_matrix) <= 50,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 6),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title = "Samples",
    row_title = "Genes",
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    ),
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_split = sample_annotations$Condition
  )

  # Save heatmap
  pdf(output_pdf, width = max(14, ncol(expr_matrix) * 0.4),
      height = max(10, nrow(expr_matrix) * 0.12))
  draw(ht,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       padding = unit(c(2, 2, 2, 2), "mm"))

  # Add title
  grid::grid.text(
    paste0(gene_set_name, " DEGs Heatmap (", project, ", ", tool, ")"),
    x = unit(0.5, "npc"),
    y = unit(0.98, "npc"),
    gp = grid::gpar(fontsize = 16, fontface = "bold")
  )

  grid::grid.text(
    paste0("n=", nrow(expr_matrix), " genes | log10(TPM+1) | log2FC>=", log2fc_threshold,
           ", FDR<", padj_threshold),
    x = unit(0.5, "npc"),
    y = unit(0.95, "npc"),
    gp = grid::gpar(fontsize = 12)
  )

  dev.off()

  cat("  Heatmap saved to:", output_pdf, "\n")
  cat("  Dimensions:", nrow(expr_matrix), "x", ncol(expr_matrix), "\n")

  return(list(
    n_genes = length(gene_ids),
    n_genes_plotted = nrow(expr_matrix),
    n_genes_in_tpm = length(genes_in_tpm)
  ))
}

# ============================================================================
# Create INTERSECTION heatmap
# ============================================================================
intersection_result <- create_deg_heatmap(
  gene_ids = intersection_genes,
  gene_set_name = "INTERSECTION",
  output_pdf = intersection_heatmap_pdf,
  gene_list_tsv = intersection_gene_list,
  annotation_tsv = intersection_annotation,
  tpm_matrix = tpm_subset,
  samples_info = samples_df,
  gene_name_map_vec = gene_name_map
)

# ============================================================================
# Create UNION heatmap
# ============================================================================
union_result <- create_deg_heatmap(
  gene_ids = union_genes,
  gene_set_name = "UNION",
  output_pdf = union_heatmap_pdf,
  gene_list_tsv = union_gene_list,
  annotation_tsv = union_annotation,
  tpm_matrix = tpm_subset,
  samples_info = samples_df,
  gene_name_map_vec = gene_name_map
)

# ============================================================================
# Generate summary statistics
# ============================================================================
cat("\n--- Generating summary statistics ---\n")

summary_df <- data.frame(
  metric = c(
    "Project",
    "Dataset",
    "Quantification tool",
    "Discovery sample type",
    "Log2FC threshold",
    "Padj threshold",
    "DESeq2 significant",
    "edgeR significant",
    "limma-trend significant",
    "limma-voom significant",
    "Intersection (all 4 tools)",
    "Union (at least 1 tool)",
    "Intersection genes in TPM",
    "Union genes in TPM",
    "Intersection plotted",
    "Union plotted",
    "Analysis date"
  ),
  value = c(
    project,
    dataset,
    tool,
    discovery_sample_type,
    as.character(log2fc_threshold),
    as.character(padj_threshold),
    length(deseq2_sig_genes),
    length(edger_sig_genes),
    length(limma_trend_sig_genes),
    length(limma_voom_sig_genes),
    length(intersection_genes),
    length(union_genes),
    if (!is.null(intersection_result)) intersection_result$n_genes_in_tpm else 0,
    if (!is.null(union_result)) union_result$n_genes_in_tpm else 0,
    if (!is.null(intersection_result)) intersection_result$n_genes_plotted else 0,
    if (!is.null(union_result)) union_result$n_genes_plotted else 0,
    as.character(Sys.time())
  ),
  stringsAsFactors = FALSE
)

write_tsv(summary_df, summary_stats)
cat("Summary statistics saved to:", summary_stats, "\n")

# Print summary
cat("\n=== Summary Statistics ===\n")
print(summary_df)

cat("\n=== DEG Intersection/Union Heatmap Analysis Complete ===\n")
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units="secs")
cat("Total runtime:", round(as.numeric(elapsed_time), 2), "seconds\n")
cat("Completed at:", end_time, "\n")

sink()
