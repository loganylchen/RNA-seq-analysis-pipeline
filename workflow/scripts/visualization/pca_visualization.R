#!/usr/bin/env Rscript

# PCA Visualization using VST transformed count matrix
# Performs PCA analysis on discovery samples using DESeq2 VST data
# Generates comprehensive PCA plots with PCAtools

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

start_time <- Sys.time()
cat("PCA Analysis started at:", start_time, "\n")

suppressPackageStartupMessages({
  library(DESeq2)
  library(PCAtools)
  library(ggplot2)
  library(cowplot)
  library(ggplotify)
  library(dplyr)
  library(readr)
})

# Setup parallelization
parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
  cat("Using", snakemake@threads, "threads for DESeq2\n")
}

# Get parameters
project <- snakemake@params[["project"]]
discovery_sample_type <- snakemake@params[["discovery_sample_type"]]
model <- snakemake@params[["model"]]
count_threshold <- snakemake@params[["count_threshold"]]
color_by <- snakemake@params[["color_by"]]
shape_by <- snakemake@params[["shape_by"]]
removeVar <- snakemake@params[["removeVar"]]

# Input files
counts_file <- snakemake@input[["counts"]]
samples_file <- snakemake@input[["samples"]]

# Output files
pca_png <- snakemake@output[["png"]]
pca_pdf <- snakemake@output[["pdf"]]
pca_data <- snakemake@output[["pca_data"]]
pca_variance <- snakemake@output[["variance"]]

cat("=== PCA Visualization using DESeq2 VST ===\n")
cat("Project:", project, "\n")
cat("Discovery sample type:", discovery_sample_type, "\n")
cat("Model formula:", model, "\n")
cat("Count threshold:", count_threshold, "\n")
cat("Color by:", color_by, "\n")
cat("Shape by:", shape_by, "\n")
cat("Remove variance threshold:", removeVar, "\n")

# Check input files exist
cat("\n--- Checking input files ---\n")
input_files <- list(
  counts = counts_file,
  samples = samples_file
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
  png = pca_png,
  pdf = pca_pdf,
  pca_data = pca_data,
  variance = pca_variance
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

# Read count matrix
cat("\n--- Reading count matrix ---\n")
cat("  File:", counts_file, "\n")
cts <- read.table(counts_file, header=TRUE, row.names="gene", check.names=FALSE, sep='\t')
cat("  Count matrix dimensions:", nrow(cts), "genes x", ncol(cts), "samples\n")
cat("  Samples:", paste(colnames(cts), collapse=", "), "\n")

# Order samples alphabetically
cts <- cts[, order(names(cts))]
cat("  Samples ordered alphabetically\n")

# Read sample information
cat("\n--- Reading sample information ---\n")
cat("  File:", samples_file, "\n")
coldata <- read.table(samples_file, header=TRUE, row.names="sample_name", check.names=FALSE, sep='\t', comment.char="#")
cat("  Total samples in file:", nrow(coldata), "\n")
cat("  Columns:", paste(colnames(coldata), collapse=", "), "\n")

# Filter to project samples
coldata <- coldata[coldata$project == project, ]
cat("  Samples in project:", nrow(coldata), "\n")

# Filter to discovery samples only
cat("\n  Filtering to discovery samples only...\n")
cat("  Discovery sample type:", discovery_sample_type, "\n")
all_samples <- rownames(coldata)
discovery_coldata <- coldata[coldata$sample_type == discovery_sample_type, , drop=FALSE]
cat("  Discovery samples:", nrow(discovery_coldata), "\n")
cat("  Excluded samples:", nrow(coldata) - nrow(discovery_coldata), "\n")

# Update coldata to only discovery samples
coldata <- discovery_coldata

# Check condition distribution for discovery samples
if ("condition" %in% colnames(coldata)) {
  condition_counts <- table(coldata$condition)
  cat("  Discovery condition distribution:\n")
  for (cond in names(condition_counts)) {
    cat("    ", cond, ":", as.character(condition_counts[cond]), "samples\n")
  }
}

# Order samples
coldata <- coldata[order(rownames(coldata)), , drop=FALSE]

# Filter count matrix to only discovery samples
discovery_sample_names <- rownames(coldata)
cat("\n  Filtering count matrix to discovery samples...\n")
cat("  Original count matrix:", nrow(cts), "x", ncol(cts), "\n")

# Check which discovery samples are in the count matrix
valid_discovery_samples <- intersect(discovery_sample_names, colnames(cts))
cat("  Discovery samples in count matrix:", length(valid_discovery_samples), "\n")

if (length(valid_discovery_samples) == 0) {
  stop("No discovery samples found in count matrix!")
}

# Subset count matrix to discovery samples
cts <- cts[, valid_discovery_samples, drop=FALSE]
cat("  Filtered count matrix:", nrow(cts), "x", ncol(cts), "\n")

# Update coldata to only include samples that are in the count matrix
coldata <- coldata[rownames(coldata) %in% valid_discovery_samples, , drop=FALSE]
cat("  Final samples for PCA:", nrow(coldata), "\n")

# Reorder coldata to match count matrix column order
sample_order <- colnames(cts)
coldata <- coldata[match(sample_order, rownames(coldata)), , drop=FALSE]

# Create DESeq2 dataset
cat("\n--- Creating DESeq2 dataset ---\n")
cat("  Count dimensions:", nrow(cts), "x", ncol(cts), "\n")
cat("  Coldata dimensions:", nrow(coldata), "x", ncol(coldata), "\n")
cat("  Design formula:", model, "\n")

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=as.formula(model))

# Remove uninformative genes (low counts)
cat("\n--- Filtering low count genes ---\n")
cat("  Original genes:", nrow(dds), "\n")
cat("  Count threshold:", count_threshold, "\n")
dds <- dds[rowSums(counts(dds)) > count_threshold, ]
cat("  Genes after filtering:", nrow(dds), "\n")

# Run DESeq2 for normalization
cat("\n--- Running DESeq2 for normalization ---\n")
cat("  Parallel:", parallel, "\n")
dds <- DESeq(dds, parallel=parallel)
cat("  DESeq2 completed\n")

# Apply VST transformation
cat("\n--- Applying VST transformation ---\n")
vst_data <- assay(vst(dds))
cat("  VST data dimensions:", nrow(vst_data), "x", ncol(vst_data), "\n")
cat("  VST range:", round(min(vst_data), 4), "-", round(max(vst_data), 4), "\n")

# Prepare metadata for PCAtools
cat("\n--- Preparing metadata for PCA ---\n")
cat("  Metadata dimensions:", nrow(coldata), "samples x", ncol(coldata), "columns\n")
cat("  Metadata columns:", paste(colnames(coldata), collapse=", "), "\n")

# Check if color_by and shape_by columns exist
if (!color_by %in% colnames(coldata)) {
  cat("  WARNING: color_by column '", color_by, "' not found in metadata\n")
  cat("  Available columns:", paste(colnames(coldata), collapse=", "), "\n")
  color_by <- colnames(coldata)[1]
  cat("  Using", color_by, "for coloring\n")
}

if (!shape_by %in% colnames(coldata)) {
  cat("  WARNING: shape_by column '", shape_by, "' not found in metadata\n")
  cat("  Available columns:", paste(colnames(coldata), collapse=", "), "\n")
  # Use first non-color_by column for shape
  shape_candidates <- setdiff(colnames(coldata), color_by)
  if (length(shape_candidates) > 0) {
    shape_by <- shape_candidates[1]
    cat("  Using", shape_by, "for shape\n")
  } else {
    shape_by <- NULL
    cat("  No alternative column found for shape\n")
  }
}

# Perform PCA
cat("\n--- Performing PCA ---\n")
cat("  Input dimensions:", nrow(vst_data), "genes x", ncol(vst_data), "samples\n")
cat("  Remove variance threshold:", removeVar, "\n")

p <- pca(vst_data, metadata = coldata, removeVar = removeVar)
cat("  PCA completed\n")
cat("  Number of PCs:", length(p$components), "\n")
cat("  Variance explained by PC1:", round(p$variance[1] * 100, 2), "%\n")
cat("  Variance explained by PC2:", round(p$variance[2] * 100, 2), "%\n")
cat("  Cumulative variance PC1+PC2:", round(sum(p$variance[1:2]) * 100, 2), "%\n")

# Calculate Horn's parallel analysis and elbow point
cat("\n--- Calculating optimal number of PCs ---\n")
horn <- parallelPCA(vst_data)
cat("  Horn's parallel analysis suggests:", horn$n, "components\n")

elbow <- findElbowPoint(p$variance)
cat("  Elbow method suggests:", elbow, "components\n")

# Save PCA data
cat("\n--- Saving PCA data ---\n")
pca_df <- as.data.frame(p$principal)
pca_df$sample_name <- rownames(pca_df)
# Merge with metadata
pca_df <- merge(pca_df, coldata, by="row.names")
rownames(pca_df) <- pca_df$Row.names
pca_df$Row.names <- NULL
write_tsv(pca_df, pca_data)
cat("  PCA data saved to:", pca_data, "\n")

# Save variance data
variance_df <- data.frame(
  PC = 1:length(p$variance),
  Variance = p$variance,
  CumulativeVariance = cumsum(p$variance)
)
write_tsv(variance_df, pca_variance)
cat("  Variance data saved to:", pca_variance, "\n")

# Generate PCA plots
cat("\n--- Generating PCA plots ---\n")

# Scree plot
cat("  Creating scree plot...\n")
vline <- as.numeric(which(cumsum(p$variance) > 80)[1])
if (is.na(vline)) vline <- length(p$variance)

pscree <- screeplot(p,
    components = getComponents(p, 1:min(20, length(p$components))),
    vline = c(horn$n, elbow)) +
    geom_label(aes(x = horn$n + 1, y = 50,
      label = 'Horn\'s', vjust = -1, size = 8)) +
    geom_label(aes(x = elbow + 1, y = 60,
      label = 'Elbow method', vjust = -1, size = 8))

# Pairs plot
cat("  Creating pairs plot...\n")
n_pairs <- min(5, length(p$components))
ppairs <- pairsplot(p, components = getComponents(p, 1:n_pairs),
triangle = TRUE, trianglelabSize = 12,
hline = 0, vline = 0,
pointSize = 0.8, gridlines.major = FALSE, gridlines.minor = FALSE,
colby = color_by,
title = '', plotaxes = FALSE,
margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
returnPlot = FALSE)

# Biplot
cat("  Creating biplot...\n")
pbiplot <- biplot(p,
  # loadings parameters
  showLoadings = TRUE,
  lengthLoadingsArrowsFactor = 1.5,
  sizeLoadingsNames = 4,
  colLoadingsNames = 'red4',
  # other parameters
  lab = NULL,
  colby = color_by,
  hline = 0, vline = 0,
  gridlines.major = FALSE, gridlines.minor = FALSE,
  pointSize = 5,
  legendLabSize = 16, legendIconSize = 8.0,
  shape = ifelse(is.null(shape_by), NULL, shape_by),
  drawConnectors = FALSE,
  title = 'PCA',
  subtitle = 'PC1 versus PC2',
  returnPlot = FALSE)

# Loadings plot
cat("  Creating loadings plot...\n")
ploadings <- plotloadings(p, rangeRetain = 0.01, labSize = 4,
  title = 'Loadings plot', axisLabSize = 12,
  subtitle = 'PC1, PC2, PC3, PC4, PC5',
  caption = 'Top 1% variables',
  shape = 24, shapeSizeRange = c(4, 8),
  col = c('limegreen', 'black', 'red3'),
  drawConnectors = FALSE,
  returnPlot = FALSE)

# Eigencorplot
cat("  Creating eigencorplot...\n")
# Get metadata columns that are numeric or factor for correlation
metavars <- colnames(coldata)
# Remove sample_name, project, sample_type from metavars
metavars <- setdiff(metavars, c("sample_name", "project", "sample_type"))
if (length(metavars) == 0) {
  metavars <- colnames(coldata)
}

n_eigencor <- max(horn$n, elbow)
if (n_eigencor > length(p$components)) {
  n_eigencor <- length(p$components)
}

peigencor <- eigencorplot(p,
  components = getComponents(p, 1:n_eigencor),
  metavars = metavars,
  col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
  cexCorval = 1.0,
  fontCorval = 2,
  posLab = 'all',
  rotLabX = 45,
  scale = FALSE,
  main = "PC clinical correlates",
  cexMain = 1.5,
  plotRsquared = FALSE,
  corFUN = 'pearson',
  corUSE = 'pairwise.complete.obs',
  signifSymbols = c('****', '***', '**', '*', ''),
  signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),
  returnPlot = FALSE)

# Combine plots
cat("  Combining plots...\n")
top_row <- plot_grid(pscree, ppairs, pbiplot,
  ncol = 3,
  labels = c('A', 'B ', 'C'),
  label_fontfamily = 'serif',
  label_fontface = 'bold',
  label_size = 22,
  align = 'h',
  rel_widths = c(1.10, 0.80, 1.10))

bottom_row <- plot_grid(ploadings,
  as.grob(peigencor),
  ncol = 2,
  labels = c('D', 'E'),
  label_fontfamily = 'serif',
  label_fontface = 'bold',
  label_size = 22,
  align = 'h',
  rel_widths = c(0.8, 1.2))

final_plot <- plot_grid(top_row, bottom_row, ncol = 1,
  rel_heights = c(1.1, 0.9))

# Save plots
cat("\n--- Saving plots ---\n")
cat("  Saving PNG:", pca_png, "\n")
ggsave(pca_png, final_plot, width=20, height=15)
cat("  Saved PNG:", round(file.info(pca_png)$size / 1024, 2), "KB\n")

cat("  Saving PDF:", pca_pdf, "\n")
ggsave(pca_pdf, final_plot, width=20, height=15)
cat("  Saved PDF:", round(file.info(pca_pdf)$size / 1024, 2), "KB\n")

cat("\n=== PCA Visualization Complete ===\n")
cat("Output files:\n")
cat("  PNG:", pca_png, "\n")
cat("  PDF:", pca_pdf, "\n")
cat("  PCA data:", pca_data, "\n")
cat("  Variance data:", pca_variance, "\n")

end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units="secs")
cat("Total runtime:", round(as.numeric(elapsed_time), 2), "seconds\n")
cat("Analysis completed at:", end_time, "\n")

sink()
