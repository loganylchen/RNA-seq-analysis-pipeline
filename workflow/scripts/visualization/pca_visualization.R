#!/usr/bin/env Rscript

# PCA Visualization using log10(TPM+1)
# Performs PCA analysis on discovery samples using TPM data
# Generates comprehensive PCA plots with PCAtools

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

start_time <- Sys.time()
cat("PCA Analysis started at:", start_time, "\n")

suppressPackageStartupMessages({
  library(PCAtools)
  library(ggplot2)
  library(cowplot)
  library(ggplotify)
  library(dplyr)
  library(readr)
})

# Get parameters
project <- snakemake@params[["project"]]
discovery_sample_type <- snakemake@params[["discovery_sample_type"]]
color_by <- snakemake@params[["color_by"]]
shape_by <- snakemake@params[["shape_by"]]
removeVar <- snakemake@params[["removeVar"]]

# Input files
tpm_file <- snakemake@input[["tpm"]]
samples_file <- snakemake@input[["samples"]]

# Output files
pca_png <- snakemake@output[["png"]]
pca_pdf <- snakemake@output[["pdf"]]
pca_data <- snakemake@output[["pca_data"]]
pca_variance <- snakemake@output[["variance"]]

cat("=== PCA Visualization using log10(TPM+1) ===\n")
cat("Project:", project, "\n")
cat("Discovery sample type:", discovery_sample_type, "\n")
cat("Color by:", color_by, "\n")
cat("Shape by:", shape_by, "\n")
cat("Remove variance threshold:", removeVar, "\n")

# Check input files exist
cat("\n--- Checking input files ---\n")
input_files <- list(
  TPM = tpm_file,
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

# Apply log10(TPM+1) transformation
cat("\n--- Applying log10(TPM+1) transformation ---\n")
expr_matrix <- as.matrix(tpm_data)
log10_tpm <- log10(expr_matrix + 1)
cat("  Applied log10(TPM+1) transformation\n")
cat("  Expression matrix dimensions:", nrow(log10_tpm), "x", ncol(log10_tpm), "\n")
cat("  Expression range:", round(min(log10_tpm), 4), "-", round(max(log10_tpm), 4), "\n")

# Read sample information
cat("\n--- Reading sample information ---\n")
cat("  File:", samples_file, "\n")
samples_df <- read.csv(samples_file, sep='\t', comment.char="#")
cat("  Total samples in file:", nrow(samples_df), "\n")
cat("  Columns:", paste(colnames(samples_df), collapse=", "), "\n")

# Filter to project samples
samples_df <- samples_df[samples_df$project == project, ]
cat("  Samples in project:", nrow(samples_df), "\n")

# Filter to discovery samples only
cat("\n  Filtering to discovery samples only...\n")
cat("  Discovery sample type:", discovery_sample_type, "\n")
all_samples <- samples_df$sample_name
discovery_samples_df <- samples_df[samples_df$sample_type == discovery_sample_type, ]
cat("  Discovery samples:", nrow(discovery_samples_df), "\n")
cat("  Excluded samples:", nrow(samples_df) - nrow(discovery_samples_df), "\n")

# Update samples_df to only discovery samples
samples_df <- discovery_samples_df

# Check condition distribution for discovery samples
if ("condition" %in% colnames(samples_df)) {
  condition_counts <- table(samples_df$condition)
  cat("  Discovery condition distribution:\n")
  for (cond in names(condition_counts)) {
    cat("    ", cond, ":", as.character(condition_counts[cond]), "samples\n")
  }
}

# Filter expression matrix to only discovery samples
discovery_sample_names <- samples_df$sample_name
cat("\n  Filtering expression matrix to discovery samples...\n")
cat("  Original expression matrix:", nrow(log10_tpm), "x", ncol(log10_tpm), "\n")

# Check which discovery samples are in the expression matrix
valid_discovery_samples <- intersect(discovery_sample_names, colnames(log10_tpm))
cat("  Discovery samples in expression matrix:", length(valid_discovery_samples), "\n")

if (length(valid_discovery_samples) == 0) {
  stop("No discovery samples found in expression matrix!")
}

# Subset expression matrix to discovery samples
log10_tpm <- log10_tpm[, valid_discovery_samples, drop=FALSE]
cat("  Filtered expression matrix:", nrow(log10_tpm), "x", ncol(log10_tpm), "\n")

# Update samples_df to only include samples that are in the expression matrix
samples_df <- samples_df[samples_df$sample_name %in% valid_discovery_samples, ]
cat("  Final samples for PCA:", nrow(samples_df), "\n")

# Reorder samples_df to match expression matrix column order
sample_order <- colnames(log10_tpm)
samples_df <- samples_df[match(sample_order, samples_df$sample_name), ]

# Prepare metadata for PCAtools
cat("\n--- Preparing metadata for PCA ---\n")
# Set row names to sample_name for PCAtools
coldata <- samples_df
rownames(coldata) <- coldata$sample_name
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

# Filter low variance genes
cat("\n--- Filtering low variance genes ---\n")
cat("  Original genes:", nrow(log10_tpm), "\n")
gene_vars <- apply(log10_tpm, 1, var)
cat("  Variance range:", round(min(gene_vars), 6), "-", round(max(gene_vars), 6), "\n")
cat("  Genes with zero variance:", sum(gene_vars == 0), "\n")

# Keep only genes with variance > 0
log10_tpm <- log10_tpm[gene_vars > 0, , drop=FALSE]
cat("  Genes after variance filter:", nrow(log10_tpm), "\n")

# Perform PCA
cat("\n--- Performing PCA ---\n")
cat("  Input dimensions:", nrow(log10_tpm), "genes x", ncol(log10_tpm), "samples\n")
cat("  Remove variance threshold:", removeVar, "\n")

p <- pca(log10_tpm, metadata = coldata, removeVar = removeVar)
cat("  PCA completed\n")
cat("  Number of PCs:", length(p$components), "\n")
cat("  Variance explained by PC1:", round(p$variance[1] * 100, 2), "%\n")
cat("  Variance explained by PC2:", round(p$variance[2] * 100, 2), "%\n")
cat("  Cumulative variance PC1+PC2:", round(sum(p$variance[1:2]) * 100, 2), "%\n")

# Calculate Horn's parallel analysis and elbow point
cat("\n--- Calculating optimal number of PCs ---\n")
horn <- parallelPCA(log10_tpm)
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
