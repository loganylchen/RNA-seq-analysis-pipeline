#!/usr/bin/env Rscript

# Dark Channel Biomarker (DCB) Analysis
# Based on: "Cell-free RNA biomarkers for cancer detection"
# (Circulating Cell-free Genome Atlas study)
#
# DCB genes are defined as genes that:
# 1. Are highly expressed in tumor tissue (discovery cohort)
# 2. Are detected in cancer cfRNA but absent in normal cfRNA (validation cohort)
# 3. Reside in "low-noise" genomic regions (not detected in normal cfRNA)

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
})

# Get parameters from Snakemake
tpm_file <- snakemake@input[["tpm"]]
count_file <- snakemake@input[["count_matrix"]]
samples_file <- snakemake@params[["samples"]]
project <- snakemake@params[["project"]]
case_condition <- snakemake@params[["case_condition"]]
control_condition <- snakemake@params[["control_condition"]]
discovery_sample_type <- snakemake@params[["discovery_sample_type"]]

# DCB detection thresholds
tissue_tpm_threshold <- as.numeric(snakemake@params[["case_tpm_threshold"]])      # Min TPM in tumor tissue
tissue_fc_threshold <- as.numeric(snakemake@params[["log2fc"]])                  # Min log2FC in tissue
normal_cfrna_detection_rate <- as.numeric(snakemake@params[["control_detection_rate"]])  # Max detection rate in normal cfRNA
cancer_cfrna_detection_rate <- as.numeric(snakemake@params[["case_detection_rate"]])    # Min detection rate in cancer cfRNA
cfrna_tpm_threshold <- as.numeric(snakemake@params[["control_tpm_threshold"]])  # TPM threshold for cfRNA detection

# Output files
discovery_dcb_tsv <- snakemake@output[["discovery_dcb_tsv"]]
discovery_dcb_rds <- snakemake@output[["discovery_dcb_rds"]]
validation_dcb_tsv <- snakemake@output[["validation_dcb_tsv"]]
validation_dcb_rds <- snakemake@output[["validation_dcb_rds"]]
discovery_summary <- snakemake@output[["discovery_summary"]]
validation_summary <- snakemake@output[["validation_summary"]]
discovery_plot <- snakemake@output[["discovery_plot"]]
validation_plot <- snakemake@output[["validation_plot"]]

cat("=== Dark Channel Biomarker (DCB) Analysis ===\n")
cat("\nStudy Design:\n")
cat("  Discovery cohort: Tissue (tumor vs normal)\n")
cat("  Validation cohort: cfRNA (cancer vs normal)\n")
cat("\nDCB Criteria:\n")
cat("  1. Highly expressed in tumor tissue (vs normal tissue)\n")
cat("  2. Detected in cancer cfRNA but absent in normal cfRNA\n")
cat("\nParameters:\n")
cat(sprintf("  Tissue TPM threshold: %.2f\n", tissue_tpm_threshold))
cat(sprintf("  Tissue log2FC threshold: %.2f\n", tissue_fc_threshold))
cat(sprintf("  cfRNA TPM detection threshold: %.2f\n", cfrna_tpm_threshold))
cat(sprintf("  Normal cfRNA max detection rate: %.0f%%\n", normal_cfrna_detection_rate * 100))
cat(sprintf("  Cancer cfRNA min detection rate: %.0f%%\n", cancer_cfrna_detection_rate * 100))

# Read sample information
samples <- read_tsv(samples_file, show_col_types = FALSE) %>%
  filter(project == !!project)

# Read TPM and count matrices
tpm_data <- read_tsv(tpm_file, show_col_types = FALSE)
count_data <- read_tsv(count_file, show_col_types = FALSE)

# Extract gene names column
gene_col <- colnames(tpm_data)[1]

# Split samples by sample_type (tissue vs cfRNA)
tissue_samples <- samples %>%
  filter(sample_type == discovery_sample_type)

cfrna_samples <- samples %>%
  filter(sample_type != discovery_sample_type)

# Split tissue samples by condition (discovery cohort)
tissue_tumor_samples <- tissue_samples %>%
  filter(condition == case_condition) %>%
  pull(sample_name)

tissue_normal_samples <- tissue_samples %>%
  filter(condition == control_condition) %>%
  pull(sample_name)

# Split cfRNA samples by condition (validation cohort)
cfrna_cancer_samples <- cfrna_samples %>%
  filter(condition == case_condition) %>%
  pull(sample_name)

cfrna_normal_samples <- cfrna_samples %>%
  filter(condition == control_condition) %>%
  pull(sample_name)

cat(sprintf("\n=== Sample Summary ===\n"))
cat(sprintf("Discovery cohort (tissue): %d tumor, %d normal\n",
            length(tissue_tumor_samples), length(tissue_normal_samples)))
cat(sprintf("Validation cohort (cfRNA): %d cancer, %d normal\n",
            length(cfrna_cancer_samples), length(cfrna_normal_samples)))

# Get sample columns that exist in the data
tissue_tumor_cols <- intersect(tissue_tumor_samples, colnames(tpm_data))
tissue_normal_cols <- intersect(tissue_normal_samples, colnames(tpm_data))
cfrna_cancer_cols <- intersect(cfrna_cancer_samples, colnames(tpm_data))
cfrna_normal_cols <- intersect(cfrna_normal_samples, colnames(tpm_data))

cat(sprintf("\nAvailable samples in TPM data:\n"))
cat(sprintf("  Tissue tumor: %d\n", length(tissue_tumor_cols)))
cat(sprintf("  Tissue normal: %d\n", length(tissue_normal_cols)))
cat(sprintf("  cfRNA cancer: %d\n", length(cfrna_cancer_cols)))
cat(sprintf("  cfRNA normal: %d\n", length(cfrna_normal_cols)))

# Check if we have sufficient samples
if (length(tissue_tumor_cols) == 0 || length(tissue_normal_cols) == 0) {
  stop("Insufficient tissue samples for discovery analysis")
}
if (length(cfrna_cancer_cols) == 0 || length(cfrna_normal_cols) == 0) {
  stop("Insufficient cfRNA samples for validation analysis")
}

#############################################################################
# Step 1: Discovery Analysis (Tissue)
# Find genes highly expressed in tumor tissue vs normal tissue
#############################################################################
cat("\n=== Step 1: Discovery Analysis (Tissue) ===\n")

# Calculate tissue expression statistics
tissue_tumor_mean_tpm <- rowMeans(tpm_data[, tissue_tumor_cols, drop = FALSE], na.rm = TRUE)
tissue_normal_mean_tpm <- rowMeans(tpm_data[, tissue_normal_cols, drop = FALSE], na.rm = TRUE)
tissue_log2fc <- log2((tissue_tumor_mean_tpm + 0.01) / (tissue_normal_mean_tpm + 0.01))

# Identify tissue-upregulated genes (candidate DCBs)
discovery_genes <- tpm_data %>%
  select(1) %>%
  mutate(
    gene = !!gene_col,
    tissue_tumor_mean_tpm = tissue_tumor_mean_tpm,
    tissue_normal_mean_tpm = tissue_normal_mean_tpm,
    tissue_log2fc = tissue_log2fc,
    is_tissue_upregulated = (tissue_tumor_mean_tpm >= tissue_tpm_threshold) &
                             (abs(tissue_log2fc) >= tissue_fc_threshold)
  )

# Summary statistics
tissue_upregulated_count <- sum(discovery_genes$is_tissue_upregulated, na.rm = TRUE)
cat(sprintf("Tissue-upregulated genes (TPM >= %.2f, |log2FC| >= %.2f): %d\n",
            tissue_tpm_threshold, tissue_fc_threshold, tissue_upregulated_count))

#############################################################################
# Step 2: Validation Analysis (cfRNA)
# Filter tissue-upregulated genes by cfRNA detection pattern
#############################################################################
cat("\n=== Step 2: Validation Analysis (cfRNA) ===\n")

# Calculate cfRNA detection statistics
# A gene is "detected" in cfRNA if TPM > threshold
cfrna_detection_matrix <- tpm_data %>%
  select(all_of(c(cfrna_cancer_cols, cfrna_normal_cols))) > cfrna_tpm_threshold

# Detection rate in normal cfRNA (should be LOW for DCB genes)
cfrna_normal_detection_rate <- rowSums(cfrna_detection_matrix[, cfrna_normal_cols, drop = FALSE], na.rm = TRUE) /
  length(cfrna_normal_cols)

# Detection rate in cancer cfRNA (should be HIGH for DCB genes)
cfrna_cancer_detection_rate <- rowSums(cfrna_detection_matrix[, cfrna_cancer_cols, drop = FALSE], na.rm = TRUE) /
  length(cfrna_cancer_cols)

# Mean TPM in cfRNA samples
cfrna_cancer_mean_tpm <- rowMeans(tpm_data[, cfrna_cancer_cols, drop = FALSE], na.rm = TRUE)
cfrna_normal_mean_tpm <- rowMeans(tpm_data[, cfrna_normal_cols, drop = FALSE], na.rm = TRUE)

# Combine all metrics
full_dcb_analysis <- discovery_genes %>%
  mutate(
    cfrna_cancer_detection_rate = cfrna_cancer_detection_rate,
    cfrna_normal_detection_rate = cfrna_normal_detection_rate,
    cfrna_cancer_mean_tpm = cfrna_cancer_mean_tpm,
    cfrna_normal_mean_tpm = cfrna_normal_mean_tpm,
    # DCB criteria:
    # 1. Upregulated in tumor tissue
    # 2. Low detection rate in normal cfRNA (dark region)
    # 3. High detection rate in cancer cfRNA
    is_dcb = is_tissue_upregulated &
            (cfrna_normal_detection_rate < normal_cfrna_detection_rate) &
            (cfrna_cancer_detection_rate >= cancer_cfrna_detection_rate)
  ) %>%
  arrange(desc(is_dcb), desc(tissue_log2fc), desc(cfrna_cancer_detection_rate))

# Count DCB genes
dcb_count <- sum(full_dcb_analysis$is_dcb, na.rm = TRUE)
dark_region_count <- sum(full_dcb_analysis$cfrna_normal_detection_rate < normal_cfrna_detection_rate, na.rm = TRUE)

cat(sprintf("\n=== DCB Summary ===\n"))
cat(sprintf("Genes in dark regions (normal cfRNA DR < %.0f%%): %d\n",
            normal_cfrna_detection_rate * 100, dark_region_count))
cat(sprintf("Final DCB genes: %d\n", dcb_count))

# Create separate discovery and validation output tables
discovery_output <- full_dcb_analysis %>%
  select(gene, tissue_tumor_mean_tpm, tissue_normal_mean_tpm, tissue_log2fc,
         is_tissue_upregulated, is_dcb) %>%
  arrange(desc(is_tissue_upregulated), desc(tissue_log2fc))

validation_output <- full_dcb_analysis %>%
  select(gene, cfrna_cancer_detection_rate, cfrna_normal_detection_rate,
         cfrna_cancer_mean_tpm, cfrna_normal_mean_tpm,
         is_dcb) %>%
  arrange(desc(is_dcb), desc(cfrna_cancer_detection_rate))

# Write summary files
discovery_summary_text <- c(
  "=== Dark Channel Biomarker (DCB) - Discovery Summary ===",
  "",
  sprintf("Analysis Date: %s", Sys.time()),
  sprintf("Project: %s", project),
  "",
  "Study Design:",
  "  Discovery cohort: Tissue (tumor vs normal)",
  "  Validation cohort: cfRNA (cancer vs normal)",
  "",
  "Parameters:",
  sprintf("  Tissue TPM threshold: %.2f", tissue_tpm_threshold),
  sprintf("  Tissue log2FC threshold: %.2f", tissue_fc_threshold),
  sprintf("  cfRNA TPM detection threshold: %.2f", cfrna_tpm_threshold),
  sprintf("  Normal cfRNA max detection rate: %.0f%%", normal_cfrna_detection_rate * 100),
  sprintf("  Cancer cfRNA min detection rate: %.0f%%", cancer_cfrna_detection_rate * 100),
  "",
  "Discovery Cohort (Tissue):",
  sprintf("  Tumor samples: %d", length(tissue_tumor_cols)),
  sprintf("  Normal samples: %d", length(tissue_normal_cols)),
  sprintf("  Tissue-upregulated genes: %d", tissue_upregulated_count),
  "",
  "Validation Cohort (cfRNA):",
  sprintf("  Cancer samples: %d", length(cfrna_cancer_cols)),
  sprintf("  Normal samples: %d", length(cfrna_normal_cols)),
  sprintf("  Dark region genes: %d", dark_region_count),
  "",
  "Final DCB Genes:",
  sprintf("  Total DCB genes: %d", dcb_count)
)

validation_summary_text <- c(
  "=== Dark Channel Biomarker (DCB) - Validation Summary ===",
  "",
  sprintf("Analysis Date: %s", Sys.time()),
  sprintf("Project: %s", project),
  "",
  "Study Design:",
  "  Discovery cohort: Tissue (tumor vs normal)",
  "  Validation cohort: cfRNA (cancer vs normal)",
  "",
  "Parameters:",
  sprintf("  Tissue TPM threshold: %.2f", tissue_tpm_threshold),
  sprintf("  Tissue log2FC threshold: %.2f", tissue_fc_threshold),
  sprintf("  cfRNA TPM detection threshold: %.2f", cfrna_tpm_threshold),
  sprintf("  Normal cfRNA max detection rate: %.0f%%", normal_cfrna_detection_rate * 100),
  sprintf("  Cancer cfRNA min detection rate: %.0f%%", cancer_cfrna_detection_rate * 100),
  "",
  "Validation Cohort (cfRNA):",
  sprintf("  Cancer samples: %d", length(cfrna_cancer_cols)),
  sprintf("  Normal samples: %d", length(cfrna_normal_cols)),
  sprintf("  Dark region genes: %d", dark_region_count),
  "",
  "Final DCB Genes:",
  sprintf("  Total DCB genes: %d", dcb_count),
  "",
  "Top 10 DCB Genes (by cancer cfRNA detection rate):"
)

# Add top DCB genes to validation summary
if (dcb_count > 0) {
  top_dcb <- full_dcb_analysis %>%
    filter(is_dcb) %>%
    arrange(desc(cfrna_cancer_detection_rate)) %>%
    head(10) %>%
    mutate(
      rank = row_number(),
      info = sprintf("%d. %s (cfRNA DR: %.1f%%, TPM: %.2f, tissue log2FC: %.2f)",
                     rank, gene, cfrna_cancer_detection_rate * 100,
                     cfrna_cancer_mean_tpm, tissue_log2fc)
    ) %>%
    pull(info)

  validation_summary_text <- c(validation_summary_text, "")
  validation_summary_text <- c(validation_summary_text, top_dcb)
}

writeLines(discovery_summary_text, discovery_summary)
writeLines(validation_summary_text, validation_summary)

# Write DCB results
write_tsv(discovery_output, discovery_dcb_tsv)
write_tsv(validation_output, validation_dcb_tsv)

# Save RDS objects
saveRDS(discovery_output, discovery_dcb_rds)
saveRDS(validation_output, validation_dcb_rds)

# Visualization function for discovery (tissue)
plot_discovery_dcb <- function(dcb_df, title, output_file) {
  dcb_df <- dcb_df %>%
    mutate(
      category = case_when(
        is_dcb ~ "DCB",
        is_tissue_upregulated ~ "Tissue-upregulated (non-DCB)",
        TRUE ~ "Other"
      )
    )

  p <- ggplot(dcb_df, aes(x = tissue_normal_mean_tpm, y = tissue_tumor_mean_tpm)) +
    geom_point(aes(color = category, size = abs(tissue_log2fc)), alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = tissue_tpm_threshold, linetype = "dashed", color = "red") +
    scale_color_manual(values = c(
      "DCB" = "#D32F2F",
      "Tissue-upregulated (non-DCB)" = "#FFA726",
      "Other" = "#BDBDBD"
    )) +
    scale_size_continuous(range = c(0.5, 3), name = "|log2FC|") +
    labs(
      title = title,
      x = sprintf("Normal Tissue Mean TPM"),
      y = sprintf("Tumor Tissue Mean TPM (threshold >= %.2f)", tissue_tpm_threshold),
      color = "Category"
    ) +
    coord_fixed(xlim = c(0, max(dcb_df$tissue_tumor_mean_tpm, na.rm = TRUE) * 1.1),
                ylim = c(0, max(dcb_df$tissue_tumor_mean_tpm, na.rm = TRUE) * 1.1)) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 8, dpi = 300)

  # Create top genes plot
  if (sum(dcb_df$is_dcb, na.rm = TRUE) > 0) {
    top_dcb <- dcb_df %>%
      filter(is_dcb) %>%
      arrange(desc(tissue_log2fc)) %>%
      head(20)

    if (nrow(top_dcb) > 0) {
      p2 <- ggplot(top_dcb, aes(x = reorder(gene, tissue_log2fc), y = tissue_log2fc)) +
        geom_col(aes(fill = tissue_tumor_mean_tpm), show.legend = FALSE) +
        coord_flip() +
        scale_fill_gradient(low = "#FFF3E0", high = "#FF6F00") +
        geom_hline(yintercept = tissue_fc_threshold, linetype = "dashed", color = "red") +
        labs(
          title = sprintf("Top %d Tissue-Upregulated Genes by log2FC", min(20, sum(dcb_df$is_dcb, na.rm = TRUE))),
          x = "Gene",
          y = sprintf("Tissue log2FC (threshold >= %.2f)", tissue_fc_threshold)
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold")
        )

      ggsave(gsub("\\.png$", "_top_genes.png", output_file), p2, width = 10, height = 8, dpi = 300)
    }
  }
}

# Visualization function for validation (cfRNA)
plot_validation_dcb <- function(dcb_df, title, output_file) {
  dcb_df <- dcb_df %>%
    mutate(
      category = case_when(
        is_dcb ~ "DCB",
        cfrna_normal_detection_rate < normal_cfrna_detection_rate ~ "Dark region (non-DCB)",
        cfrna_cancer_detection_rate >= cancer_cfrna_detection_rate ~ "Detected in cancer cfRNA",
        TRUE ~ "Other"
      )
    )

  p <- ggplot(dcb_df, aes(x = cfrna_normal_detection_rate, y = cfrna_cancer_detection_rate)) +
    geom_point(aes(color = category, size = cfrna_cancer_mean_tpm), alpha = 0.6) +
    geom_hline(yintercept = cancer_cfrna_detection_rate, linetype = "dashed", color = "red") +
    geom_vline(xintercept = normal_cfrna_detection_rate, linetype = "dashed", color = "red") +
    scale_color_manual(values = c(
      "DCB" = "#D32F2F",
      "Dark region (non-DCB)" = "#FFA726",
      "Detected in cancer cfRNA" = "#42A5F5",
      "Other" = "#BDBDBD"
    )) +
    scale_size_continuous(range = c(0.5, 3), name = "Cancer cfRNA Mean TPM") +
    labs(
      title = title,
      x = sprintf("Normal cfRNA Detection Rate (threshold < %.0f%%)", normal_cfrna_detection_rate * 100),
      y = sprintf("Cancer cfRNA Detection Rate (threshold >= %.0f%%)", cancer_cfrna_detection_rate * 100),
      color = "Category"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 10, height = 8, dpi = 300)

  # Create top DCB genes plot by cfRNA detection
  if (sum(dcb_df$is_dcb, na.rm = TRUE) > 0) {
    top_dcb <- full_dcb_analysis %>%
      filter(is_dcb) %>%
      arrange(desc(cfrna_cancer_detection_rate)) %>%
      head(20)

    if (nrow(top_dcb) > 0) {
      p2 <- ggplot(top_dcb, aes(x = reorder(gene, cfrna_cancer_detection_rate), y = cfrna_cancer_detection_rate)) +
        geom_col(aes(fill = cfrna_cancer_mean_tpm), show.legend = FALSE) +
        coord_flip() +
        scale_fill_gradient(low = "#FFF3E0", high = "#FF6F00") +
        labs(
          title = sprintf("Top %d DCB Genes by Cancer cfRNA Detection", min(20, sum(dcb_df$is_dcb, na.rm = TRUE))),
          x = "Gene",
          y = "Cancer cfRNA Detection Rate"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold")
        )

      ggsave(gsub("\\.png$", "_top_genes.png", output_file), p2, width = 10, height = 8, dpi = 300)
    }
  }
}

# Generate plots
plot_discovery_dcb(
  discovery_output,
  sprintf("Discovery - Tissue Expression (%s)", project),
  discovery_plot
)

plot_validation_dcb(
  validation_output,
  sprintf("Validation - cfRNA Detection (%s)", project),
  validation_plot
)

cat("\n=== DCB Analysis Complete ===\n")
cat(sprintf("Output files written to:\n"))
cat(sprintf("  Discovery: %s\n", discovery_dcb_tsv))
cat(sprintf("  Validation: %s\n", validation_dcb_tsv))
