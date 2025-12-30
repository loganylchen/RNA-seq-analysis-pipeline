#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(purrr)
  library(broom)
  library(stats)
})

# Get parameters from Snakemake
merged_file <- snakemake@input[["merged"]]

samples_file <- snakemake@params[["samples"]]
project <- snakemake@params[["project"]]
case_condition <- snakemake@params[["case_condition"]]
control_condition <- snakemake@params[["control_condition"]]
discovery_sample_type <- snakemake@params[["discovery_sample_type"]]

# Output files
results_output <- snakemake@output[["results"]]
significant_output <- snakemake@output[["significant"]]
summary_output <- snakemake@output[["summary"]]

cat("=== ModTect Statistical Analysis ===\n")
cat("Project:", project, "\n")
cat("Case condition:", case_condition, "\n")
cat("Control condition:", control_condition, "\n")
cat("Discovery sample type:", discovery_sample_type, "\n")

# Load merged modTect data
cat("\n--- Loading merged modTect data ---\n")
modtect_data <- read_tsv(merged_file, show_col_types = FALSE)

# Load sample information
cat("--- Loading sample information ---\n")
sample_info <- read_tsv(samples_file, show_col_types = FALSE) %>%
  filter(project == !!project, sample_type == !!discovery_sample_type) %>%
  column_to_rownames(var = "sample_name")

cat("Discovery samples:", nrow(sample_info), "\n")
cat("Case samples:", sum(sample_info$condition == case_condition), "\n")
cat("Control samples:", sum(sample_info$condition == control_condition), "\n")

# Get sample columns (all columns except chrom, position, reference_nt)
sample_cols <- colnames(modtect_data)[!colnames(modtect_data) %in% c("chrom", "position", "reference_nt")]

# Filter to only discovery samples
discovery_samples <- intersect(sample_cols, rownames(sample_info))
cat("Samples in merged data:", length(sample_cols), "\n")
cat("Discovery samples found:", length(discovery_samples), "\n")

if (length(discovery_samples) == 0) {
  stop("No discovery samples found in the merged modTect data!")
}

# Get case and control sample names
case_samples <- sample_info %>%
  filter(condition == case_condition) %>%
  rownames()
control_samples <- sample_info %>%
  filter(condition == control_condition) %>%
  rownames()

cat("Case samples in data:", length(intersect(case_samples, discovery_samples)), "\n")
cat("Control samples in data:", length(intersect(control_samples, discovery_samples)), "\n")

# Prepare data for statistical analysis
cat("\n--- Preparing data for statistical analysis ---\n")
analysis_data <- modtect_data %>%
  select(all_of(c("chrom", "position", "reference_nt", discovery_samples))) %>%
  mutate(
    site_id = paste(chrom, position, reference_nt, sep = ":")
  ) %>%
  select(site_id, all_of(discovery_samples))

cat("Total sites:", nrow(analysis_data), "\n")

# Add condition information
case_cols <- intersect(discovery_samples, case_samples)
control_cols <- intersect(discovery_samples, control_samples)

# Function to perform statistical test for each site
perform_site_test <- function(site_data) {
  # Extract values for case and control
  case_values <- site_data[case_cols]
  control_values <- site_data[control_cols]

  # Convert to numeric
  case_values <- as.numeric(case_values)
  control_values <- as.numeric(control_values)

  # Handle missing values (NA)
  case_values <- case_values[!is.na(case_values)]
  control_values <- control_values[!is.na(control_values)]

  # Calculate summary statistics
  n_case <- length(case_values)
  n_control <- length(control_values)

  # Check if we have enough data
  if (n_case < 2 || n_control < 2) {
    return(list(
      p_value = NA,
      statistic = NA,
      method = "insufficient_data",
      n_case = n_case,
      n_control = n_control
    ))
  }

  # Calculate means
  mean_case <- mean(case_values, na.rm = TRUE)
  mean_control <- mean(control_values, na.rm = TRUE)

  # Calculate fold change (log2)
  # Add small pseudocount to avoid division by zero
  pseudocount <- 1e-6
  log2fc <- log2((mean_case + pseudocount) / (mean_control + pseudocount))

  # Perform Wilcoxon rank-sum test (non-parametric)
  tryCatch({
    test_result <- wilcox.test(case_values, control_values, exact = FALSE, correct = TRUE)

    return(list(
      p_value = test_result$p.value,
      statistic = test_result$statistic,
      method = "wilcoxon_rank_sum",
      n_case = n_case,
      n_control = n_control,
      mean_case = mean_case,
      mean_control = mean_control,
      log2fc = log2fc,
      sd_case = sd(case_values),
      sd_control = sd(control_values)
    ))
  }, error = function(e) {
    return(list(
      p_value = NA,
      statistic = NA,
      method = "error",
      n_case = n_case,
      n_control = n_control,
      mean_case = mean_case,
      mean_control = mean_control,
      log2fc = log2fc
    ))
  })
}

# Apply statistical test to each site
cat("\n--- Performing statistical tests ---\n")
cat("Testing", nrow(analysis_data), "sites...\n")

# Initialize results data frame
results_list <- list()

for (i in 1:nrow(analysis_data)) {
  if (i %% 1000 == 0) {
    cat("Processed", i, "sites...\n")
  }

  site_data <- analysis_data[i, ]
  test_result <- perform_site_test(site_data)

  results_list[[i]] <- data.frame(
    site_id = site_data$site_id,
    chrom = modtect_data$chrom[i],
    position = modtect_data$position[i],
    reference_nt = modtect_data$reference_nt[i],
    p_value = test_result$p_value,
    statistic = test_result$statistic,
    test_method = test_result$method,
    n_case = test_result$n_case,
    n_control = test_result$n_control,
    mean_case = ifelse(is.null(test_result$mean_case), NA, test_result$mean_case),
    mean_control = ifelse(is.null(test_result$mean_control), NA, test_result$mean_control),
    log2fc = ifelse(is.null(test_result$log2fc), NA, test_result$log2fc),
    sd_case = ifelse(is.null(test_result$sd_case), NA, test_result$sd_case),
    sd_control = ifelse(is.null(test_result$sd_control), NA, test_result$sd_control)
  )
}

results_df <- bind_rows(results_list)

# Adjust p-values for multiple testing (FDR using Benjamini-Hochberg)
cat("\n--- Adjusting p-values (FDR) ---\n")
results_df <- results_df %>%
  mutate(
    padj = p.adjust(p_value, method = "BH"),
    significant = padj < 0.05
  )

cat("Sites tested:", nrow(results_df), "\n")
cat("Sites with valid p-values:", sum(!is.na(results_df$p_value)), "\n")
cat("Significant sites (FDR < 0.05):", sum(results_df$significant, na.rm = TRUE), "\n")

# Write full results
cat("\n--- Writing results ---\n")
write_tsv(results_df, results_output)
cat("Full results written to:", results_output, "\n")

# Write significant sites only
significant_df <- results_df %>%
  filter(significant == TRUE) %>%
  arrange(padj)

write_tsv(significant_df, significant_output)
cat("Significant modifications written to:", significant_output, "\n")
cat("Number of significant modifications:", nrow(significant_df), "\n")

# Write summary statistics
summary_text <- c(
  "=== ModTect Statistical Analysis Summary ===",
  "",
  paste("Project:", project),
  paste("Case condition:", case_condition),
  paste("Control condition:", control_condition),
  paste("Discovery sample type:", discovery_sample_type),
  "",
  "--- Sample Counts ---",
  paste("Total discovery samples:", length(discovery_samples)),
  paste("Case samples:", length(intersect(case_samples, discovery_samples))),
  paste("Control samples:", length(intersect(control_samples, discovery_samples))),
  "",
  "--- Site Statistics ---",
  paste("Total sites tested:", nrow(results_df)),
  paste("Sites with valid p-values:", sum(!is.na(results_df$p_value))),
  paste("Sites with insufficient data:", sum(results_df$test_method == "insufficient_data", na.rm = TRUE)),
  "",
  "--- Significance Testing ---",
  paste("Significant sites (FDR < 0.05):", sum(results_df$significant, na.rm = TRUE)),
  paste("Percentage significant:", round(sum(results_df$significant, na.rm = TRUE) / sum(!is.na(results_df$p_value)) * 100, 2), "%"),
  "",
  "--- Top Significant Sites ---",
  if (nrow(significant_df) > 0) {
    head(significant_df) %>%
      transmute(Site = site_id,
                Chrom = chrom,
                Pos = position,
                `log2FC` = round(log2fc, 3),
                `P-value` = formatC(p_value, format = "e", digits = 2),
                FDR = formatC(padj, format = "e", digits = 2)) %>%
      as.data.frame() %>%
      capture.output(type = "output", cat(., sep = "\n"))
  } else {
    "No significant sites found"
  }
)

writeLines(summary_text, summary_output)
cat("Summary written to:", summary_output, "\n")

cat("\n=== ModTect Statistical Analysis Complete ===\n")

sink()
