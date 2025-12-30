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

# The merged data has columns: chrom, position, reference_nt, then for each sample:
# - SampleName_ModTect_score
# - SampleName_variant_proportion
# We need to separate these

cat("\n--- Parsing merged ModTect data ---\n")
cat("Total columns:", ncol(modtect_data), "\n")
cat("Column names:", head(colnames(modtect_data), 10), "\n")

# Get all column names except the first 3 (chrom, position, reference_nt)
all_cols <- colnames(modtect_data)[-(1:3)]

# Extract sample names from column names
# Sample columns have prefix ModTect_score_ or variant_proportion_
score_cols <- all_cols[grepl("^ModTect_score_", all_cols)]
prop_cols <- all_cols[grepl("^variant_proportion_", all_cols)]

# Extract sample names by removing the prefix
sample_names_score <- gsub("^ModTect_score_", "", score_cols)
sample_names_prop <- gsub("^variant_proportion_", "", prop_cols)

# Get unique sample names
all_samples <- unique(c(sample_names_score, sample_names_prop))
cat("Unique samples found in data:", length(all_samples), "\n")

# Filter to only discovery samples
discovery_samples <- intersect(all_samples, rownames(sample_info))
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
# We need to work with variant_proportion for statistics
# And ModTect_score for counting > 10

cat("\n--- Preparing data for statistical analysis ---\n")

# Create site_id
modtect_data <- modtect_data %>%
  mutate(site_id = paste(chrom, position, reference_nt, sep = ":"))

# Get case and control sample columns
case_cols <- intersect(discovery_samples, case_samples)
control_cols <- intersect(discovery_samples, control_samples)

# Build the proportion column names for case and control
case_prop_cols <- paste0("variant_proportion_", case_cols)
control_prop_cols <- paste0("variant_proportion_", control_cols)

# Build the score column names for case and control
case_score_cols <- paste0("ModTect_score_", case_cols)
control_score_cols <- paste0("ModTect_score_", control_cols)

# Check which columns exist in data
case_prop_cols <- intersect(case_prop_cols, colnames(modtect_data))
control_prop_cols <- intersect(control_prop_cols, colnames(modtect_data))
case_score_cols <- intersect(case_score_cols, colnames(modtect_data))
control_score_cols <- intersect(control_score_cols, colnames(modtect_data))

cat("Case proportion columns found:", length(case_prop_cols), "\n")
cat("Control proportion columns found:", length(control_prop_cols), "\n")
cat("Case score columns found:", length(case_score_cols), "\n")
cat("Control score columns found:", length(control_score_cols), "\n")

if (length(case_prop_cols) == 0 || length(control_prop_cols) == 0) {
  stop("Insufficient proportion columns for statistical analysis!")
}

# Extract the data we need
analysis_data <- modtect_data %>%
  select(site_id, all_of(c(case_prop_cols, control_prop_cols, case_score_cols, control_score_cols)))

cat("Total sites:", nrow(analysis_data), "\n")

# Extract sample names from column names (remove prefix)
case_sample_names <- gsub("^variant_proportion_", "", case_prop_cols)
control_sample_names <- gsub("^variant_proportion_", "", control_prop_cols)

# Function to perform statistical test for each site
perform_site_test <- function(site_data) {
  # Extract proportion values for case and control
  case_prop_values <- site_data[case_prop_cols]
  control_prop_values <- site_data[control_prop_cols]

  # Extract score values for counting > 10
  case_score_values <- site_data[case_score_cols]
  control_score_values <- site_data[control_score_cols]

  # Convert to numeric
  case_prop_values <- as.numeric(case_prop_values)
  control_prop_values <- as.numeric(control_prop_values)
  case_score_values <- as.numeric(case_score_values)
  control_score_values <- as.numeric(control_score_values)

  # Handle missing values (NA) - remove from both prop and score together
  # For statistics, use only non-NA proportion values
  case_na <- is.na(case_prop_values)
  control_na <- is.na(control_prop_values)

  case_prop_values <- case_prop_values[!case_na]
  control_prop_values <- control_prop_values[!control_na]
  case_score_values <- case_score_values[!case_na]
  control_score_values <- control_score_values[!control_na]

  # Calculate summary statistics
  n_case <- length(case_prop_values)
  n_control <- length(control_prop_values)

  # Check if we have enough data
  if (n_case < 2 || n_control < 2) {
    return(list(
      p_value = NA,
      statistic = NA,
      method = "insufficient_data",
      n_case = n_case,
      n_control = n_control,
      mean_prop_case = NA,
      mean_prop_control = NA,
      case_score_gt10_str = NA,
      control_score_gt10_str = NA,
      case_score_gt10_prop = NA,
      control_score_gt10_prop = NA,
      log2fc = NA
    ))
  }

  # Calculate means of variant proportion
  mean_prop_case <- mean(case_prop_values, na.rm = TRUE)
  mean_prop_control <- mean(control_prop_values, na.rm = TRUE)

  # Calculate proportion of samples with ModTect_score > 10
  case_score_gt10 <- sum(case_score_values > 10, na.rm = TRUE)
  control_score_gt10 <- sum(control_score_values > 10, na.rm = TRUE)
  case_score_gt10_str <- sprintf("%d/%d", case_score_gt10, n_case)
  control_score_gt10_str <- sprintf("%d/%d", control_score_gt10, n_control)
  case_score_gt10_prop <- case_score_gt10 / n_case
  control_score_gt10_prop <- control_score_gt10 / n_control

  # Calculate fold change (log2) using variant proportion
  # Add small pseudocount to avoid division by zero
  pseudocount <- 1e-6
  log2fc <- log2((mean_prop_case + pseudocount) / (mean_prop_control + pseudocount))

  # Perform Wilcoxon rank-sum test on variant proportion (non-parametric)
  tryCatch({
    test_result <- wilcox.test(case_prop_values, control_prop_values, exact = FALSE, correct = TRUE)

    return(list(
      p_value = test_result$p.value,
      statistic = test_result$statistic,
      method = "wilcoxon_rank_sum",
      n_case = n_case,
      n_control = n_control,
      mean_prop_case = mean_prop_case,
      mean_prop_control = mean_prop_control,
      case_score_gt10_str = case_score_gt10_str,
      control_score_gt10_str = control_score_gt10_str,
      case_score_gt10_prop = case_score_gt10_prop,
      control_score_gt10_prop = control_score_gt10_prop,
      log2fc = log2fc,
      sd_prop_case = sd(case_prop_values),
      sd_prop_control = sd(control_prop_values)
    ))
  }, error = function(e) {
    return(list(
      p_value = NA,
      statistic = NA,
      method = "error",
      n_case = n_case,
      n_control = n_control,
      mean_prop_case = mean_prop_case,
      mean_prop_control = mean_prop_control,
      case_score_gt10_str = case_score_gt10_str,
      control_score_gt10_str = control_score_gt10_str,
      case_score_gt10_prop = case_score_gt10_prop,
      control_score_gt10_prop = control_score_gt10_prop,
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
    mean_prop_case = ifelse(is.null(test_result$mean_prop_case), NA, test_result$mean_prop_case),
    mean_prop_control = ifelse(is.null(test_result$mean_prop_control), NA, test_result$mean_prop_control),
    case_score_gt10_str = ifelse(is.null(test_result$case_score_gt10_str), NA, test_result$case_score_gt10_str),
    control_score_gt10_str = ifelse(is.null(test_result$control_score_gt10_str), NA, test_result$control_score_gt10_str),
    case_score_gt10_prop = ifelse(is.null(test_result$case_score_gt10_prop), NA, test_result$case_score_gt10_prop),
    control_score_gt10_prop = ifelse(is.null(test_result$control_score_gt10_prop), NA, test_result$control_score_gt10_prop),
    log2fc = ifelse(is.null(test_result$log2fc), NA, test_result$log2fc),
    sd_prop_case = ifelse(is.null(test_result$sd_prop_case), NA, test_result$sd_prop_case),
    sd_prop_control = ifelse(is.null(test_result$sd_prop_control), NA, test_result$sd_prop_control)
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
  "--- Column Descriptions ---",
  "mean_prop_case: Mean variant proportion in case samples",
  "mean_prop_control: Mean variant proportion in control samples",
  "case_score_gt10_str: Case samples with ModTect_score > 10 (e.g., '1/11')",
  "control_score_gt10_str: Control samples with ModTect_score > 10 (e.g., '0/10')",
  "case_score_gt10_prop: Proportion of case samples with ModTect_score > 10 (float, e.g., 0.09)",
  "control_score_gt10_prop: Proportion of control samples with ModTect_score > 10 (float, e.g., 0.00)",
  "log2fc: log2 fold change of variant proportion (case/control)",
  "",
  "--- Top Significant Sites ---",
  if (nrow(significant_df) > 0) {
    head(significant_df) %>%
      transmute(Site = site_id,
                Chrom = chrom,
                Pos = position,
                `Mean Prop Case` = round(mean_prop_case, 4),
                `Mean Prop Control` = round(mean_prop_control, 4),
                `Case Score>10` = case_score_gt10_str,
                `Control Score>10` = control_score_gt10_str,
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
