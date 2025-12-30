#!/usr/bin/env Rscript

# DEG-based LASSO Classifier
# Builds a LASSO logistic regression model using DEGs from discovery dataset
# Tests performance on validation dataset
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(glmnet)
  library(pROC)
  library(caret)
  library(ggplot2)
  library(ggsci)
})

# Get parameters from Snakemake
discovery_deg_rds <- snakemake@input[["discovery_deg_rds"]]
validation_deg_rds <- snakemake@input[["validation_deg_rds"]]
discovery_tpm_file <- snakemake@input[["discovery_tpm"]]
validation_tpm_file <- snakemake@input[["validation_tpm"]]
samples_file <- snakemake@params[["samples"]]
project <- snakemake@params[["project"]]
case_condition <- snakemake@params[["case_condition"]]
control_condition <- snakemake@params[["control_condition"]]
discovery_sample_type <- snakemake@params[["discovery_sample_type"]]
log2fc_threshold <- snakemake@params[["log2fc_threshold"]]
padj_threshold <- snakemake@params[["padj_threshold"]]

# Output files
signature_output <- snakemake@output[["signature"]]
coefficients_output <- snakemake@output[["coefficients"]]
discovery_predictions <- snakemake@output[["discovery_predictions"]]
validation_predictions <- snakemake@output[["validation_predictions"]]
roc_plot <- snakemake@output[["roc_plot"]]
summary_output <- snakemake@output[["summary"]]

cat("=== DEG-based LASSO Classifier ===\n")
cat("Project:", project, "\n")
cat("Case condition:", case_condition, "\n")
cat("Control condition:", control_condition, "\n")
cat("Discovery sample type:", discovery_sample_type, "\n")
cat("Log2FC threshold:", log2fc_threshold, "\n")
cat("Padj threshold:", padj_threshold, "\n")

# Load data
cat("\n--- Loading data ---\n")
discovery_deg <- readRDS(discovery_deg_rds)
validation_deg <- readRDS(validation_deg_rds)

# Load TPM matrices
discovery_tpm <- read_tsv(discovery_tpm_file, show_col_types = FALSE)
validation_tpm <- read_tsv(validation_tpm_file, show_col_types = FALSE)

# Get gene names column (first column)
gene_col_discovery <- colnames(discovery_tpm)[1]
gene_col_validation <- colnames(validation_tpm)[1]

# Set row names as gene IDs
discovery_tpm <- discovery_tpm %>%
  column_to_rownames(var = gene_col_discovery)
validation_tpm <- validation_tpm %>%
  column_to_rownames(var = gene_col_validation)

sample_info <- read_tsv(samples_file, show_col_types = FALSE) %>%
  filter(project == !!project) %>%
  column_to_rownames(var = "sample_name")

cat("Discovery DEGs:", nrow(discovery_deg), "\n")
cat("Validation DEGs:", nrow(validation_deg), "\n")
cat("Discovery TPM:", nrow(discovery_tpm), "genes x", ncol(discovery_tpm), "samples\n")
cat("Validation TPM:", nrow(validation_tpm), "genes x", ncol(validation_tpm), "samples\n")

# Get significant DEGs from discovery
cat("\n--- Selecting significant DEGs from discovery ---\n")
discovery_sig_genes <- discovery_deg %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  filter(padj < padj_threshold, abs(log2FoldChange) >= log2fc_threshold) %>%
  pull(gene_id)

cat("Significant DEGs in discovery:", length(discovery_sig_genes), "\n")

# Get intersection with available TPM genes
discovery_available <- intersect(discovery_sig_genes, rownames(discovery_tpm))
validation_available <- intersect(discovery_available, rownames(validation_tpm))

cat("DEGs available in discovery TPM:", length(discovery_available), "\n")
cat("DEGs available in validation TPM:", length(validation_available), "\n")

if (length(validation_available) < 5) {
  stop("Too few DEGs available for LASSO classification!")
}

# Use genes available in both datasets
feature_genes <- validation_available
cat("Using", length(feature_genes), "genes for classification\n")

# Prepare discovery data
cat("\n--- Preparing discovery dataset ---\n")
discovery_samples <- colnames(discovery_tpm)
discovery_sample_info <- sample_info[discovery_samples, , drop = FALSE]

# Filter to discovery sample type only
discovery_samples_filtered <- discovery_sample_info %>%
  filter(sample_type == !!discovery_sample_type) %>%
  rownames()

cat("Discovery samples (filtered):", length(discovery_samples_filtered), "\n")

# Create design matrix and response
# Use log2(TPM + 1) for better numerical stability
x_discovery <- t(log2(discovery_tpm[feature_genes, discovery_samples_filtered, drop = FALSE] + 1))
y_discovery <- ifelse(discovery_sample_info[discovery_samples_filtered, "condition"] == case_condition, 1, 0)

cat("Case samples:", sum(y_discovery), "\n")
cat("Control samples:", sum(1 - y_discovery), "\n")

# Prepare validation data
cat("\n--- Preparing validation dataset ---\n")
validation_samples <- colnames(validation_tpm)
validation_sample_info <- sample_info[validation_samples, , drop = FALSE]

# Filter to NON-discovery sample type (validation cohort)
validation_samples_filtered <- validation_sample_info %>%
  filter(sample_type != !!discovery_sample_type) %>%
  rownames()

cat("Validation samples (filtered):", length(validation_samples_filtered), "\n")

x_validation <- t(log2(validation_tpm[feature_genes, validation_samples_filtered, drop = FALSE] + 1))
y_validation <- ifelse(validation_sample_info[validation_samples_filtered, "condition"] == case_condition, 1, 0)

cat("Case samples:", sum(y_validation), "\n")
cat("Control samples:", sum(1 - y_validation), "\n")

# Fit LASSO logistic regression with cross-validation
cat("\n--- Fitting LASSO model with cross-validation ---\n")
set.seed(42)

# Use cv.glmnet for cross-validation
cv_fit <- cv.glmnet(
  x = x_discovery,
  y = y_discovery,
  family = "binomial",
  alpha = 1,  # LASSO
  nfolds = min(5, sum(y_discovery), sum(1 - y_discovery)),  # Ensure both classes in each fold
  type.measure = "deviance"
)

cat("Cross-validation complete\n")
cat("Optimal lambda:", cv_fit$lambda.min, "\n")
cat("Lambda within 1 SE:", cv_fit$lambda.1se, "\n")

# Get coefficients at optimal lambda
coefs <- coef(cv_fit, s = "lambda.min")
selected_genes <- rownames(coefs)[coefs[,1] != 0][-1]  # Exclude intercept

cat("Number of selected genes:", length(selected_genes), "\n")
if (length(selected_genes) > 0) {
  cat("Selected genes:", paste(selected_genes, collapse = ", "), "\n")
}

# Make predictions on discovery data
cat("\n--- Evaluating on discovery dataset ---\n")
discovery_pred_prob <- predict(cv_fit, newx = x_discovery, s = "lambda.min", type = "response")
discovery_pred_class <- ifelse(discovery_pred_prob > 0.5, 1, 0)

# Discovery performance
discovery_cm <- table(Predicted = discovery_pred_class, Actual = y_discovery)
discovery_accuracy <- sum(diag(discovery_cm)) / sum(discovery_cm)
discovery_sensitivity <- discovery_cm["1", "1"] / sum(y_discovery == 1)
discovery_specificity <- discovery_cm["0", "0"] / sum(y_discovery == 0)

cat("Discovery Accuracy:", round(discovery_accuracy, 3), "\n")
cat("Discovery Sensitivity:", round(discovery_sensitivity, 3), "\n")
cat("Discovery Specificity:", round(discovery_specificity, 3), "\n")

# Discovery AUC
discovery_roc <- roc(y_discovery, as.numeric(discovery_pred_prob))
discovery_auc <- auc(discovery_roc)
cat("Discovery AUC:", round(discovery_auc, 3), "\n")

# Make predictions on validation data
cat("\n--- Evaluating on validation dataset ---\n")
validation_pred_prob <- predict(cv_fit, newx = x_validation, s = "lambda.min", type = "response")
validation_pred_class <- ifelse(validation_pred_prob > 0.5, 1, 0)

# Validation performance
validation_cm <- table(Predicted = validation_pred_class, Actual = y_validation)
validation_accuracy <- sum(diag(validation_cm)) / sum(validation_cm)
validation_sensitivity <- ifelse("1" %in% rownames(validation_cm) && "1" %in% colnames(validation_cm),
                                 validation_cm["1", "1"] / sum(y_validation == 1), NA)
validation_specificity <- ifelse("0" %in% rownames(validation_cm) && "0" %in% colnames(validation_cm),
                                 validation_cm["0", "0"] / sum(y_validation == 0), NA)

cat("Validation Accuracy:", round(validation_accuracy, 3), "\n")
if (!is.na(validation_sensitivity)) {
  cat("Validation Sensitivity:", round(validation_sensitivity, 3), "\n")
}
if (!is.na(validation_specificity)) {
  cat("Validation Specificity:", round(validation_specificity, 3), "\n")
}

# Validation AUC
validation_roc <- roc(y_validation, as.numeric(validation_pred_prob))
validation_auc <- auc(validation_roc)
cat("Validation AUC:", round(validation_auc, 3), "\n")

# Write signature gene list
cat("\n--- Writing signature gene list ---\n")
signature_df <- data.frame(
  gene_id = selected_genes,
  stringsAsFactors = FALSE
)
write_tsv(signature_df, signature_output)
cat("Signature written to:", signature_output, "\n")

# Write coefficients
cat("\n--- Writing coefficients ---\n")
coef_df <- data.frame(
  gene_id = rownames(coefs),
  coefficient = as.numeric(coefs),
  stringsAsFactors = FALSE
)
coef_df <- coef_df[coef_df$gene_id != "(Intercept)", ]  # Remove intercept
coef_df <- coef_df[order(abs(coef_df$coefficient), decreasing = TRUE), ]
write_tsv(coef_df, coefficients_output)
cat("Coefficients written to:", coefficients_output, "\n")

# Write prediction results
cat("\n--- Writing prediction results ---\n")
discovery_pred_df <- data.frame(
  sample = discovery_samples_filtered,
  condition = discovery_sample_info[discovery_samples_filtered, "condition"],
  actual = y_discovery,
  predicted_prob = as.numeric(discovery_pred_prob),
  predicted_class = discovery_pred_class
)
write_tsv(discovery_pred_df, discovery_predictions)

validation_pred_df <- data.frame(
  sample = validation_samples_filtered,
  condition = validation_sample_info[validation_samples_filtered, "condition"],
  actual = y_validation,
  predicted_prob = as.numeric(validation_pred_prob),
  predicted_class = validation_pred_class
)
write_tsv(validation_pred_df, validation_predictions)
cat("Predictions written\n")

# Create ROC plot
cat("\n--- Creating ROC plot ---\n")
roc_data <- data.frame(
  TPR = c(discovery_roc$sensitivities, validation_roc$sensitivities),
  FPR = c(1 - discovery_roc$specificities, 1 - validation_roc$specificities),
  Dataset = c(rep("Discovery", length(discovery_roc$sensitivities)),
              rep("Validation", length(validation_roc$sensitivities)))
)

p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Dataset)) +
  geom_line(size = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  annotate("text", x = 0.75, y = 0.25,
           label = sprintf("Discovery AUC = %.3f\nValidation AUC = %.3f",
                          discovery_auc, validation_auc),
           hjust = 0, size = 4) +
  labs(x = "False Positive Rate", y = "True Positive Rate",
       title = sprintf("LASSO Classifier ROC Curve (%s)", project),
       color = "Dataset") +
  scale_color_manual(values = c("Discovery" = "#377EB8", "Validation" = "#4DAF4A")) +
  coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "bottom"
  )

ggsave(roc_plot, p_roc, width = 8, height = 6, dpi = 300)
cat("ROC plot saved to:", roc_plot, "\n")

# Write summary
cat("\n--- Writing summary ---\n")
summary_text <- c(
  "=== LASSO Classifier Summary ===",
  "",
  paste("Project:", project),
  paste("Case condition:", case_condition),
  paste("Control condition:", control_condition),
  paste("Discovery sample type:", discovery_sample_type),
  "",
  paste("Analysis Date:", Sys.time()),
  "",
  "--- Feature Selection ---",
  paste("Total DEGs tested:", length(feature_genes)),
  paste("Signature genes selected:", length(selected_genes)),
  if (length(selected_genes) > 0) {
    c(paste("Selected genes:", paste(selected_genes, collapse = ", ")),
      "")
  } else {
    ""
  },
  "--- Discovery Performance ---",
  paste("Samples:", length(discovery_samples_filtered)),
  paste("  Case:", sum(y_discovery)),
  paste("  Control:", sum(1 - y_discovery)),
  sprintf("Accuracy: %.3f", discovery_accuracy),
  sprintf("Sensitivity: %.3f", discovery_sensitivity),
  sprintf("Specificity: %.3f", discovery_specificity),
  sprintf("AUC: %.3f", discovery_auc),
  "",
  "--- Validation Performance ---",
  paste("Samples:", length(validation_samples_filtered)),
  paste("  Case:", sum(y_validation)),
  paste("  Control:", sum(1 - y_validation)),
  sprintf("Accuracy: %.3f", validation_accuracy),
  if (!is.na(validation_sensitivity)) {
    sprintf("Sensitivity: %.3f", validation_sensitivity)
  } else {
    "Sensitivity: NA (no positive predictions)"
  },
  if (!is.na(validation_specificity)) {
    sprintf("Specificity: %.3f", validation_specificity)
  } else {
    "Specificity: NA (no negative predictions)"
  },
  sprintf("AUC: %.3f", validation_auc)
)

writeLines(summary_text, summary_output)
cat("Summary written to:", summary_output, "\n")

cat("\n=== LASSO Classifier Analysis Complete ===\n")

sink()
