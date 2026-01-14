#!/usr/bin/env Rscript
# Mime Response Prediction Analysis with Model Saving
# This script:
# 1. Trains models on discovery dataset (Dataset1)
# 2. Validates/benchmarks on all datasets (training + validation)
# 3. Saves trained models for future implementation on new datasets

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
    library(Mime)
    library(dplyr)
})

cat("==============================================================\n")
cat("Mime Response Prediction Analysis\n")
cat("==============================================================\n\n")

# Get parameters from Snakemake
combined_rds <- snakemake@input[["combined_rds"]]
genelist_file <- snakemake@input[["genelist"]]
output_dir <- snakemake@output[["output_dir"]]
project <- snakemake@params[["project"]]
tool <- snakemake@params[["tool"]]
methods <- snakemake@params[["methods"]]
seed <- snakemake@params[["seed"]]

cat("Parameters:\n")
cat("  Project:", project, "\n")
cat("  Tool:", tool, "\n")
cat("  Combined RDS:", combined_rds, "\n")
cat("  Gene list:", genelist_file, "\n")
cat("  Output directory:", output_dir, "\n")
cat("  Methods:", paste(methods, collapse=", "), "\n")
cat("  Seed:", seed, "\n\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading combined datasets...\n")
load(combined_rds)  # Creates mime_datasets
cat("  Loaded:", length(mime_datasets), "datasets\n")
for (name in names(mime_datasets)) {
    cat("    ", name, ":", nrow(mime_datasets[[name]]), "samples\n")
}

# Rename for Mime compatibility
list_train_vali_Data <- mime_datasets

# Load gene list
cat("\nLoading gene list...\n")
genelist <- read.delim(genelist_file, header=FALSE, stringsAsFactors = FALSE)$V1
cat("  Loaded", length(genelist), "genes\n")

# ============================================================================
# TRAIN RESPONSE PREDICTION MODELS
# ============================================================================

cat("\n==============================================================\n")
cat("Step 1: Training Response Prediction Models\n")
cat("==============================================================\n\n")

cat("Training on:", names(list_train_vali_Data)[1], "\n")
cat("Validating on:", paste(names(list_train_vali_Data), collapse=", "), "\n\n")

# Train models
res.ici <- ML.Dev.Pred.Category.Sig(
    train_data = list_train_vali_Data$Dataset1,
    list_train_vali_Data = list_train_vali_Data,
    candidate_genes = genelist,
    methods = methods,
    seed = seed,
    cores_for_parallel = 4
)

cat("\nTraining complete!\n")
cat("Models trained:", paste(names(res.ici), collapse=", "), "\n")

# ============================================================================
# SAVE TRAINED MODELS
# ============================================================================

cat("\n==============================================================\n")
cat("Step 2: Saving Trained Models\n")
cat("==============================================================\n\n")

model_save_path <- file.path(output_dir, "trained_models.rds")
saveRDS(res.ici, file = model_save_path)
cat("  Saved models to:", model_save_path, "\n")

# Also save model metadata
model_metadata <- list(
    project = project,
    tool = tool,
    train_dataset = names(list_train_vali_Data)[1],
    validation_datasets = names(list_train_vali_Data)[-1],
    all_datasets = names(list_train_vali_Data),
    n_genes = length(genelist),
    methods = methods,
    seed = seed,
    train_date = Sys.Date(),
    model_file = model_save_path
)

metadata_path <- file.path(output_dir, "model_metadata.rds")
saveRDS(model_metadata, file = metadata_path)
cat("  Saved metadata to:", metadata_path, "\n")

# ============================================================================
# BENCHMARK ON ALL DATASETS
# ============================================================================

cat("\n==============================================================\n")
cat("Step 3: Benchmarking on All Datasets\n")
cat("==============================================================\n\n")

# Calculate AUC for each dataset
auc_results <- list()
for (model_name in names(res.ici)) {
    cat("\nModel:", model_name, "\n")

    auc_results[[model_name]] <- list()

    for (dataset_name in names(list_train_vali_Data)) {
        # Get predictions for this dataset
        pred_obj <- res.ici[[model_name]]$prediction_list[[dataset_name]]

        if (!is.null(pred_obj)) {
            auc_results[[model_name]][[dataset_name]] <- pred_obj$auc
            cat("  ", dataset_name, "AUC:", round(pred_obj$auc, 3), "\n")
        }
    }
}

# Save AUC results
auc_path <- file.path(output_dir, "benchmark_auc_results.rds")
saveRDS(auc_results, file = auc_path)
cat("\n  Saved AUC results to:", auc_path, "\n")

# Create summary table
auc_summary <- data.frame(
    Model = character(),
    Dataset = character(),
    AUC = numeric(),
    stringsAsFactors = FALSE
)

for (model_name in names(auc_results)) {
    for (dataset_name in names(auc_results[[model_name]])) {
        auc_summary <- rbind(auc_summary, data.frame(
            Model = model_name,
            Dataset = dataset_name,
            AUC = auc_results[[model_name]][[dataset_name]],
            stringsAsFactors = FALSE
        ))
    }
}

auc_summary_file <- file.path(output_dir, "auc_summary.tsv")
write.table(auc_summary, auc_summary_file, sep="\t", row.names=FALSE, quote=FALSE)
cat("  Saved AUC summary to:", auc_summary_file, "\n")

# ============================================================================
# GENERATE VISUALIZATIONS
# ============================================================================

cat("\n==============================================================\n")
cat("Step 4: Generating Visualizations\n")
cat("==============================================================\n\n")

# 1. AUC distribution across all datasets
cat("  Creating AUC distribution plot...\n")
pdf(file.path(output_dir, "auc_distribution_all_datasets.pdf"), width = 12, height = 6)
auc_vis_category_all(
    res.ici,
    dataset = names(list_train_vali_Data),
    order = names(list_train_vali_Data)
)
dev.off()

# 2. ROC curves for each method
cat("  Creating ROC curves...\n")
pdf(file.path(output_dir, "roc_curves_all_models.pdf"), width = 15, height = 12)
plot_list <- list()
for (model_name in methods) {
    if (model_name %in% names(res.ici)) {
        plot_list[[model_name]] <- roc_vis_category(
            res.ici,
            model_name = model_name,
            dataset = names(list_train_vali_Data),
            order = names(list_train_vali_Data),
            anno_position = c(0.4, 0.25)
        )
    }
}
if (length(plot_list) > 0) {
    gridExtra::grid.arrange(grobs = plot_list, ncol = 3)
}
dev.off()

# 3. Performance comparison bar plot
cat("  Creating performance comparison plot...\n")
pdf(file.path(output_dir, "performance_comparison.pdf"), width = 10, height = 6)
ggplot2::ggplot(auc_summary, ggplot2::aes(x = Model, y = AUC, fill = Dataset)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_text(ggplot2::aes(label = round(AUC, 3)),
                       position = ggplot2::position_dodge(width = 0.9),
                       vjust = -0.25, size = 3) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
        legend.position = "right",
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    ) +
    ggplot2::labs(
        title = paste("Model Performance Comparison -", project, "-", tool),
        y = "AUC",
        fill = "Dataset"
    ) +
    ggplot2::ylim(0, 1)
dev.off()

# ============================================================================
# SAVE BEST MODEL SELECTION
# ============================================================================

cat("\n==============================================================\n")
cat("Step 5: Selecting Best Model\n")
cat("==============================================================\n\n")

# Calculate mean AUC across all datasets for each model
mean_auc <- sapply(names(auc_results), function(model_name) {
    aucs <- sapply(auc_results[[model_name]], function(x) x)
    mean(aucs, na.rm = TRUE)
})

# Sort by mean AUC
mean_auc_sorted <- sort(mean_auc, decreasing = TRUE)
best_model <- names(mean_auc_sorted)[1]

cat("Model ranking (by mean AUC across all datasets):\n")
for (i in seq_along(mean_auc_sorted)) {
    cat("  ", i, ".", names(mean_auc_sorted)[i], ":", round(mean_auc_sorted[i], 3), "\n")
}
cat("\nBest model:", best_model, "\n")

# Save best model info
best_model_info <- list(
    best_model = best_model,
    mean_auc = mean_auc_sorted[1],
    all_model_performance = mean_auc_sorted
)

best_model_path <- file.path(output_dir, "best_model_info.rds")
saveRDS(best_model_info, file = best_model_path)
cat("  Saved best model info to:", best_model_path, "\n")

# ============================================================================
# CREATE IMPLEMENTATION GUIDE
# ============================================================================

cat("\n==============================================================\n")
cat("Step 6: Creating Implementation Guide\n")
cat("==============================================================\n\n")

implementation_script <- paste0(
    "#!/usr/bin/env Rscript\n",
    "# Implementation script for trained Mime model\n",
    "# Use this script to apply the trained model to new datasets\n\n",
    "# Load the trained model\n",
    "trained_model <- readRDS('", model_save_path, "')\n",
    "model_metadata <- readRDS('", metadata_path, "')\n\n",
    "# Load your new dataset\n",
    "# new_data <- read.delim('your_new_data.tsv')\n",
    "# Ensure columns: ID, Var (Y/N), and gene expression columns\n\n",
    "# Make predictions using the best model: ", best_model, "\n",
    "# predictions <- predict_Mime_model(\n",
    "#   model = trained_model[['", best_model, "']],\n",
    "#   new_data = new_data,\n",
    "#   genes = model_metadata$n_genes\n",
    "# )\n\n",
    "# The trained model was:\n",
    "# - Project: ", project, "\n",
    "# - Tool: ", tool, "\n",
    "# - Training dataset: ", names(list_train_vali_Data)[1], "\n",
    "# - Validation datasets: ", paste(names(list_train_vali_Data)[-1], collapse=", "), "\n",
    "# - Number of genes: ", length(genelist), "\n",
    "# - Training date: ", as.character(Sys.Date()), "\n",
    "# - Mean AUC: ", round(mean_auc_sorted[1], 3), "\n\n"
)

implementation_script_path <- file.path(output_dir, "implement_model.R")
writeLines(implementation_script, implementation_script_path)
cat("  Created implementation script:", implementation_script_path, "\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat("\n==============================================================\n")
cat("Analysis Complete!\n")
cat("==============================================================\n\n")

cat("Output files:\n")
cat("  1. Trained models:", model_save_path, "\n")
cat("  2. Model metadata:", metadata_path, "\n")
cat("  3. AUC results:", auc_path, "\n")
cat("  4. AUC summary:", auc_summary_file, "\n")
cat("  5. Best model info:", best_model_path, "\n")
cat("  6. Implementation script:", implementation_script_path, "\n")
cat("  7. Visualizations:\n")
cat("     - auc_distribution_all_datasets.pdf\n")
cat("     - roc_curves_all_models.pdf\n")
cat("     - performance_comparison.pdf\n\n")

cat("Best model:", best_model, "with mean AUC =", round(mean_auc_sorted[1], 3), "\n\n")

cat("To use the trained model on new datasets:\n")
cat("  1. Load the model: readRDS('", model_save_path, "')\n")
cat("  2. Follow the implementation script: ", implementation_script_path, "\n")
cat("  3. Apply to new data using the predict function\n\n")

sink()
sink(type="message")
close(log)
