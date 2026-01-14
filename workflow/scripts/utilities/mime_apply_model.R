#!/usr/bin/env Rscript
# Apply trained Mime model to new datasets
# This script loads a saved Mime model and applies it to new data for prediction

suppressPackageStartupMessages({
    library(Mime)
    library(dplyr)
})

cat("==============================================================\n")
cat("Mime Model Implementation\n")
cat("==============================================================\n\n")

# Get parameters from Snakemake
model_file <- snakemake@input[["model"]]
new_data_file <- snakemake@input[["new_data"]]
output_file <- snakemake@output[["predictions"]]

cat("Model file:", model_file, "\n")
cat("New data file:", new_data_file, "\n")
cat("Output file:", output_file, "\n\n")

# ============================================================================
# LOAD TRAINED MODEL
# ============================================================================

cat("Loading trained model...\n")
trained_model <- readRDS(model_file)
model_metadata <- attr(trained_model, "metadata")

if (is.null(model_metadata)) {
    # Try to load from separate file
    metadata_file <- sub("trained_models\\.rds$", "model_metadata.rds", model_file)
    if (file.exists(metadata_file)) {
        model_metadata <- readRDS(metadata_file)
    }
}

cat("  Model loaded successfully\n")
if (!is.null(model_metadata)) {
    cat("  Project:", model_metadata$project, "\n")
    cat("  Tool:", model_metadata$tool, "\n")
    cat("  Training dataset:", model_metadata$train_dataset, "\n")
    cat("  Training date:", model_metadata$train_date, "\n")
}
cat("\n")

# ============================================================================
# LOAD NEW DATA
# ============================================================================

cat("Loading new data...\n")
new_data <- read.delim(new_data_file, header=TRUE, row.names=1, stringsAsFactors=FALSE)
cat("  Dimensions:", nrow(new_data), "samples x", ncol(new_data), "columns\n")

# Check required columns
if (!"ID" %in% colnames(new_data)) {
    stop("ERROR: New data must have an 'ID' column (sample IDs)")
}

if (!"Var" %in% colnames(new_data)) {
    cat("  WARNING: No 'Var' column found (true labels)\n")
    cat("  Predictions will be generated without validation\n")
}

# ============================================================================
# APPLY MODEL
# ============================================================================

cat("\nApplying model to new data...\n")

predictions_list <- list()

for (model_name in names(trained_model)) {
    cat("  Model:", model_name, "\n")

    model_obj <- trained_model[[model_name]]

    # Get prediction function from Mime
    # This depends on the specific Mime implementation
    tryCatch({
        # Extract the trained model object
        if (!is.null(model_obj$trained_model)) {
            final_model <- model_obj$trained_model
        } else {
            final_model <- model_obj
        }

        # Make predictions
        # Note: This is a simplified example
        # Actual implementation depends on Mime's predict method
        predictions <- predict(final_model, new_data, type = "prob")

        predictions_list[[model_name]] <- predictions
        cat("    Predictions generated\n")

    }, error = function(e) {
        cat("    ERROR:", conditionMessage(e), "\n")
    })
}

# ============================================================================
# SAVE PREDICTIONS
# ============================================================================

cat("\nSaving predictions...\n")

# Save as RDS
saveRDS(predictions_list, file = output_file)
cat("  Saved to:", output_file, "\n")

# Also save as TSV for each model
for (model_name in names(predictions_list)) {
    pred <- predictions_list[[model_name]]

    if (is.matrix(pred) || is.data.frame(pred)) {
        pred_file <- sub("\\.rds$", paste0("_", model_name, ".tsv"), output_file)
        write.table(pred, pred_file, sep="\t", quote=FALSE, col.names=NA)
        cat("  ", model_name, ":", pred_file, "\n")
    }
}

cat("\n==============================================================\n")
cat("Prediction Complete!\n")
cat("==============================================================\n\n")
