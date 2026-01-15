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
    library(Mime1)
    library(dplyr)
    library(readr)
})

cat("==============================================================\n")
cat("Mime Response Prediction Analysis\n")
cat("==============================================================\n\n")

# Get parameters from Snakemake
combined_rds <- snakemake@input[["combined_rds"]]
genelist_file <- snakemake@input[["genelist"]]
output_dir <- snakemake@output[["directory"]]



tool <- snakemake@params[["tool"]]
seed <- snakemake@params[["seed"]]


cat("Parameters:\n")
cat("  Tool:", tool, "\n")
cat("  Combined RDS:", combined_rds, "\n")
cat("  Gene list file:", genelist_file, "\n")
cat("  Output directory:", output_dir, "\n")
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

# Load and filter gene list - get common up-regulated genes from discovery dataset
cat("\nLoading and filtering gene list from discovery dataset...\n")
combined_deg <-  fread(gene_list_file, header=TRUE, stringsAsFactors=FALSE) 
cat("  Loaded", nrow(combined_deg), "genes from combined DEG file\n")

# Filter for common up-regulated genes detected by ALL 4 tools
# up_regulated_count == 4 means detected as up-regulated by all tools
genelist <- combined_deg %>%
    filter(up_regulated_count == 4) %>%
    pull(gene_id)

cat("  Common up-regulated genes (detected by all 4 tools):", length(genelist), "\n")

# Check if we have enough genes
if (length(genelist) == 0) {
    cat("  WARNING: No common up-regulated genes found!\n")
    cat("  Using top 100 up-regulated genes instead...\n")
    genelist <- combined_deg %>%
        arrange(desc(up_regulated_count)) %>%
        head(100) %>%
        pull(gene_id)
    cat("  Using", length(genelist), "genes\n")
}

# Save gene list for reference
genelist_output <- file.path(output_dir, "genelist_used.txt")
writeLines(genelist, genelist_output)
cat("  Saved gene list to:", genelist_output, "\n")

# ============================================================================
# TRAIN RESPONSE PREDICTION MODELS
# ============================================================================

cat("\n==============================================================\n")
cat("Step 1: Training Response Prediction Models\n")
cat("==============================================================\n\n")

cat("Training on:", names(list_train_vali_Data)[1], "\n")
cat("Validating on:", paste(names(list_train_vali_Data), collapse=", "), "\n")
cat("Using", length(genelist), "genes as features\n\n")

# Train models
res.ici <- ML.Dev.Pred.Category.Sig(
    train_data = list_train_vali_Data[[1]],
    test_data_list = list_train_vali_Data,
    sig = genelist,
    methods = c("nb", "svmRadialWeights", "rf", "kknn", "adaboost", "LogitBoost", "cancerclass"),
    seed = seed
)

# Save results
saveRDS(res.ici, file.path(output_dir, "mime_results.rds"))
cat("  Saved results to:", file.path(output_dir, "mime_results.rds"), "\n")




# ============================================================================
# GENERATE SUMMARY STATISTICS
# ============================================================================

# cat("\n==============================================================\n")
# cat("Step 2: Generating Summary Statistics\n")
# cat("==============================================================\n\n")

# # Extract performance metrics
# summary_stats <- data.frame(
#     Dataset = character(),
#     Method = character(),
#     AUC = numeric(),
#     Accuracy = numeric(),
#     Sensitivity = numeric(),
#     Specificity = numeric(),
#     stringsAsFactors = FALSE
# )

# for (dataset_name in names(res.ici$test_performance)) {
#     perf <- res.ici$test_performance[[dataset_name]]
#     for (method_name in names(perf)) {
#         if (!is.null(perf[[method_name]])) {
#             metrics <- perf[[method_name]]
#             summary_stats <- rbind(summary_stats, data.frame(
#                 Dataset = dataset_name,
#                 Method = method_name,
#                 AUC = ifelse(!is.null(metrics$auc), metrics$auc, NA),
#                 Accuracy = ifelse(!is.null(metrics$acc), metrics$acc, NA),
#                 Sensitivity = ifelse(!is.null(metrics$sensitivity), metrics$sensitivity, NA),
#                 Specificity = ifelse(!is.null(metrics$specificity), metrics$specificity, NA),
#                 stringsAsFactors = FALSE
#             ))
#         }
#     }
# }

# # Save summary table
# summary_file <- file.path(output_dir, "benchmark_summary.tsv")
# write.table(summary_stats, summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
# cat("  Saved summary to:", summary_file, "\n")

# # Print summary
# cat("\nPerformance Summary:\n")
# print(summary_stats)
# cat("\n")

# # ============================================================================
# # SAVE TRAINED MODELS
# # ============================================================================

# cat("\n==============================================================\n")
# cat("Step 3: Saving Trained Models\n")
# cat("==============================================================\n\n")

# # Save trained models for each method
# model_dir <- file.path(output_dir, "models")
# dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# trained_models <- list()
# for (method in methods) {
#     if (!is.null(res.ici$trained_model[[method]])) {
#         model_file <- file.path(model_dir, paste0(method, "_model.rds"))
#         saveRDS(res.ici$trained_model[[method]], model_file)
#         trained_models[[method]] <- model_file
#         cat("  Saved", method, "model to:", model_file, "\n")
#     }
# }

# # Save model list for easy loading
# model_list_file <- file.path(output_dir, "trained_models.rds")
# saveRDS(trained_models, model_list_file)
# cat("\n  Saved model list to:", model_list_file, "\n")

# # ============================================================================
# # CREATE IMPLEMENTATION GUIDE
# # ============================================================================

# cat("\n==============================================================\n")
# cat("Step 4: Creating Implementation Guide\n")
# cat("==============================================================\n\n")

# guide_file <- file.path(output_dir, "implementation_guide.txt")

# guide_content <- paste0(
#     "Mime Model Implementation Guide\n",
#     "================================\n\n",
#     "Project: ", project, "\n",
#     "Quantification Tool: ", tool, "\n",
#     "Training Dataset: Discovery (Dataset1)\n",
#     "Number of Features (Genes): ", length(genelist), "\n",
#     "Methods Trained: ", paste(methods, collapse = ", "), "\n",
#     "Seed: ", seed, "\n\n",
#     "Feature Selection:\n",
#     "  - Common up-regulated genes from discovery dataset\n",
#     "  - Detected by all 4 DEG tools (DESeq2, edgeR, limma-trend, limma-voom)\n",
#     "  - log2FC >= ", log2fc_threshold, ", padj < ", padj_threshold, "\n\n",
#     "Trained Models:\n"
# )

# for (method in names(trained_models)) {
#     guide_content <- paste0(guide_content, "  - ", method, ": ", trained_models[[method]], "\n")
# }

# guide_content <- paste0(guide_content, "\n",
#     "Best Performing Model:\n",
#     "  See benchmark_summary.tsv for performance metrics\n\n",
#     "To Apply Models to New Data:\n",
#     "  1. Load the trained model:\n",
#     "       model <- readRDS('", model_list_file, "')\n",
#     "       trained_model <- readRDS(model$<method_name>)\n\n",
#     "  2. Prepare your new data:\n",
#     "       - Must have the same ", length(genelist), " genes as features\n",
#     "       - Expression values should be TPM or normalized counts\n",
#     "       - Samples should be in rows, genes in columns\n\n",
#     "  3. Use Mime::ML.Pred.Category.Sig.Single to predict:\n",
#     "       predictions <- ML.Pred.Category.Sig.Single(\n",
#     "           test_data = your_data,\n",
#     "           sig = genelist,\n",
#     "           model = trained_model\n",
#     "       )\n\n"
# )

# writeLines(guide_content, guide_file)
# cat("  Saved implementation guide to:", guide_file, "\n")

# # ============================================================================
# # CREATE VISUALIZATIONS
# # ============================================================================

# cat("\n==============================================================\n")
# cat("Step 5: Creating Visualizations\n")
# cat("==============================================================\n\n")

# # Try to create plots (may fail if ggplot2 not available)
# tryCatch({
#     library(ggplot2)
#     library(gridExtra)

#     # AUC comparison plot
#     auc_data <- summary_stats[!is.na(summary_stats$AUC), ]
#     if (nrow(auc_data) > 0) {
#         p_auc <- ggplot(auc_data, aes(x = Method, y = AUC, fill = Dataset)) +
#             geom_bar(stat = "identity", position = "dodge") +
#             geom_text(aes(label = round(AUC, 3)), position = position_dodge(width = 0.9), vjust = -0.25) +
#             theme_minimal() +
#             labs(title = "AUC Comparison Across Methods and Datasets",
#                  y = "AUC") +
#             ylim(0.5, 1)

#         ggsave(file.path(output_dir, "auc_comparison.pdf"), p_auc, width = 10, height = 6)
#         ggsave(file.path(output_dir, "auc_comparison.png"), p_auc, width = 10, height = 6, dpi = 300)
#         cat("  Saved AUC comparison plot\n")
#     }

#     cat("\n")
# }, error = function(e) {
#     cat("  Note: Could not create visualizations (ggplot2 may not be available)\n")
#     cat("  Error:", conditionMessage(e), "\n\n")
# })



# auc_vis_category_all(res.ici,dataset = c("training","validation"),
#                      order= c("training","validation"))
# plot_list<-list()
# methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
# for (i in methods) {
#   plot_list[[i]]<-roc_vis_category(res.ici,model_name = i,dataset = c("training","validation"),
#                                    order= c("training","validation"),
#                                    anno_position=c(0.4,0.25))
# }
# aplot::plot_list(gglist=plot_list,ncol=3)
# auc.other.pre <- cal_auc_previous_sig(list_train_vali_Data = list_train_vali_Data,seed = 5201314,
#                                       train_data = list_train_vali_Data$training,
#                                       cores_for_parallel = 32)
# auc_category_comp(res.ici,
#                   auc.other.pre,
#                   model_name="svmRadialWeights",
#                   dataset=names(list_train_vali_Data))
# load("./Example.cohort.Rdata")
# load("./genelist.Rdata")
# res.feature.all <- ML.Corefeature.Prog.Screen(InputMatrix = list_train_vali_Data$Dataset1,
#                                             candidate_genes = genelist,
#                                             mode = "all",nodesize =5,seed = 5201314 )
# core_feature_select(res.feature.all)
# core_feature_rank(res.feature.all, top=20)
# dataset_col<-c("#3182BDFF","#E6550DFF")
# corplot <- list()
# for (i in c(1:2)) {
#   print(corplot[[i]]<-cor_plot(list_train_vali_Data[[i]],
#                                dataset=names(list_train_vali_Data)[i],
#                                color = dataset_col[i],
#                                feature1="PSEN2",
#                                feature2="WNT5B",
#                                method="pearson"))
# }
# aplot::plot_list(gglist=corplot,ncol=2)
# survplot <- vector("list",2) 
# for (i in c(1:2)) {
#   print(survplot[[i]]<-core_feature_sur("PSEN2", 
#                                         InputMatrix=list_train_vali_Data[[i]],
#                                         dataset = names(list_train_vali_Data)[i],
#                                         #color=c("blue","green"),
#                                         median.line = "hv",
#                                         cutoff = 0.5,
#                                         conf.int = T,
#                                         xlab="Day",pval.coord=c(1000,0.9)))
# }
# aplot::plot_list(gglist=survplot,ncol=2)

cat("==============================================================\n")
cat("Mime Analysis Complete!\n")
cat("==============================================================\n\n")

sink()
sink(type="message")
close(log)
