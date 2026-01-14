#!/usr/bin/env Rscript
# Example script for running Mime analysis with prepared datasets
# This script demonstrates how to use Mime for survival analysis and response prediction

# ============================================================================
# SETUP
# ============================================================================

# Install Mime if not already installed
if (!requireNamespace("Mime", quietly = TRUE)) {
    # Install dependencies
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    dependons <- c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow',
                    'sva', 'ComplexHeatmap')
    for (depen in dependons) {
        if (!requireNamespace(depen, quietly = TRUE))
            BiocManager::install(depen, update = FALSE)
    }

    # Install GitHub packages
    if (!requireNamespace("CoxBoost", quietly = TRUE))
        devtools::install_github("binderh/CoxBoost")

    if (!requireNamespace("fastAdaboost", quietly = TRUE))
        devtools::install_github("souravc83/fastAdaboost")

    # Install Mime
    devtools::install_github("l-magnificence/Mime")
}

library(Mime)

# ============================================================================
# SURVIVAL ANALYSIS EXAMPLE
# ============================================================================

cat("==============================================================\n")
cat("Mime Survival Analysis Example\n")
cat("==============================================================\n\n")

# Load prepared dataset
# This file is created by prepare_mime_dataset.R
load("./survival_dataset.rds")  # Creates list_train_vali_Data

# Load gene list (example: MSigDB Wnt/beta-catenin signaling pathway)
# Replace with your DEG results or pathway of interest
genelist <- c("MYC", "CTNNB1", "JAG2", "NOTCH1", "DLL1", "AXIN2", "PSEN2",
              "FZD1", "NOTCH4", "LEF1", "AXIN1", "NKD1", "WNT5B", "CUL1",
              "JAG1", "MAML1", "KAT2A", "GNAI1", "WNT6", "PTCH1", "NCOR2",
              "DKK4", "HDAC2", "DKK1", "TCF7", "WNT1", "NUMB", "ADAM17",
              "DVL2", "PPARD", "NCSTN", "HDAC5", "CCND2", "FRAT1", "CSNK1E",
              "RBPJ", "FZD8", "TP53", "SKP2", "HEY2", "HEY1", "HDAC11")

# Alternatively, use DEG results from your pipeline
# deg_data <- read.delim("./DEG/deseq2/STAR_FC/discovery_deg.tsv")
# genelist <- deg_data$gene_id[deg_data$padj < 0.05 & deg_data$log2FoldChange > 0]

cat("Dataset structure:\n")
print(str(list_train_vali_Data))
cat("\nGene list length:", length(genelist), "genes\n\n")

# 1. Develop predictive models for prognosis
cat("Step 1: Developing ML models for prognosis...\n")
res <- ML.Dev.Prog.Sig(
    train_data = list_train_vali_Data$Dataset1,
    list_train_vali_Data = list_train_vali_Data,
    unicox.filter.for.candi = TRUE,
    unicox_p_cutoff = 0.05,
    candidate_genes = genelist,
    mode = 'all',
    nodesize = 5,
    seed = 5201314
)

# 2. Visualize C-index across models
cat("\nStep 2: Plotting C-index distribution...\n")
pdf("mime_cindex_distribution.pdf", width = 12, height = 6)
cindex_dis_all(
    res,
    validate_set = names(list_train_vali_Data)[-1],
    order = names(list_train_vali_Data),
    width = 0.35
)
dev.off()

# 3. Plot specific model C-index
cat("\nStep 3: Plotting C-index for best model...\n")
pdf("mime_cindex_best_model.pdf", width = 8, height = 6)
cindex_dis_select(
    res,
    model = "StepCox[forward] + plsRcox",  # Replace with your best model
    order = names(list_train_vali_Data)
)
dev.off()

# 4. Survival curves
cat("\nStep 4: Generating survival curves...\n")
pdf("mime_survival_curves.pdf", width = 12, height = 6)
survplot <- vector("list", length(list_train_vali_Data))
for (i in seq_along(list_train_vali_Data)) {
    survplot[[i]] <- rs_sur(
        res,
        model_name = "StepCox[forward] + plsRcox",
        dataset = names(list_train_vali_Data)[i],
        median.line = "hv",
        cutoff = 0.5,
        conf.int = TRUE,
        xlab = "Day",
        pval.coord = c(1000, 0.9)
    )
}
aplot::plot_list(gglist = survplot, ncol = 2)
dev.off()

# 5. Calculate AUC scores
cat("\nStep 5: Calculating time-dependent AUC...\n")
all.auc.1y <- cal_AUC_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    train_data = list_train_vali_Data[["Dataset1"]],
    inputmatrix.list = list_train_vali_Data,
    mode = 'all',
    AUC_time = 1,
    auc_cal_method = "KM"
)

all.auc.3y <- cal_AUC_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    train_data = list_train_vali_Data[["Dataset1"]],
    inputmatrix.list = list_train_vali_Data,
    mode = 'all',
    AUC_time = 3,
    auc_cal_method = "KM"
)

all.auc.5y <- cal_AUC_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    train_data = list_train_vali_Data[["Dataset1"]],
    inputmatrix.list = list_train_vali_Data,
    mode = 'all',
    AUC_time = 5,
    auc_cal_method = "KM"
)

# 6. Visualize AUC distribution
cat("\nStep 6: Plotting AUC distribution...\n")
pdf("mime_auc_distribution.pdf", width = 12, height = 6)
auc_dis_all(
    all.auc.1y,
    dataset = names(list_train_vali_Data),
    validate_set = names(list_train_vali_Data)[-1],
    order = names(list_train_vali_Data),
    width = 0.35,
    year = 1
)
dev.off()

# 7. ROC curves
cat("\nStep 7: Generating ROC curves...\n")
pdf("mime_roc_curves.pdf", width = 8, height = 6)
roc_vis(
    all.auc.1y,
    model_name = "StepCox[forward] + plsRcox",
    dataset = names(list_train_vali_Data),
    order = names(list_train_vali_Data),
    anno_position = c(0.65, 0.55),
    year = 1
)
dev.off()

# 8. Meta-analysis
cat("\nStep 8: Performing meta-analysis...\n")
unicox.rs.res <- cal_unicox_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    optimal.model = "StepCox[forward] + plsRcox",
    type = 'categorical'
)

metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)

pdf("mime_meta_analysis.pdf", width = 10, height = 6)
meta_unicox_vis(
    metamodel,
    dataset = names(list_train_vali_Data)
)
dev.off()

# ============================================================================
# RESPONSE PREDICTION EXAMPLE
# ============================================================================

cat("\n==============================================================\n")
cat("Mime Response Prediction Example\n")
cat("==============================================================\n\n")

# Load response dataset if available
if (file.exists("./response_dataset.rds")) {
    load("./response_dataset.rds")  # Creates list_train_vali_Data

    cat("Developing response prediction models...\n")
    res.ici <- ML.Dev.Pred.Category.Sig(
        train_data = list_train_vali_Data$training,
        list_train_vali_Data = list_train_vali_Data,
        candidate_genes = genelist,
        methods = c('nb', 'svmRadialWeights', 'rf', 'kknn', 'adaboost',
                    'LogitBoost', 'cancerclass'),
        seed = 5201314,
        cores_for_parallel = 4  # Adjust based on your system
    )

    # Visualize AUC across methods
    pdf("mime_response_auc.pdf", width = 10, height = 6)
    auc_vis_category_all(
        res.ici,
        dataset = c("training", "validation"),
        order = c("training", "validation")
    )
    dev.off()

    # ROC curves for each method
    pdf("mime_response_roc_all.pdf", width = 15, height = 10)
    plot_list <- list()
    methods <- c('nb', 'svmRadialWeights', 'rf', 'kknn', 'adaboost',
                 'LogitBoost', 'cancerclass')
    for (i in methods) {
        plot_list[[i]] <- roc_vis_category(
            res.ici,
            model_name = i,
            dataset = c("training", "validation"),
            order = c("training", "validation"),
            anno_position = c(0.4, 0.25)
        )
    }
    aplot::plot_list(gglist = plot_list, ncol = 3)
    dev.off()
}

# ============================================================================
# CORE FEATURE SELECTION EXAMPLE
# ============================================================================

cat("\n==============================================================\n")
cat("Mime Core Feature Selection Example\n")
cat("==============================================================\n\n")

cat("Performing core feature selection...\n")
res.feature.all <- ML.Corefeature.Prog.Screen(
    InputMatrix = list_train_vali_Data$Dataset1,
    candidate_genes = genelist,
    mode = "all",
    nodesize = 5,
    seed = 5201314
)

# Upset plot of selected features
pdf("mime_core_features_upset.pdf", width = 10, height = 6)
core_feature_select(res.feature.all)
dev.off()

# Feature rank plot
pdf("mime_core_features_rank.pdf", width = 12, height = 8)
core_feature_rank(res.feature.all, top = 20)
dev.off()

cat("\n==============================================================\n")
cat("Mime Analysis Complete!\n")
cat("==============================================================\n\n")

cat("Output files:\n")
cat("  - mime_cindex_distribution.pdf\n")
cat("  - mime_cindex_best_model.pdf\n")
cat("  - mime_survival_curves.pdf\n")
cat("  - mime_auc_distribution.pdf\n")
cat("  - mime_roc_curves.pdf\n")
cat("  - mime_meta_analysis.pdf\n")
cat("  - mime_response_auc.pdf (if response data available)\n")
cat("  - mime_response_roc_all.pdf (if response data available)\n")
cat("  - mime_core_features_upset.pdf\n")
cat("  - mime_core_features_rank.pdf\n\n")

cat("Next steps:\n")
cat("  1. Review the generated plots\n")
cat("  2. Select the best performing model\n")
cat("  3. Extract core features for biomarker discovery\n")
cat("  4. Validate findings in independent datasets\n\n")
