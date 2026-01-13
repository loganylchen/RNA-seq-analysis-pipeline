#!/usr/bin/env Rscript
# Differential expression analysis using edgeR

# Set up log file sink BEFORE loading libraries
log_file <- snakemake@log[[1]]
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
  library(edgeR)
  library(dplyr)
  library(limma)
})

# Read inputs
cat("Reading parameters and inputs...\n")
project<- snakemake@params[["project"]]
dataset<- snakemake@params[["dataset"]]
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
samples<-snakemake@params[["samples"]]
counts <- snakemake@input[["counts"]]
design_string<-snakemake@params[["design"]]
# Read outputs


deg_rds<-snakemake@output[["deg_rds"]]
deg_tsv<-snakemake@output[["deg_tsv"]]


cat("=== Limma Trend Analysis ===\n")
cat("Preparing coldata...\n")
cat("Preparing coldata...\n")
coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t',) %>%
            dplyr::filter(dataset_id==dataset)


cts <- read.table(counts, header=TRUE, check.names=FALSE,sep='\t',row.names=1)

limma_trend_pipeline <- function(design_string,count,coldata,
                            condition_col,
                            case_condition,
                            control_condition, parallel=TRUE) {
    cat("\n=== Limma Trend Pipeline Start ===\n")
    cat("Sample information:\n")
    cat("  Total samples:", nrow(coldata), "\n")
    cat("  Columns available:", paste(names(coldata), collapse=", "), "\n")
    cat("Conditions:", case_condition, "vs", control_condition, "\n")

    condition <- factor(coldata[[condition_col]], levels=c(control_condition, case_condition))
    cat("Condition levels:", levels(condition), "\n")
    cat("Condition counts:\n")
    print(table(condition))

    # Build design formula with optional covariates
    if(design_string == ""){
        formula_str <- "~ condition"
    }else{
        # design_string contains the covariate column name
        formula_str <- paste0("~", design_string, " + condition")
    }
    cat("Design formula:", formula_str, "\n")

    cts <- count[,rownames(coldata)]
    cat("Count matrix dimensions:", nrow(cts), "genes x", ncol(cts), "samples\n")
    cat("Count matrix summary:\n")
    print(summary(as.vector(cts)))

    y <- DGEList(counts = cts)
    cat("Original DGEList:", nrow(y), "genes\n")

    # Add condition to coldata for model.matrix to use
    coldata$condition <- condition

    design <- model.matrix(as.formula(formula_str), data=coldata)
    cat("Design matrix:\n")
    print(design)

    keep <- filterByExpr(y, design)
    cat("Genes passing filterByExpr:", sum(keep), "/", length(keep), "\n")
    cat("Filtering percentage:", round(sum(keep)/length(keep)*100, 2), "%\n")

    y <- y[keep, , keep.lib.sizes = FALSE]
    cat("DGEList after filtering:", nrow(y), "genes\n")

    y <- calcNormFactors(y)
    cat("Normalization factors calculated\n")
    cat("Library sizes:\n")
    print(y$samples$lib.size)
    cat("Normalized library sizes (effective):\n")
    print(y$samples$lib.size * y$samples$norm.factors)

    logCPM <- cpm(y,log=TRUE,prior.count=2)
    cat("logCPM matrix dimensions:", nrow(logCPM), "genes x", ncol(logCPM), "samples\n")
    cat("logCPM summary:\n")
    print(summary(as.vector(logCPM)))
    cat("logCPM matrix (first few rows):\n")
    print(head(logCPM))

    fit <- lmFit(logCPM, design)
    cat("Linear model fit completed\n")

    fit <- eBayes(fit, trend = TRUE)
    cat("Empirical Bayes moderation with trend completed\n")

    res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    res$gene_id <- rownames(res)

    # Reorder columns
    res_df <- res[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    colnames(res_df) <- c("gene_id", "log2FoldChange", "average_expression", "t_statistic", "pvalue", "padj", "B_statistic")

    cat("Results summary:\n")
    cat("  Total genes tested:", nrow(res_df), "\n")
    cat("  Significant genes (p_adj < 0.05):", sum(res_df$padj < 0.05), "\n")
    cat("  Significant genes (p_adj < 0.01):", sum(res_df$padj < 0.01), "\n")
    cat("  Significant genes (p_adj < 0.001):", sum(res_df$padj < 0.001), "\n")
    cat("  |log2FC| > 1:", sum(abs(res_df$log2FoldChange) > 1), "\n")
    cat("  |log2FC| > 2:", sum(abs(res_df$log2FoldChange) > 2), "\n")
    cat("Top 10 upregulated genes:\n")
    print(head(res_df[order(res_df$log2FoldChange, decreasing=TRUE), c("gene_id", "log2FoldChange", "padj")], 10))
    cat("Top 10 downregulated genes:\n")
    print(head(res_df[order(res_df$log2FoldChange), c("gene_id", "log2FoldChange", "padj")], 10))
    cat("=== Limma Trend Pipeline End ===\n\n")

    return(res_df)
}

cat("Processing  set...\n")
res <- limma_trend_pipeline(design_string, cts, coldata,
                                condition_col= "condition",
                                case_condition,
                                control_condition)

cat("Saving  results...\n")


                      
saveRDS(res, deg_rds)
write.table(res, deg_tsv, sep='\t', quote=FALSE, row.names=FALSE)



# Close sink connections
sink()
sink(type = "message")
close(log)