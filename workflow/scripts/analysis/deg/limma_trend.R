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
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
discovery_sample_type<-snakemake@params[["discovery_sample_type"]]
samples<-snakemake@params[["samples"]]
counts <- snakemake@input[["counts"]]

# Read outputs

discovery_deg_rds<-snakemake@output[["discovery_deg_rds"]]
validation_deg_rds<-snakemake@output[["validation_deg_rds"]]
discovery_deg_tsv<-snakemake@output[["discovery_deg_tsv"]]
validation_deg_tsv<-snakemake@output[["validation_deg_tsv"]]


cat("=== Limma Trend Analysis ===\n")
cat("Preparing coldata...\n")
coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t',)
coldata_discovery <- coldata %>% 
                    dplyr::filter(sample_type==discovery_sample_type)
coldata_validation <- coldata %>% 
                    dplyr::filter(sample_type != discovery_sample_type)

cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')

limma_trend_pipeline <- function(count,coldata,
                            condition_col,
                            case_condition,
                            control_condition, parallel=TRUE) {
    cat("\n=== Limma Trend Pipeline Start ===\n")
    cat("Sample information:\n")
    cat("  Total samples:", nrow(coldata), "\n")
    cat("  Columns available:", paste(names(coldata), collapse=", "), "\n")
    cat("Conditions:", case_condition, "vs", control_condition, "\n")

    group <- coldata %>% mutate(condition = factor({{ condition_col }},levels=c(case_condition,control_condition))) %>% pull(condition)
    cat("Group levels:", levels(group), "\n")
    cat("Group counts:\n")
    print(table(group))

    cts <- count[,rownames(coldata)]
    cat("Count matrix dimensions:", nrow(cts), "genes x", ncol(cts), "samples\n")
    cat("Count matrix summary:\n")
    print(summary(as.vector(cts)))

    y <- DGEList(counts = cts,
                group = group)
    cat("Original DGEList:", nrow(y), "genes\n")

    design <- model.matrix(~ group)
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

cat("Processing discovery set...\n")
res_discovery <- limma_trend_pipeline(cts, coldata_discovery,
                                condition_col= "condition",
                                case_condition,
                                control_condition)
saveRDS(res_discovery, discovery_deg_rds)
write.table(res_discovery, discovery_deg_tsv, sep='\t', quote=FALSE, row.names=FALSE)

cat("Processing validation set...\n")
res_validation <- limma_trend_pipeline(cts, coldata_validation,
                                condition_col= "condition",
                                case_condition,
                                control_condition)
saveRDS(res_validation, validation_deg_rds)
write.table(res_validation, validation_deg_tsv, sep='\t', quote=FALSE, row.names=FALSE)



# Close sink connections
sink()
sink(type = "message")
close(log)