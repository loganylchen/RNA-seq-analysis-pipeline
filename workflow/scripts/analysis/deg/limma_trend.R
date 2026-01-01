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
    group <- coldata %>% mutate(condition = factor({{ condition_col }},levels=c(case_condition,control_condition))) %>% pull(condition)
    cts <- count[,rownames(coldata)]
    y <- DGEList(counts = cts,
                group = group)
    design <- model.matrix(~ group)
    keep <- filterByExpr(y, design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    logCPM <- cpm(y,log=TRUE,prior.count=2)
    head(logCPM)
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend = TRUE)
    res <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
    res$gene_id <- rownames(res)

    # Reorder columns
    res_df <- res[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    colnames(res_df) <- c("gene_id", "log2FoldChange", "average_expression", "t_statistic", "pvalue", "padj", "B_statistic")
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