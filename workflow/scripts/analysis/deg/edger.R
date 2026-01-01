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


cat("=== edgeR Analysis ===\n")
cat("Preparing coldata...\n")
coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t',)
coldata_discovery <- coldata %>% 
                    dplyr::filter(sample_type==discovery_sample_type)
coldata_validation <- coldata %>% 
                    dplyr::filter(sample_type != discovery_sample_type)

cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')

edger_pipeline <- function(count,coldata,
                            condition_col,
                            case_condition, 
                            control_condition, parallel=TRUE) {
    group <- coldata %>% mutate(condition = factor({{ condition_col }},levels=c(case_condition,control_condition))) %>% pull(condition)
    cts <- count[,rownames(coldata)]
    y <- DGEList(counts = cts, 
                group = group)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    head(y)
    design <- model.matrix(~ group)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    res <- topTags(qlf, n = Inf, adjust.method = "BH")
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- res_df[, c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")]
    colnames(res_df) <- c("gene_id", "log2FoldChange", "logCPM", "F_statistic", "pvalue", "padj")
    res_df <- res_df[order(res_df$pvalue), ]
    return(res_df)
}

cat("Processing discovery set...\n")
res_discovery <- edger_pipeline(cts, coldata_discovery,
                                condition_col= "condition",
                                case_condition,
                                control_condition)
saveRDS(res_discovery, discovery_deg_rds)
write.table(res_discovery, discovery_deg_tsv, sep='\t', quote=FALSE, row.names=FALSE)

cat("Processing validation set...\n")
res_validation <- edger_pipeline(cts, coldata_validation,
                                condition_col= "condition",
                                case_condition,
                                control_condition)
saveRDS(res_validation, validation_deg_rds)
write.table(res_validation, validation_deg_tsv, sep='\t', quote=FALSE, row.names=FALSE)



# Close sink connections
sink()
sink(type = "message")
close(log)