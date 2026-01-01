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
    cat("\n=== edgeR Pipeline Start ===\n")
    cat("Sample information:\n")
    cat("  Total samples:", nrow(coldata), "\n")
    cat("  Columns available:", paste(names(coldata), collapse=", "), "\n")
    cat("Conditions:", case_condition, "vs", control_condition, "\n")

    group <- factor(coldata[[condition_col]], levels=c(case_condition, control_condition))
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

    cat("DGEList object (first few rows):\n")
    print(head(y))

    y <- estimateDisp(y, design)
    cat("Dispersion estimated\n")
    cat("Dispersion summary:\n")
    print(summary(y$tagwise.dispersion))
    cat("Common dispersion:", y$common.dispersion, "\n")

    fit <- glmQLFit(y, design)
    cat("GLM fit completed\n")

    qlf <- glmQLFTest(fit, coef = 2)
    cat("QLF test completed\n")

    res <- topTags(qlf, n = Inf, adjust.method = "BH")
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- res_df[, c("gene_id", "logFC", "logCPM", "F", "PValue", "FDR")]
    colnames(res_df) <- c("gene_id", "log2FoldChange", "logCPM", "F_statistic", "pvalue", "padj")
    res_df <- res_df[order(res_df$pvalue), ]

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
    cat("=== edgeR Pipeline End ===\n\n")

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