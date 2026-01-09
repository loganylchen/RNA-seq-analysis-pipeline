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
dataset<- snakemake@params[["dataset"]]
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
samples<-snakemake@params[["samples"]]
counts <- snakemake@input[["counts"]]
design_string<-snakemake@params[["design"]]
# Read outputs

deg_rds<-snakemake@output[["deg_rds"]]
deg_tsv<-snakemake@output[["deg_tsv"]]


cat("=== edgeR Analysis ===\n")
cat("Preparing coldata...\n")
cat("Preparing coldata...\n")
coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t',) %>%
            dplyr::filter(dataset_id==dataset)


cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')

edger_pipeline <- function(design_string,count,coldata,
                            condition_col,
                            case_condition,
                            control_condition, parallel=TRUE) {
    cat("\n=== edgeR Pipeline Start ===\n")
    cat("Sample information:\n")
    cat("  Total samples:", nrow(coldata), "\n")
    cat("  Columns available:", paste(names(coldata), collapse=", "), "\n")
    cat("Conditions:", case_condition, "vs", control_condition, "\n")

    condition <- factor(coldata[[condition_col]], levels=c(control_condition, case_condition))
    cat("Condition levels:", levels(condition), "\n")
    cat("Condition counts:\n")
    print(table(condition))

    cts <- count[,rownames(coldata)]
    cat("Count matrix dimensions:", nrow(cts), "genes x", ncol(cts), "samples\n")
    cat("Count matrix summary:\n")
    print(summary(as.vector(cts)))

    y <- DGEList(counts = cts)
    cat("Original DGEList:", nrow(y), "genes\n")

    design <- model.matrix(as.formula(design_string))
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

cat("Processing  set...\n")
res <- edger_pipeline(design_string,cts, coldata,
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