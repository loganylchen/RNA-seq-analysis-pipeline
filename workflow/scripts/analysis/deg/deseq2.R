# Set up log file sink BEFORE loading libraries
log_file <- snakemake@log[[1]]
log <- file(log_file, open = "wt")
sink(log)
sink(log, type = "message")

suppressPackageStartupMessages({
    library(dplyr)
    library(DESeq2)
})

cat("Reading parameters and inputs...\n")
project<- snakemake@params[["project"]]
dataset<- snakemake@params[["dataset"]]
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
samples<-snakemake@params[["samples"]]
design_string<-snakemake@params[["design"]]
counts <- snakemake@input[["counts"]]

# output
cat("Preparing outputs...\n")
count_rds<-snakemake@output[["count_rds"]]
vst_rds<-snakemake@output[["vst_rds"]]
deg_rds<-snakemake@output[["deg_rds"]]
deg_tsv<-snakemake@output[["deg_tsv"]]

cat("=== DESeq2 Analysis ===\n")
cat("Preparing coldata...\n")
coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t',) %>%
            dplyr::filter(dataset_id==dataset)

cts <- read.table(counts, header=TRUE, row.names=1, check.names=FALSE,sep='\t')

deseq2_pipeline <- function(design_string,count,
                            coldata,
                            condition_col,
                            case_condition,
                            control_condition, parallel=TRUE){
    cat("\n=== DESeq2 Pipeline Start ===\n")
    cat("Sample information:\n")
    cat("  Total samples:", nrow(coldata), "\n")
    cat("  Columns available:", paste(names(coldata), collapse=", "), "\n")
    cat("Conditions:", case_condition, "vs", control_condition, "\n")

    # Create condition factor with proper levels
    condition <- factor(coldata[[condition_col]], levels=c(control_condition, case_condition))
    cat("Condition levels:", levels(condition), "\n")
    cat("Condition counts:\n")
    print(table(condition))

    coldata$condition <- condition

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

    design <- as.formula(formula_str)
    dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=design)

    cat("DESeqDataSet created:\n")
    cat("  Genes:", nrow(dds), "\n")
    cat("  Samples:", ncol(dds), "\n")

    dds <- DESeq(dds, parallel=parallel)
    cat("DESeq completed\n")
    cat("Dispersion estimates:\n")
    print(summary(mcols(dds)$dispGeneEst))
    cat("Dispersion outlier:", sum(mcols(dds)$dispOutlier), "\n")

    vsd <- vst(dds)
    cat("VST transformation completed\n")
    cat("VST dimensions:", nrow(vsd), "genes x", ncol(vsd), "samples\n")

    res_dds <- results(dds, contrast=c("condition",case_condition,control_condition),parallel=parallel)
    res_df <- as.data.frame(res_dds)
    res_df$gene_id <- rownames(res_df)

    # Reorder columns to match edgeR output format
    res_df <- res_df[, c("gene_id", "log2FoldChange", "baseMean", "lfcSE", "stat", "pvalue", "padj")]
    colnames(res_df) <- c("gene_id", "log2FoldChange", "baseMean", "lfcSE", "stat", "pvalue", "padj")
    res_df <- res_df[order(res_df$pvalue), ]

    cat("Results summary:\n")
    cat("  Total genes tested:", nrow(res_df), "\n")
    cat("  Significant genes (p_adj < 0.05):", sum(res_df$padj < 0.05, na.rm=TRUE), "\n")
    cat("  Significant genes (p_adj < 0.01):", sum(res_df$padj < 0.01, na.rm=TRUE), "\n")
    cat("  Significant genes (p_adj < 0.001):", sum(res_df$padj < 0.001, na.rm=TRUE), "\n")
    cat("  |log2FC| > 1:", sum(abs(res_df$log2FoldChange) > 1, na.rm=TRUE), "\n")
    cat("  |log2FC| > 2:", sum(abs(res_df$log2FoldChange) > 2, na.rm=TRUE), "\n")
    cat("Top 10 upregulated genes:\n")
    print(head(res_df[order(res_df$log2FoldChange, decreasing=TRUE), c("gene_id", "log2FoldChange", "padj")], 10))
    cat("Top 10 downregulated genes:\n")
    print(head(res_df[order(res_df$log2FoldChange), c("gene_id", "log2FoldChange", "padj")], 10))
    cat("=== DESeq2 Pipeline End ===\n\n")

    return(list(dds=dds, vsd=vsd, res_dds=res_dds))
}

save_list <- function(deseq2_list,
                          dds_file,
                          vsd_file,
                          res_file,
                          tsv_file){
    cat("\nSaving results...\n")
    cat("  Saving DESeqDataSet to:", dds_file, "\n")
    saveRDS(deseq2_list$dds, file=dds_file)

    cat("  Saving VST to:", vsd_file, "\n")
    saveRDS(deseq2_list$vsd, file=vsd_file)

    cat("  Saving results to:", res_file, "\n")
    saveRDS(deseq2_list$res_dds, file=res_file)

    cat("  Writing TSV to:", tsv_file, "\n")
    write.table(deseq2_list$res_dds,tsv_file,sep='\t',quote=FALSE)
}

# Setup parallelization
parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
    cat("Parallel execution enabled with", snakemake@threads, "threads\n")
} else {
    cat("Single-threaded execution\n")
}

cat("\nRunning DESeq2...\n")
deseq2_data <- deseq2_pipeline(design_string=design_string,
                                count=cts,
                                   coldata=coldata,
                                   condition_col="condition",
                                   case_condition=case_condition,
                                   control_condition=control_condition,
                                   parallel=parallel)

cat("Saving results...\n")
save_list(deseq2_data,
          dds_file=count_rds,
          vsd_file=vst_rds,
          res_file=deg_rds,
          tsv_file=deg_tsv)

cat("\nDESeq2 analysis completed successfully.\n")
cat("=== DESeq2 Analysis Complete ===\n")

# Close sink connections
sink()
sink(type="message")
close(log)