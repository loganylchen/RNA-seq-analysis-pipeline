#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(data.table)
    library(readr)
})



making_TPM_matrix <- function(tpm_matrix,case_samples,control_samples, discovery_tpm_matrix) {
    
    df_merge <- read.csv(tpm_matrix,sep='\t')
    df_merge <- df_merge[,c(control_samples,case_samples)]
    write.table(df_merge,discovery_tpm_matrix,quote=F,sep='\t',row.names=T)
}

making_TPM_matrix(
snakemake@input[['tpm_matrix']], 
unname(unlist(snakemake@params[['case_samples']])),
unname(unlist(snakemake@params[['control_samples']])),
snakemake@output[['discovery_tpm_matrix']]
                    )


