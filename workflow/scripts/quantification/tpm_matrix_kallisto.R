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


making_TPM_matrix <- function(kallisto_files, tpm_matrix) {
    df_list <- list()
    for(f in kallisto_files){
        sample_name <- basename(dirname(f))
        message(sample_name,':',f)
        tmp_df <- read_tsv(f, comment = "#", progress = FALSE) %>%
                        dplyr::select(Name, TPM) %>%
                        dplyr::mutate(Sample=sample_name)
        message('reading:',f)
        df_list[[sample_name]] <- tmp_df
    }
    message('merging')
    df_merge <- data.table::rbindlist(df_list) %>%
                rename(target_id = Name) %>%
                tidyr::pivot_wider(id_cols=target_id, names_from=Sample, values_from=TPM)
    message(head(df_merge))
    write.table(df_merge, tpm_matrix, quote=F, sep='\t', row.names=F)
}

making_TPM_matrix(unname(unlist(snakemake@input)),
                    snakemake@output[['tpm_matrix']]
                    )
