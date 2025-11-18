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



making_TPM_matrix <- function(fc_count_files, count_matrix) {
    df_list <- list()
    for(f in fc_count_files){
        sample_name <- gsub('.txt','',basename(f))
        message(sample_name,':',f)
        tmp_df <- read_tsv(f, comment = "#", progress = FALSE) %>% 
                        dplyr::mutate(Sample=sample_name)
        colnames(tmp_df) <- c('Geneid','Chr','Start','End','Strand','Length','Count','Sample')
        tmp_df %>% 
        dplyr::mutate(length_kb=Length/1000) %>%
        dplyr::mutate(rpk=Count/length_kb) %>%
        dplyr::mutate(tpm=rpk/(sum(rpk)/1e6))
        message('reading:',f)
        df_list[[sample_name]] <- tmp_df
    }
    message('merging')
    df_merge <- data.table::rbindlist(df_list) %>% 
                dplyr::mutate(TPM=tpm)%>%
                tidyr::pivot_wider(id_cols=Geneid,names_from=Sample,values_from=TPM)
    message(head(df_merge))
    write.table(df_merge,count_matrix,quote=F,sep='\t',row.names=F)
}

making_TPM_matrix(unname(unlist(snakemake@input)), 
                    snakemake@output[['count_matrix']]
                    )


