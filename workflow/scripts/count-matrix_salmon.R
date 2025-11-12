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


making_count_matrix <- function(fc_count_files, count_matrix) {
    df_list <- list()
    for(f in fc_count_files){
        sample_name <- basename(dirname(f))
        message(sample_name,':',f)
        tmp_df <- read_tsv(f, comment = "#", progress = FALSE) %>% 
                        dplyr::mutate(sample=sample_name)
        message('reading:',f)
        df_list[[sample_name]] <- tmp_df
    }
    message('merging')
    df_merge <- data.table::rbindlist(df_list) %>% 
                dplyr::mutate(Count=as.integer(NumReads))%>%
                tidyr::pivot_wider(id_cols=Name,names_from=Sample,values_from=Count)
    message(head(df_merge))
    write.table(df_merge,count_matrix,quote=F,sep='\t',row.names=F)
}

making_count_matrix(unname(unlist(snakemake@input)), 
                    snakemake@output[['count_matrix']],
                    )


