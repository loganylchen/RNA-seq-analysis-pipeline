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


making_count_matrix <- function(fc_count_files, count_matrix,puree_count_matrix) {
    df_list <- list()
    for(f in fc_count_files){
        sample_name <- gsub('.star_counts.txt','',basename(f))
        message(sample_name,':',f)
        tmp_df <- read_tsv(f, comment = "#", progress = FALSE) %>% 
                        dplyr::mutate(sample=sample_name)
        colnames(tmp_df) <- c('Geneid','Chr','Start','End','Strand','Length','Count','Sample')
        message('reading:',f)
        df_list[[sample_name]] <- tmp_df
    }
    message('merging')
    df_merge <- data.table::rbindlist(df_list) %>% 
                tidyr::pivot_wider(id_cols=Geneid,names_from=Sample,values_from=Count)
    write.table(df_merge,count_matrix,quote=F,sep='\t',row.names=F)
    write.table(t(df_merge),puree_count_matrix,quote=F,sep='\t',row.names=T,col.names=F)

}

making_count_matrix(unname(unlist(snakemake@input)), 
                    snakemake@output[['count_matrix']],
                    snakemake@output[['puree_count_matrix']]
                    )


