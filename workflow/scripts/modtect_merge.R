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


merge_modtect_files <- function(mod_tect_files, merged_file) {
    df_list <- list()
    for(f in mod_tect_files){
        sample_name <- basename(dirname(f))
        message(sample_name,':',f)
        tmp_df <- read_tsv(f, comment = "#", progress = FALSE,header=FALSE) %>% 
                        dplyr::mutate(Sample=sample_name)
        colnames(tmp_df)<-c(
            'chrom',	'position',	'reference_nt',
            'mono-alleleLogP'	,'bi-alleleLogP',	'tri-alleleLogP'	,
            'tetra-alleleLogP',	'variant_proportion',	'depth'	,
            'deletionCount',	'deletion_proportion',	'refUpperCount',
            'refLowerCount'	,'altUpperCount'	,'altLowerCount',
            'refPlusProp',	'altPlusProp'	,'A_count'	,
            'T_count'	,'G_count'	,'C_count',
            'a_count',	't_count',	'g_count',
            'c_count',	'has_reference_nt'	,'types_of_nt',
            'ModTect_score',	'variant_readposn_median_fwd',	'variant_readposn_median_rev',
            'median_abs_dev_fwd',	'median_abs_dev_rev', 'Sample'
        )
        message('reading:',f)
        df_list[[sample_name]] <- tmp_df
    }
    message('merging')
    df_merge <- data.table::rbindlist(df_list) %>% 
                tidyr::pivot_wider(id_cols=c(chrom,position,reference_nt),names_from=Sample,values_from=ModTect_score)
    message(head(df_merge))
    write.table(df_merge,merged_file,quote=F,sep='\t',row.names=F)
}

making_count_matrix(unname(unlist(snakemake@input)), 
                    snakemake@output[['output']]
                    )


