#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


suppressPackageStartupMessages({
    library(clusterProfiler)
    library(DOSE)
    library(dplyr)
})



enrichment_rds <- readRDS(snakemake@input[['enrichment']])
discovery_go_output <- snakemake@output[['discovery_go']]
discovery_kegg_output <- snakemake@output[['discovery_kegg']]
discovery_others_output <- snakemake@output[['discovery_others']]



discovery_others <- rbind(
    enrichment_rds$up_ora$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='WikiPathway'),
    enrichment_rds$up_ora$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Disease Ontology'),
    enrichment_rds$up_ora$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Network of Cancer Gene'),
    enrichment_rds$up_ora$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Disease Gene Network'),
    enrichment_rds$down_ora$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='WikiPathway'),
    enrichment_rds$down_ora$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Disease Ontology'),
    enrichment_rds$down_ora$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Network of Cancer Gene'),
    enrichment_rds$down_ora$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Disease Gene Network')
)

discovery_kegg <- rbind(
    enrichment_rds$up_ora$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='KEGG'),
    enrichment_rds$down_ora$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='KEGG')
)

discovery_go <- rbind(
    enrichment_rds$up_ora$go_bp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO BP'),
    enrichment_rds$up_ora$go_mf %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO MF'),
    enrichment_rds$up_ora$go_cc %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO CC'),
    enrichment_rds$down_ora$go_bp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO BP'),
    enrichment_rds$down_ora$go_mf %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO MF'),
    enrichment_rds$down_ora$go_cc %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO CC')
)






write.table(validation_go,validation_go_output,quote=FALSE,sep='\t',row.names=F)
write.table(validation_kegg,validation_kegg_output,quote=FALSE,sep='\t',row.names=F)
write.table(validation_others,validation_others_output,quote=FALSE,sep='\t',row.names=F)
