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
discovery_gsea_output <- snakemake@output[['discovery_gsea']]

validation_go_output <- snakemake@output[['validation_go']]
validation_kegg_output <- snakemake@output[['validation_kegg']]
validation_others_output <- snakemake@output[['validation_others']]
validation_gsea_output <- snakemake@output[['validation_gsea']]


discovery_others <- rbind(
    enrichment_rds$discovery$up_ora$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='WikiPathway'),
    enrichment_rds$discovery$up_ora$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Disease Ontology'),
    enrichment_rds$discovery$up_ora$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Network of Cancer Gene'),
    enrichment_rds$discovery$up_ora$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Disease Gene Network'),
    enrichment_rds$discovery$down_ora$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='WikiPathway'),
    enrichment_rds$discovery$down_ora$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Disease Ontology'),
    enrichment_rds$discovery$down_ora$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Network of Cancer Gene'),
    enrichment_rds$discovery$down_ora$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Disease Gene Network')
)

discovery_kegg <- rbind(
    enrichment_rds$discovery$up_ora$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='KEGG'),
    enrichment_rds$discovery$down_ora$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='KEGG')
)

discovery_go <- rbind(
    enrichment_rds$discovery$up_ora$go_bp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO BP'),
    enrichment_rds$discovery$up_ora$go_mf %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO MF'),
    enrichment_rds$discovery$up_ora$go_cc %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO CC'),
    enrichment_rds$discovery$down_ora$go_bp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO BP'),
    enrichment_rds$discovery$down_ora$go_mf %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO MF'),
    enrichment_rds$discovery$down_ora$go_cc %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO CC')
)

validation_others <- rbind(
    enrichment_rds$validation$up_ora$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='WikiPathway'),
    enrichment_rds$validation$up_ora$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Disease Ontology'),
    enrichment_rds$validation$up_ora$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Network of Cancer Gene'),
    enrichment_rds$validation$up_ora$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='Disease Gene Network'),
    enrichment_rds$validation$down_ora$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='WikiPathway'),
    enrichment_rds$validation$down_ora$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Disease Ontology'),
    enrichment_rds$validation$down_ora$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Network of Cancer Gene'),
    enrichment_rds$validation$down_ora$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='Disease Gene Network')
)

validation_kegg <- rbind(
    enrichment_rds$validation$up_ora$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='KEGG'),
    enrichment_rds$validation$down_ora$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='KEGG')
)

validation_go <- rbind(
    enrichment_rds$validation$up_ora$go_bp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO BP'),
    enrichment_rds$validation$up_ora$go_mf %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO MF'),
    enrichment_rds$validation$up_ora$go_cc %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Up-regulated',Database='GO CC'),
    enrichment_rds$validation$down_ora$go_bp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO BP'),
    enrichment_rds$validation$down_ora$go_mf %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO MF'),
    enrichment_rds$validation$down_ora$go_cc %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='Down-regulated',Database='GO CC')
)


discovery_gsea <- rbind(
    enrichment_rds$discovery$gsea$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='KEGG'),
    enrichment_rds$discovery$gsea$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='WikiPathway'),
    enrichment_rds$discovery$gsea$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='Disease Ontology'),
    enrichment_rds$discovery$gsea$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='Network of Cancer Gene'),
    enrichment_rds$discovery$gsea$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='Disease Gene Network')
)

validation_gsea <- rbind(
    enrichment_rds$validation$gsea$kegg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='KEGG'),
    enrichment_rds$validation$gsea$wp %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='WikiPathway'),
    enrichment_rds$validation$gsea$do %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='Disease Ontology'),
    enrichment_rds$validation$gsea$ncg %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='Network of Cancer Gene'),
    enrichment_rds$validation$gsea$dgn %>% as.data.frame() %>% dplyr::mutate(TargetGeneSet='GSEA',Database='Disease Gene Network')
)





write.table(validation_go,validation_go_output,quote=FALSE,sep='\t',row.names=F)
write.table(validation_kegg,validation_kegg_output,quote=FALSE,sep='\t',row.names=F)
write.table(validation_others,validation_others_output,quote=FALSE,sep='\t',row.names=F)
write.table(validation_gsea,validation_gsea_output,quote=FALSE,sep='\t',row.names=F)

write.table(discovery_go,discovery_go_output,quote=FALSE,sep='\t',row.names=F)
write.table(discovery_kegg,discovery_kegg_output,quote=FALSE,sep='\t',row.names=F)
write.table(discovery_others,discovery_others_output,quote=FALSE,sep='\t',row.names=F)
write.table(discovery_gsea,discovery_gsea_output,quote=FALSE,sep='\t',row.names=F)