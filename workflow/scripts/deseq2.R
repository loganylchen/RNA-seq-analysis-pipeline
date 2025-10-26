log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}





discovery_dds<- readRDS(snakemake@input[['discovery_count_rds']])
validation_dds<- readRDS(snakemake@input[['validation_count_rds']])

case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]





# output
discovery_deg_rds<-snakemake@output[["discovery_deg_rds"]]
validation_deg_rds<-snakemake@output[["validation_deg_rds"]]
discovery_deg_tsv<-snakemake@output[["discovery_deg_tsv"]]
validation_deg_tsv<-snakemake@output[["validation_deg_tsv"]]


dds_discovery <- DESeq(discovery_dds)
dds_validation <- DESeq(validation_dds)



res_dds_discovery <- results(dds_discovery, contrast=c("condition",case_condition,control_condition),parallel=parallel)
res_dds_validation <- results(dds_validation, contrast=c("condition",case_condition,control_condition),parallel=parallel)


saveRDS(res_dds_discovery,file=discovery_deg_rds)
saveRDS(res_dds_validation,file=validation_deg_rds)

write.table(res_dds_discovery,discovery_deg_tsv,sep='\t',quote=FALSE)
write.table(res_dds_validation,validation_deg_tsv,sep='\t',quote=FALSE)