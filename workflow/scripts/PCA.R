log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


suppressPackageStartupMessages({
    library(dplyr)
    library(DESeq2)
    library(PCAtools)
})


parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}


project<- snakemake@params[["project"]]
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
discovery_sample_type<-snakemake@params[["discovery_sample_type"]]
samples<-snakemake@params[["samples"]]
counts <- snakemake@input[["counts"]]


# output
discovery_count_rds<-snakemake@output[["discovery_count_rds"]]
validation_count_rds<-snakemake@output[["validation_count_rds"]]
discovery_vst_rds<-snakemake@output[["discovery_vst_rds"]]
validation_vst_rds<-snakemake@output[["validation_vst_rds"]]


coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t')
coldata_discovery <- coldata %>% 
                    dplyr::filter(sample_type==discovery_sample_type)

coldata_validation <- coldata %>% 
                    dplyr::filter(sample_type != discovery_sample_type)

cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')
cts_discovery <- cts[,rownames(coldata_discovery)]
cts_validation <- cts[,rownames(coldata_validation)]


dds_discovery <- DESeqDataSetFromMatrix(countData=cts_discovery,
                              colData=coldata_discovery,
                              design=~condition)
dds_validation <- DESeqDataSetFromMatrix(countData=cts_validation,
                              colData=coldata_validation,
                              design=~condition)


dds_discovery <- DESeq(dds_discovery, parallel=parallel)
dds_validation <- DESeq(dds_validation, parallel=parallel)

# Write dds object as RDS
saveRDS(dds_discovery, file=discovery_count_rds)
saveRDS(dds_validation, file=validation_count_rds)

vsd_discovery <- vst(dds_discovery)
vsd_validation <- vst(dds_validation)


saveRDS(vsd_discovery, file=discovery_vst_rds)
saveRDS(vsd_validation, file=validation_vst_rds)