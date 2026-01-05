log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


suppressPackageStartupMessages({
    library(dplyr)
    library(DESeq2)
})


parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

cat("Reading parameters and inputs...\n")
project<- snakemake@params[["project"]]
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
discovery_sample_type<-snakemake@params[["discovery_sample_type"]]
samples<-snakemake@params[["samples"]]
counts <- snakemake@input[["counts"]]


# output
cat("Preparing outputs...\n")
discovery_count_rds<-snakemake@output[["discovery_count_rds"]]
validation_count_rds<-snakemake@output[["validation_count_rds"]]
discovery_vst_rds<-snakemake@output[["discovery_vst_rds"]]
validation_vst_rds<-snakemake@output[["validation_vst_rds"]]
discovery_deg_rds<-snakemake@output[["discovery_deg_rds"]]
validation_deg_rds<-snakemake@output[["validation_deg_rds"]]
discovery_deg_tsv<-snakemake@output[["discovery_deg_tsv"]]
validation_deg_tsv<-snakemake@output[["validation_deg_tsv"]]


cat("Preparing coldata...\n")
coldata <- read.table(samples, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t',)
coldata_discovery <- coldata %>% 
                    dplyr::filter(sample_type==discovery_sample_type)
coldata_validation <- coldata %>% 
                    dplyr::filter(sample_type != discovery_sample_type)

cts <- read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')

deseq2_pipeline <- function(count,coldata,
                            condition,
                            case_condition, 
                            control_condition, parallel=TRUE){
    cts <- count[,rownames(coldata)]
    design <- as.formula(paste("~",condition))
    dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=design)
    dds <- DESeq(dds, parallel=parallel)
    vsd <- vst(dds)
    res_dds <- results(dds, contrast=c("condition",case_condition,control_condition),parallel=parallel)
    return(list(dds=dds,vsd=vsd,res_dds=res_dds))
}

save_list <- function(deseq2_list,
                          dds_file,
                          vsd_file,
                          res_file,
                          tsv_file){
    saveRDS(deseq2_list$dds, file=dds_file)
    saveRDS(deseq2_list$vsd, file=vsd_file)
    saveRDS(deseq2_list$res_dds, file=res_file)
    write.table(deseq2_list$res_dds,tsv_file,sep='\t',quote=FALSE)
}

cat("Running DESeq2 for discovery dataset...\n")
deseq2_discovery <- deseq2_pipeline(count=cts,
                                   coldata=coldata_discovery,
                                   condition="condition",
                                   case_condition=case_condition,
                                   control_condition=control_condition,
                                   parallel=parallel)
cat("Saving discovery results...\n")
save_list(deseq2_discovery,
          dds_file=discovery_count_rds,
          vsd_file=discovery_vst_rds,
          res_file=discovery_deg_rds,
          tsv_file=discovery_deg_tsv)   
cat("Running DESeq2 for validation dataset...\n")
deseq2_validation <- deseq2_pipeline(count=cts,
                                   coldata=coldata_validation,
                                   condition="condition",
                                   case_condition=case_condition,
                                   control_condition=control_condition,
                                   parallel=parallel)
cat("Saving validation results...\n")
save_list(deseq2_validation,
          dds_file=validation_count_rds,
          vsd_file=validation_vst_rds,
          res_file=validation_deg_rds,
          tsv_file=validation_deg_tsv)
cat("DESeq2 analysis completed successfully.\n")