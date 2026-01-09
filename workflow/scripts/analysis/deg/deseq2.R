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
dataset<- snakemake@params[["dataset"]]
case_condition<-snakemake@params[["case_condition"]]
control_condition<-snakemake@params[["control_condition"]]
samples<-snakemake@params[["samples"]]
design_string<-snakemake@params[["design"]]
counts <- snakemake@input[["counts"]]


# output
cat("Preparing outputs...\n")
count_rds<-snakemake@output[["count_rds"]]
vst_rds<-snakemake@output[["vst_rds"]]
deg_rds<-snakemake@output[["deg_rds"]]
deg_tsv<-snakemake@output[["deg_tsv"]]


cat("Preparing coldata...\n")
coldata <- read.table(samples, header=TRUE, row.names=1, check.names=FALSE,sep='\t',) %>%
            dplyr::filter(dataset_id==dataset)
cts <- read.table(counts, header=TRUE, row.names=1, check.names=FALSE,sep='\t')

deseq2_pipeline <- function(design_string,count,
                            coldata,
                            condition,
                            case_condition, 
                            control_condition, parallel=TRUE){
    cts <- count[,rownames(coldata)]
    if(design_string == ""){
        design_string <- "~ condition"
    }else{
        design_string <- paste0("~",design_string,"+ condition")
    }
    design <- as.formula(design_string)
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
deseq2_data <- deseq2_pipeline(design_string=design_string,
                                count=cts,
                                   coldata=coldata,
                                   condition="condition",
                                   case_condition=case_condition,
                                   control_condition=control_condition,
                                   parallel=parallel)
cat("Saving discovery results...\n")
save_list(deseq2_data,
          dds_file=count_rds,
          vsd_file=vst_rds,
          res_file=deg_rds,
          tsv_file=deg_tsv)   
    
cat("DESeq2 analysis completed successfully.\n")