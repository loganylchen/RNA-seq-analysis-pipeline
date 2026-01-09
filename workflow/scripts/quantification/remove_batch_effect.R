#!/usr/bin/env Rscript
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


suppressPackageStartupMessages({
    library(dplyr)
    library(sva)
    library(readr)
})


batch_vars <- unname(unlist(snakemake@params[["batch_vars"]]))
counts <- snakemake@input[["counts"]]
coldata_file <- snakemake@input[["coldata"]]
output_counts <- snakemake@output[["counts"]]


count_matrix <-read.table(counts, header=TRUE, row.names="Geneid", check.names=FALSE,sep='\t')
coldata <- read.table(coldata_file, header=TRUE, row.names="sample_name", check.names=FALSE,sep='\t')
batch_info <- as.matrix(coldata[, batch_vars, drop=FALSE])

if(length(batch_vars) == 0) {
    write.table(count_matrix, file=output_counts, sep="\t", quote=FALSE, col.names=NA)}

if(ncol(batch_info) == 1) {
    batch <- as.factor(batch_info[,1])
} else {
    batch <- apply(batch_info, 1, function(x) paste(x, collapse="_"))
}
adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)
write.table(as.integer(adjusted), file=output_counts, sep="\t", quote=FALSE, col.names=NA)