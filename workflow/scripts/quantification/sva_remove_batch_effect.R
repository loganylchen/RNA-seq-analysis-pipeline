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
# more logging
cat("Starting batch effect removal using ComBat-seq...\n")
cat(paste0("Batch variables: ", paste(batch_vars, collapse=", "), "\n"))

cat("Reading count matrix and coldata...\n")
count_matrix <-read.table(counts, header=TRUE,  check.names=FALSE,sep='\t')
cat("Count matrix dimensions: ", dim(count_matrix)[1], " genes and ", dim(count_matrix)[2], " samples.\n")
cat("First few rows of count matrix:\n")
print(head(count_matrix))
cat("Preparing batch information...\n")
coldata <- read.table(coldata_file, header=TRUE, check.names=FALSE,sep='\t')
cat("Coldata dimensions: ", dim(coldata)[1], " samples and ", dim(coldata)[2], " variables.\n")
cat("First few rows of coldata:\n")
print(head(coldata))
cat("Removing batch effects using ComBat-seq...\n")
batch_info <- as.matrix(coldata[, batch_vars, drop=FALSE])
cat("Batch information dimensions: ", dim(batch_info)[1], " samples and ", dim(batch_info)[2], " batch variables.\n")
cat("First few rows of batch information:\n")
print(head(batch_info))

if(length(batch_vars) == 0) {
    write.table(count_matrix, file=output_counts, sep="\t", quote=FALSE, col.names=NA)}

if(ncol(batch_info) == 1) {
    batch <- as.factor(batch_info[,1])
} else {
    batch <- apply(batch_info, 1, function(x) paste(x, collapse="_"))
}
adjusted <- ComBat_seq(count_matrix, batch=batch, group=NULL)
write.table(as.integer(adjusted), file=output_counts, sep="\t", quote=FALSE, col.names=NA)