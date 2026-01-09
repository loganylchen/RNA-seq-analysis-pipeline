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

