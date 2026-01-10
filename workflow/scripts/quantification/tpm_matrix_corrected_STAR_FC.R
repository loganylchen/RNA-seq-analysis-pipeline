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


making_TPM_from_corrected_counts <- function(fc_count_files, corrected_count_matrix, tpm_matrix) {
    # Step 1: Extract gene lengths from featureCounts files
    cat("Extracting gene lengths from featureCounts files...\n")
    gene_length_df <- NULL

    for(f in fc_count_files){
        message('Reading gene lengths from:', f)
        tmp_df <- read_tsv(f, comment = "#", progress = FALSE) %>%
                        dplyr::select(Geneid, Length)
        if(is.null(gene_length_df)){
            gene_length_df <- tmp_df
        } else {
            # Verify consistency of gene lengths across samples
            if(!all(gene_length_df$Geneid == tmp_df$Geneid)){
                stop("Gene IDs or order differs between samples. Please check featureCounts files.")
            }
        }
        break  # Only need to read one file to get gene lengths
    }

    cat("Gene lengths extracted for", nrow(gene_length_df), "genes.\n")

    # Step 2: Read the batch-corrected count matrix
    cat("Reading batch-corrected count matrix:", corrected_count_matrix, "\n")
    corrected_counts <- read.table(corrected_count_matrix, header=TRUE, row.names=1, check.names=FALSE, sep='\t')
    cat("Corrected count matrix dimensions:", nrow(corrected_counts), "genes x", ncol(corrected_counts), "samples.\n")

    # Step 3: Calculate TPM from corrected counts
    cat("Calculating TPM from batch-corrected counts...\n")

    # Create a copy for TPM calculation
    tpm_df <- as.data.frame(matrix(0, nrow=nrow(corrected_counts), ncol=ncol(corrected_counts)))
    rownames(tpm_df) <- rownames(corrected_counts)
    colnames(tpm_df) <- colnames(corrected_counts)

    # Calculate TPM for each sample
    for(sample_name in colnames(corrected_counts)){
        cat("Processing sample:", sample_name, "\n")

        # Get counts for this sample
        counts <- corrected_counts[, sample_name]

        # Merge with gene lengths
        sample_df <- data.frame(
            Geneid = names(counts),
            Count = as.numeric(counts)
        )

        sample_df <- merge(sample_df, gene_length_df, by="Geneid", all.x=TRUE)

        # Calculate RPK (reads per kilobase)
        sample_df <- sample_df %>%
            dplyr::mutate(length_kb = Length / 1000) %>%
            dplyr::mutate(rpk = Count / length_kb)

        # Calculate TPM: RPK / (sum(RPK) / 1e6)
        total_rpk <- sum(sample_df$rpk, na.rm=TRUE)
        sample_df <- sample_df %>%
            dplyr::mutate(tpm = rpk / (total_rpk / 1e6))

        # Store TPM values
        tpm_df[sample_df$Geneid, sample_name] <- sample_df$tpm
    }

    cat("TPM calculation completed.\n")
    cat("TPM matrix dimensions:", nrow(tpm_df), "genes x", ncol(tpm_df), "samples.\n")

    # Step 4: Write TPM matrix
    cat("Writing TPM matrix to:", tpm_matrix, "\n")
    write.table(tpm_df, tpm_matrix, quote=FALSE, sep='\t', col.names=NA)
    cat("Done!\n")
}

making_TPM_from_corrected_counts(
    fc_count_files = unname(unlist(snakemake@input[["featurecounts"]])),
    corrected_count_matrix = snakemake@input[["corrected_counts"]],
    tpm_matrix = snakemake@output[['tpm_matrix']]
)
