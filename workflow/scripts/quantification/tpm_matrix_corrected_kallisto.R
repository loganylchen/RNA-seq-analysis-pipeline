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


making_TPM_from_corrected_kallisto_counts <- function(kallisto_files, corrected_count_matrix, tpm_matrix) {
    # Step 1: Extract effective lengths from Kallisto abundance.tsv files
    cat("Extracting effective lengths from Kallisto abundance.tsv files...\n")
    effective_length_list <- list()

    for(f in kallisto_files){
        sample_name <- basename(dirname(f))
        message('Reading effective lengths from:', f, 'for sample:', sample_name)
        tmp_df <- read_tsv(f, comment = "#", progress = FALSE) %>%
                        dplyr::select(target_id, eff_length) %>%
                        dplyr::mutate(Sample = sample_name)
        effective_length_list[[sample_name]] <- tmp_df
    }

    cat("Effective lengths extracted for", length(effective_length_list), "samples.\n")

    # Step 2: Read the batch-corrected count matrix
    cat("Reading batch-corrected count matrix:", corrected_count_matrix, "\n")
    corrected_counts <- read.table(corrected_count_matrix, header=TRUE, row.names=1, check.names=FALSE, sep='\t')
    cat("Corrected count matrix dimensions:", nrow(corrected_counts), "genes x", ncol(corrected_counts), "samples.\n")

    # Step 3: Calculate TPM from corrected counts using sample-specific effective lengths
    cat("Calculating TPM from batch-corrected counts using effective lengths...\n")

    # Create a data frame for TPM values
    tpm_df <- as.data.frame(matrix(0, nrow=nrow(corrected_counts), ncol=ncol(corrected_counts)))
    rownames(tpm_df) <- rownames(corrected_counts)
    colnames(tpm_df) <- colnames(corrected_counts)

    # Calculate TPM for each sample
    for(sample_name in colnames(corrected_counts)){
        cat("Processing sample:", sample_name, "\n")

        if(!sample_name %in% names(effective_length_list)){
            warning(paste("Sample", sample_name, "not found in Kallisto files. Skipping TPM calculation for this sample."))
            next
        }

        # Get counts for this sample
        counts <- corrected_counts[, sample_name]

        # Get effective lengths for this sample
        eff_length_df <- effective_length_list[[sample_name]]

        # Prepare data frame for calculation
        sample_df <- data.frame(
            target_id = names(counts),
            Count = as.numeric(counts),
            stringsAsFactors = FALSE
        )

        # Merge with effective lengths
        sample_df <- merge(sample_df, eff_length_df, by="target_id", all.x=TRUE)

        # Calculate TPM: TPM = (Count / eff_length) * 1e3 / sum(Count / eff_length) * 1e6
        # This is equivalent to: RPK = Count / (eff_length / 1000), then TPM = RPK / (sum(RPK) / 1e6)
        sample_df <- sample_df %>%
            dplyr::mutate(length_kb = eff_length / 1000) %>%
            dplyr::mutate(rpk = Count / length_kb)

        # Handle zero or negative effective lengths
        sample_df$rpk[is.na(sample_df$rpk) | is.infinite(sample_df$rpk)] <- 0

        # Calculate scaling factor
        total_rpk <- sum(sample_df$rpk, na.rm=TRUE)

        if(total_rpk > 0){
            sample_df <- sample_df %>%
                dplyr::mutate(tpm = (rpk / total_rpk) * 1e6)
        } else {
            sample_df$tpm <- 0
        }

        # Store TPM values
        rownames(sample_df) <- sample_df$target_id
        tpm_df[sample_df$target_id, sample_name] <- sample_df$tpm
    }

    cat("TPM calculation completed.\n")
    cat("TPM matrix dimensions:", nrow(tpm_df), "genes x", ncol(tpm_df), "samples.\n")

    # Step 4: Write TPM matrix
    cat("Writing TPM matrix to:", tpm_matrix, "\n")
    write.table(tpm_df, tpm_matrix, quote=FALSE, sep='\t', col.names=NA)
    cat("Done!\n")
}

making_TPM_from_corrected_kallisto_counts(
    kallisto_files = unname(unlist(snakemake@input[["kallisto"]])),
    corrected_count_matrix = snakemake@input[["corrected_counts"]],
    tpm_matrix = snakemake@output[['tpm_matrix']]
)
