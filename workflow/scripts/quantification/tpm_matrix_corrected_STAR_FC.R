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
    cat("First few rows of gene length data:\n")
    print(head(gene_length_df))

    # Step 2: Read the batch-corrected count matrix
    cat("\nReading batch-corrected count matrix:", corrected_count_matrix, "\n")
    corrected_counts <- read.table(corrected_count_matrix, header=TRUE, row.names=1, check.names=FALSE, sep='\t')
    cat("Corrected count matrix dimensions:", nrow(corrected_counts), "genes x", ncol(corrected_counts), "samples.\n")
    cat("Sample names in count matrix:", paste(colnames(corrected_counts), collapse=", "), "\n")
    cat("First few gene names:", paste(head(rownames(corrected_counts)), collapse=", "), "\n")

    # Check if genes in count matrix match genes in featureCounts
    common_genes <- intersect(rownames(corrected_counts), gene_length_df$Geneid)
    cat("Common genes between count matrix and featureCounts:", length(common_genes), "\n")

    if (length(common_genes) == 0) {
        stop("ERROR: No common genes found between count matrix and featureCounts files!")
    }

    if (length(common_genes) < nrow(corrected_counts)) {
        cat("WARNING: Only", length(common_genes), "of", nrow(corrected_counts),
            "genes in count matrix found in featureCounts.\n")
        cat("Filtering count matrix to common genes...\n")
        corrected_counts <- corrected_counts[common_genes, ]
    }

    # Step 3: Calculate TPM from corrected counts
    cat("\nCalculating TPM from batch-corrected counts...\n")

    # Create a copy for TPM calculation
    tpm_df <- as.data.frame(matrix(0, nrow=nrow(corrected_counts), ncol=ncol(corrected_counts)))
    rownames(tpm_df) <- rownames(corrected_counts)
    colnames(tpm_df) <- colnames(corrected_counts)

    # Calculate TPM for each sample
    for(sample_name in colnames(corrected_counts)){
        cat("\nProcessing sample:", sample_name, "\n")

        # Get counts for this sample as a named vector
        counts <- corrected_counts[, sample_name]
        cat("  Non-zero counts:", sum(counts != 0), "/", length(counts), "\n")

        # Create data frame with gene IDs and counts
        sample_df <- data.frame(
            Geneid = rownames(corrected_counts),
            Count = as.numeric(counts),
            stringsAsFactors = FALSE
        )

        cat("  Sample data frame dimensions:", nrow(sample_df), "x", ncol(sample_df), "\n")

        # Merge with gene lengths
        sample_df <- merge(sample_df, gene_length_df, by="Geneid", all.x=TRUE)
        cat("  After merge with gene lengths:", nrow(sample_df), "genes\n")

        # Check if we have any genes with valid lengths
        if (nrow(sample_df) == 0) {
            cat("  WARNING: No genes found for sample", sample_name, "- setting TPM to 0\n")
            next
        }

        # Calculate RPK (reads per kilobase)
        sample_df <- sample_df %>%
            dplyr::mutate(length_kb = Length / 1000) %>%
            dplyr::mutate(rpk = Count / length_kb)

        # Calculate TPM: RPK / (sum(RPK) / 1e6)
        total_rpk <- sum(sample_df$rpk, na.rm=TRUE)
        cat("  Total RPK:", total_rpk, "\n")

        if (total_rpk == 0) {
            cat("  WARNING: Total RPK is 0 for sample", sample_name, "- setting TPM to 0\n")
            next
        }

        sample_df <- sample_df %>%
            dplyr::mutate(tpm = rpk / (total_rpk / 1e6))

        # Store TPM values
        tpm_df[sample_df$Geneid, sample_name] <- sample_df$tpm
        cat("  TPM range:", min(sample_df$tpm, na.rm=TRUE), "-", max(sample_df$tpm, na.rm=TRUE), "\n")
    }

    cat("\nTPM calculation completed.\n")
    cat("TPM matrix dimensions:", nrow(tpm_df), "genes x", ncol(tpm_df), "samples.\n")

    # Summary statistics
    cat("\nTPM summary:\n")
    for (sample_name in colnames(tpm_df)) {
        cat("  ", sample_name, ": ", sum(tpm_df[, sample_name], na.rm=TRUE),
            " total TPM,", sum(tpm_df[, sample_name] > 0, na.rm=TRUE), " genes detected\n")
    }

    # Step 4: Write TPM matrix
    cat("\nWriting TPM matrix to:", tpm_matrix, "\n")
    write.table(tpm_df, tpm_matrix, quote=FALSE, sep='\t', col.names=NA)
    cat("Done!\n")
}

making_TPM_from_corrected_counts(
    fc_count_files = unname(unlist(snakemake@input[["featurecounts"]])),
    corrected_count_matrix = snakemake@input[["corrected_counts"]],
    tpm_matrix = snakemake@output[['tpm_matrix']]
)

cat("\n==============================================================\n")
cat("TPM calculation complete!\n")
cat("==============================================================\n")

sink()
sink(type="message")