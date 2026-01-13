#!/usr/bin/env Rscript

# Extract Strandness Information from QualiMap RNA-seq Results
# Reads strandness data from qualimap rnaseq_qc_results.txt files
# and adds it as a new column to samples.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# Get parameters from Snakemake
qualimap_files <- unlist(snakemake@input[["qualimap"]])
samples_file <- snakemake@input[["samples"]]
output_file <- snakemake@output[["samples_with_strandness"]]
project <- snakemake@params[["project"]]

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("=== Extract Strandness Information ===\n")
cat("Project:", project, "\n")
cat("Samples file:", samples_file, "\n")
cat("Output file:", output_file, "\n")
cat("Number of QualiMap files:", length(qualimap_files), "\n\n")

# ============================================================================
# Function to extract strandness from QualiMap RNA-seq results
# ============================================================================
extract_strandness <- function(qualimap_file) {
  # Extract sample name from file path
  # Expected path format: {project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt
  sample_name <- basename(dirname(qualimap_file))

  # Read the QualiMap results file
  # The file has sections and we need to find the strandness section
  lines <- readLines(qualimap_file, warn=FALSE)

  # Initialize strandness values
  strandness <- NA
  fwd_strandness <- NA
  rev_strandness <- NA

  # Look for SSP estimation format: "SSP estimation (fwd/rev) = 0.9 / 0.1"
  # This is the most reliable format for QualiMap RNA-seq results
  for (i in seq_along(lines)) {
    line <- lines[i]

    # Look for SSP estimation line
    if (grepl("SSP estimation", line, ignore.case=TRUE)) {
      # Extract the forward and reverse values
      # Format: "SSP estimation (fwd/rev) = 0.9 / 0.1"
      # or: "SSP estimation (fwd/rev) = 0.90 / 0.10"

      # Extract both numbers from the line
      numbers <- str_extract_all(line, "\\d+\\.?\\d*")[[1]]

      if (length(numbers) >= 2) {
        fwd_strandness <- as.numeric(numbers[1])
        rev_strandness <- as.numeric(numbers[2])

        # Strandness is the larger of the two (dominant strand)
        strandness <- max(fwd_strandness, rev_strandness)

        cat("    Found SSP estimation: fwd =", fwd_strandness, ", rev =", rev_strandness, "\n")
        break
      }
    }
  }

  # If SSP not found, try other formats
  if (is.na(strandness)) {
    # Look for strand specificity percentage
    for (line in lines) {
      # Format: "Strand specificity : 99.5%"
      if (grepl("Strand.*specificity", line, ignore.case=TRUE)) {
        # Extract numeric value
        value <- str_extract(line, "\\d+\\.?\\d*")
        if (!is.na(value)) {
          strandness <- as.numeric(value)
          cat("    Found strand specificity:", strandness, "%\n")
          break
        }
      }
    }
  }

  return(list(
    sample_name = sample_name,
    strandness = strandness,
    fwd_strandness = fwd_strandness,
    rev_strandness = rev_strandness
  ))
}

# ============================================================================
# Read samples file
# ============================================================================
cat("Reading samples file...\n")
samples_df <- read_tsv(samples_file, show_col_types = FALSE, comment = "#")
cat("  Total samples:", nrow(samples_df), "\n")
cat("  Columns:", paste(colnames(samples_df), collapse = ", "), "\n\n")

# Filter to project samples
original_samples_df <- samples_df
samples_df <- samples_df[samples_df$project_id == project, ]
cat("  Project samples:", nrow(samples_df), "\n\n")

# ============================================================================
# Extract strandness from all QualiMap files
# ============================================================================
cat("Extracting strandness from QualiMap files...\n")
strandness_list <- list()
fwd_list <- list()
rev_list <- list()

for (qualimap_file in qualimap_files) {
  result <- extract_strandness(qualimap_file)
  strandness_list[[result$sample_name]] <- result$strandness
  fwd_list[[result$sample_name]] <- result$fwd_strandness
  rev_list[[result$sample_name]] <- result$rev_strandness

  if (!is.na(result$strandness)) {
    cat("  ", result$sample_name, ": ", result$strandness, "% ",
        "(fwd: ", result$fwd_strandness, ", rev: ", result$rev_strandness, ")\n", sep = "")
  } else {
    cat("  ", result$sample_name, ": NA (strandness not found)\n", sep = "")
  }
}

cat("\nExtracted strandness for", length(strandness_list), "samples\n")

# ============================================================================
# Add strandness column to samples dataframe
# ============================================================================
cat("\nAdding strandness column to samples dataframe...\n")

# Create strandness vector matching samples_df order
strandness_values <- sapply(samples_df$sample_name, function(x) {
  if (x %in% names(strandness_list)) {
    return(strandness_list[[x]])
  } else {
    return(NA)
  }
})

# Create forward strandness vector
fwd_strandness_values <- sapply(samples_df$sample_name, function(x) {
  if (x %in% names(fwd_list)) {
    return(fwd_list[[x]])
  } else {
    return(NA)
  }
})

# Create reverse strandness vector
rev_strandness_values <- sapply(samples_df$sample_name, function(x) {
  if (x %in% names(rev_list)) {
    return(rev_list[[x]])
  } else {
    return(NA)
  }
})

# Add the strandness columns
samples_df$strandness <- strandness_values
samples_df$strandness_fwd <- fwd_strandness_values
samples_df$strandness_rev <- rev_strandness_values

# Classify strandness based on forward and reverse values
# FWD: forward > 0.9
# REV: reverse > 0.9
# UNSTRAND: everything else
samples_df$strandness_classification <- sapply(seq_len(nrow(samples_df)), function(i) {
  fwd <- fwd_strandness_values[i]
  rev <- rev_strandness_values[i]

  if (is.na(fwd) || is.na(rev)) {
    return(NA)
  } else if (fwd > 0.9) {
    return("FWD")
  } else if (rev > 0.9) {
    return("REV")
  } else {
    return("UNSTRAND")
  }
})

# Count NA values
na_count <- sum(is.na(samples_df$strandness))
cat("  Added strandness column to", nrow(samples_df), "samples\n")
cat("  NA values:", na_count, "\n")
cat("  Non-NA values:", sum(!is.na(samples_df$strandness)), "\n")

# Print classification distribution
cat("\n  Strandness classification:\n")
class_counts <- table(samples_df$strandness_classification, useNA = "ifany")
for (class_name in names(class_counts)) {
  cat("    ", class_name, ": ", as.character(class_counts[class_name]), " samples\n", sep = "")
}
cat("\n")

# ============================================================================
# Determine overall strandness for the dataset
# ============================================================================
cat("Strandness statistics:\n")
valid_strandness <- samples_df$strandness[!is.na(samples_df$strandness)]
if (length(valid_strandness) > 0) {
  cat("  Mean:", round(mean(valid_strandness), 2), "%\n")
  cat("  Median:", round(median(valid_strandness), 2), "%\n")
  cat("  Min:", round(min(valid_strandness), 2), "%\n")
  cat("  Max:", round(max(valid_strandness), 2), "%\n")

  # Determine if stranded or unstranded
  # Threshold: if mean strandness < 20%, consider unstranded
  mean_strandness <- mean(valid_strandness)
  if (mean_strandness < 20) {
    cat("\n  Overall: Unstranded (strandness < 20%)\n")
  } else if (mean_strandness < 50) {
    cat("\n  Overall: Weakly stranded\n")
  } else if (mean_strandness < 80) {
    cat("\n  Overall: Moderately stranded\n")
  } else {
    cat("\n  Overall: Highly stranded (strandness > 80%)\n")
  }
} else {
  cat("  No valid strandness data found\n")
}

# ============================================================================
# Write output file
# ============================================================================
cat("\nWriting samples with strandness to:", output_file, "\n")
write_tsv(samples_df, output_file)
cat("Done!\n")

# ============================================================================
# Update original samples file with strandness for all samples
# ============================================================================
# For samples not in this project, set strandness to NA
all_samples_strandness <- sapply(original_samples_df$sample_name, function(x) {
  if (x %in% samples_df$sample_name) {
    return(samples_df$strandness[samples_df$sample_name == x])
  } else {
    return(NA)
  }
})

all_samples_fwd <- sapply(original_samples_df$sample_name, function(x) {
  if (x %in% samples_df$sample_name) {
    return(samples_df$strandness_fwd[samples_df$sample_name == x])
  } else {
    return(NA)
  }
})

all_samples_rev <- sapply(original_samples_df$sample_name, function(x) {
  if (x %in% samples_df$sample_name) {
    return(samples_df$strandness_rev[samples_df$sample_name == x])
  } else {
    return(NA)
  }
})

all_samples_classification <- sapply(original_samples_df$sample_name, function(x) {
  if (x %in% samples_df$sample_name) {
    return(samples_df$strandness_classification[samples_df$sample_name == x])
  } else {
    return(NA)
  }
})

original_samples_df$strandness <- all_samples_strandness
original_samples_df$strandness_fwd <- all_samples_fwd
original_samples_df$strandness_rev <- all_samples_rev
original_samples_df$strandness_classification <- all_samples_classification

# Reorder columns: put strandness columns after other metadata columns
col_order <- c(
  setdiff(colnames(original_samples_df), c("strandness", "strandness_fwd", "strandness_rev", "strandness_classification")),
  "strandness",
  "strandness_fwd",
  "strandness_rev",
  "strandness_classification"
)
original_samples_df <- original_samples_df[, col_order]

# Write the complete file (overwrite original with strandness added)
write_tsv(original_samples_df, samples_file)
cat("\nUpdated original samples file with strandness columns\n")

cat("\n=== Strandness Extraction Complete ===\n")

sink()
sink(type="message")
