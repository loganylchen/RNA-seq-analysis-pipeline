#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(GenomicRanges)
})

# Get parameters from Snakemake
merged_file <- snakemake@input[["merged"]]
gtf_file <- snakemake@input[["gtf"]]

output_annotated <- snakemake@output[["annotated"]]
output_summary <- snakemake@output[["summary"]]

cat("=== Annotating ModTect Results ===\n")

# Read merged modTect data
cat("--- Loading merged modTect data ---\n")
modtect_data <- read_tsv(merged_file, show_col_types = FALSE)

cat(sprintf("Loaded %d modification sites\n", nrow(modtect_data)))

# Extract chrom, position, reference_nt
# The merged data has chrom, position, reference_nt as first 3 columns
# Then sample columns with ModTect scores
annotated_data <- modtect_data %>%
  mutate(
    # Create unique site identifier
    site_id = paste(chrom, position, reference_nt, sep = ":"),
    # Extract genomic features
    chromosome = chrom,
    start_pos = as.numeric(position),
    end_pos = as.numeric(position)
  )

# Read GTF file for gene annotations
cat("--- Loading GTF annotation ---\n")

# Read GTF and extract all feature information
gtf_data <- read_tsv(gtf_file, show_col_types = FALSE, comment = "#")

# GTF format: seqname source feature start end score strand frame attributes
# We need to parse this manually
colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end",
                        "score", "strand", "frame", "attributes")

# Extract gene_id from all features
gtf_data <- gtf_data %>%
  mutate(
    gene_id = str_extract(attributes, 'gene_id "([^"]+)"'),
    gene_name = str_extract(attributes, 'gene_name "([^"]+)"'),
    transcript_id = str_extract(attributes, 'transcript_id "([^"]+)"')
  )

# Create separate GRanges for different features
# 1. Gene features (for gene-level annotation)
gene_annotations <- gtf_data %>%
  filter(feature == "gene")

cat(sprintf("Loaded %d gene annotations\n", nrow(gene_annotations)))

# 2. CDS features
cds_annotations <- gtf_data %>%
  filter(feature == "CDS")

cat(sprintf("Loaded %d CDS annotations\n", nrow(cds_annotations)))

# 3. UTR features
utr_annotations <- gtf_data %>%
  filter(feature %in% c("five_prime_utr", "three_prime_utr", "UTR"))

cat(sprintf("Loaded %d UTR annotations\n", nrow(utr_annotations)))

# 4. Start/stop codon features
codon_annotations <- gtf_data %>%
  filter(feature %in% c("start_codon", "stop_codon"))

cat(sprintf("Loaded %d codon annotations\n", nrow(codon_annotations)))

# Create GRanges objects for overlap analysis
modtect_gr <- makeGRangesFromDataFrame(
  annotated_data,
  keep.extra.columns = TRUE,
  seqnames.field = "chromosome",
  start.field = "start_pos",
  end.field = "end_pos"
)

gene_gr <- makeGRangesFromDataFrame(
  gene_annotations,
  seqnames.field = "seqname",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

cds_gr <- makeGRangesFromDataFrame(
  cds_annotations,
  seqnames.field = "seqname",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

utr_gr <- makeGRangesFromDataFrame(
  utr_annotations,
  seqnames.field = "seqname",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

codon_gr <- makeGRangesFromDataFrame(
  codon_annotations,
  seqnames.field = "seqname",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

# Find overlaps with genes (for gene_id and gene_symbol)
cat("--- Finding gene overlaps ---\n")
gene_overlaps <- findOverlaps(modtect_gr, gene_gr)

# Find overlaps with specific genomic features
cds_overlaps <- findOverlaps(modtect_gr, cds_gr)
utr_overlaps <- findOverlaps(modtect_gr, utr_gr)
codon_overlaps <- findOverlaps(modtect_gr, codon_gr)

# Function to get comma-separated gene IDs
get_gene_ids <- function(site_idx, overlaps, annotations) {
  matching <- which(queryHits(overlaps) == site_idx)
  if (length(matching) == 0) {
    return(".")
  }
  gene_ids <- unique(annotations$gene_id[subjectHits(overlaps)[matching]])
  paste(gene_ids, collapse = ",")
}

# Function to get comma-separated gene symbols
get_gene_symbols <- function(site_idx, overlaps, annotations) {
  matching <- which(queryHits(overlaps) == site_idx)
  if (length(matching) == 0) {
    return(".")
  }
  gene_names <- annotations$gene_name[subjectHits(overlaps)[matching]]
  # Use gene_id if gene_name is NA
  gene_names <- ifelse(is.na(gene_names),
                       annotations$gene_id[subjectHits(overlaps)[matching]],
                       gene_names)
  gene_names <- unique(gene_names)
  paste(gene_names, collapse = ",")
}

# Function to determine genomic region
get_genomic_region <- function(site_idx, cds_overlaps, utr_overlaps, codon_overlaps,
                               cds_ann, utr_ann, codon_ann) {
  # Check codon overlaps first (highest priority)
  matching <- which(queryHits(codon_overlaps) == site_idx)
  if (length(matching) > 0) {
    features <- unique(codon_ann$feature[subjectHits(codon_overlaps)[matching]])
    if ("start_codon" %in% features) {
      return("start_codon")
    } else if ("stop_codon" %in% features) {
      return("stop_codon")
    }
  }

  # Check UTR overlaps
  matching <- which(queryHits(utr_overlaps) == site_idx)
  if (length(matching) > 0) {
    features <- unique(utr_ann$feature[subjectHits(utr_overlaps)[matching]])
    # Check if we have specific UTR type or generic UTR
    if ("five_prime_utr" %in% features) {
      return("5prime_UTR")
    } else if ("three_prime_utr" %in% features) {
      return("3prime_UTR")
    } else if ("UTR" %in% features) {
      # Generic UTR - need to determine if 5' or 3' based on strand
      # For now, just mark as UTR
      return("UTR")
    }
  }

  # Check CDS overlaps
  matching <- which(queryHits(cds_overlaps) == site_idx)
  if (length(matching) > 0) {
    return("CDS")
  }

  # If none of the above, check if within gene body
  return("intergenic")
}

# Annotate each modification site
cat("--- Annotating modification sites ---\n")

n_sites <- nrow(annotated_data)
gene_ids <- character(n_sites)
gene_symbols <- character(n_sites)
genomic_regions <- character(n_sites)

for (i in 1:n_sites) {
  if (i %% 10000 == 0) {
    cat(sprintf("Processed %d/%d sites...\n", i, n_sites))
  }

  gene_ids[i] <- get_gene_ids(i, gene_overlaps, gene_annotations)
  gene_symbols[i] <- get_gene_symbols(i, gene_overlaps, gene_annotations)
  genomic_regions[i] <- get_genomic_region(i, cds_overlaps, utr_overlaps, codon_overlaps,
                                           cds_annotations, utr_annotations, codon_annotations)
}

# Create annotated data frame
final_annotated <- annotated_data %>%
  mutate(
    gene_id = gene_ids,
    gene_symbol = gene_symbols,
    gene_region = genomic_regions
  )

# Reorder columns for better readability
sample_cols <- colnames(modtect_data)[!colnames(modtect_data) %in% c("chrom", "position", "reference_nt")]

output_cols <- c(
  "chrom", "position", "reference_nt",
  "gene_id", "gene_symbol", "gene_region",
  sample_cols
)

final_annotated <- final_annotated %>%
  select(all_of(output_cols))

# Write annotated results
write_tsv(final_annotated, output_annotated)

cat(sprintf("Annotated results written to: %s\n", output_annotated))

# Generate summary statistics
n_with_genes <- sum(final_annotated$gene_id != ".", na.rm = TRUE)
n_intergenic <- sum(final_annotated$gene_region == "intergenic", na.rm = TRUE)

# Count by genomic region
region_counts <- final_annotated %>%
  filter(gene_region != "intergenic") %>%
  group_by(gene_region) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

# Write summary
summary_text <- c(
  "=== ModTect Annotation Summary ===",
  "",
  sprintf("Analysis Date: %s", Sys.time()),
  "",
  "--- Sites ---",
  sprintf("Total modification sites: %d", nrow(final_annotated)),
  sprintf("Sites overlapping genes: %d (%.1f%%)",
          n_with_genes,
          100 * n_with_genes / nrow(final_annotated)),
  sprintf("Intergenic sites: %d (%.1f%%)",
          n_intergenic,
          100 * n_intergenic / nrow(final_annotated)),
  "",
  "--- Samples ---",
  sprintf("Number of samples: %d", length(sample_cols)),
  "",
  "--- Genomic Region Distribution ---",
  if (nrow(region_counts) > 0) {
    apply(region_counts, 1, function(row) {
      sprintf("  %s: %d", row["gene_region"], row["n"])
    })
  } else {
    "  No gene-overlapping sites found"
  },
  "",
  "--- Top Genes with Most Modifications ---"
)

# Get top genes by modification count
if (n_with_genes > 0) {
  top_genes <- final_annotated %>%
    filter(gene_id != ".") %>%
    separate_rows(gene_id, sep = ",") %>%
    separate_rows(gene_symbol, sep = ",") %>%
    group_by(gene_symbol) %>%
    summarise(n_sites = n(), .groups = "drop") %>%
    arrange(desc(n_sites)) %>%
    head(20) %>%
    transmute(
      rank = row_number(),
      gene = gene_symbol,
      n_modifications = n_sites
    )

  for (i in 1:nrow(top_genes)) {
    summary_text <- c(summary_text,
      sprintf("  %d. %s: %d modifications",
              top_genes$rank[i],
              top_genes$gene[i],
              top_genes$n_modifications[i]))
  }
}

writeLines(summary_text, output_summary)

cat(sprintf("Summary written to: %s\n", output_summary))
cat("\n=== ModTect Annotation Complete ===\n")

sink()
