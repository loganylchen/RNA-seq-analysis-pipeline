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
    end_pos = as.numeric(position),
    strand = "."
  )

# Read GTF file for gene annotations
cat("--- Loading GTF annotation ---\n")

# Read GTF and extract gene information
gtf_data <- read_tsv(gtf_file, show_col_types = FALSE, comment = "#")

# GTF format: seqname source feature start end score strand frame attributes
# We need to parse this manually
colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end",
                        "score", "strand", "frame", "attributes")

# Filter for gene features
gene_annotations <- gtf_data %>%
  filter(feature == "gene") %>%
  mutate(
    gene_id = str_extract(attributes, 'gene_id "([^"]+)"'),
    gene_name = str_extract(attributes, 'gene_name "([^"]+)"'),
    gene_biotype = str_extract(attributes, 'gene_biotype "([^"]+)"')
  )

cat(sprintf("Loaded %d gene annotations\n", nrow(gene_annotations)))

# Create GRanges objects for overlap analysis
modtect_gr <- makeGRangesFromDataFrame(
  annotated_data,
  keep.extra.columns = TRUE,
  seqnames.field = "chromosome",
  start.field = "start_pos",
  end.field = "end_pos",
  strand.field = "strand"
)

gene_gr <- makeGRangesFromDataFrame(
  gene_annotations,
  seqnames.field = "seqname",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

# Find overlaps with genes
cat("--- Finding gene overlaps ---\n")
overlaps <- findOverlaps(modtect_gr, gene_gr)

# Map modifications to genes
modtect_genes <- data.frame(
  modtect_idx = queryHits(overlaps),
  gene_idx = subjectHits(overlaps)
) %>%
  mutate(
    site_id = annotated_data$site_id[modtect_idx],
    gene_id = gene_annotations$gene_id[gene_idx],
    gene_name = ifelse(is.na(gene_annotations$gene_name[gene_idx]),
                      gene_annotations$gene_id[gene_idx],
                      gene_annotations$gene_name[gene_idx]),
    gene_biotype = gene_annotations$gene_biotype[gene_idx],
    gene_chrom = gene_annotations$seqname[gene_idx],
    gene_start = gene_annotations$start[gene_idx],
    gene_end = gene_annotations$end[gene_idx],
    gene_strand = gene_annotations$strand[gene_idx]
  ) %>%
  group_by(site_id) %>%
  summarise(
    gene_id = paste(unique(gene_id), collapse = ";"),
    gene_name = paste(unique(gene_name), collapse = ";"),
    gene_biotype = paste(unique(gene_biotype), collapse = ";"),
    n_genes = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Determine genomic context
    # If within gene body: exonic, intronic (we'll just call it "gene_body")
    # If within 2kb upstream: promoter
    # If within 2kb downstream: downstream
    # Otherwise: intergenic
    genomic_region = case_when(
      n_genes > 0 ~ "gene_body",
      TRUE ~ "intergenic"
    )
  )

# Merge annotations back to modTect data
final_annotated <- annotated_data %>%
  left_join(modtect_genes, by = "site_id") %>%
  mutate(
    gene_name = ifelse(is.na(gene_name), ".", gene_name),
    gene_id = ifelse(is.na(gene_id), ".", gene_id),
    genomic_region = ifelse(is.na(genomic_region), "intergenic", genomic_region),
    # Create a more descriptive annotation
    annotation = case_when(
      genomic_region == "gene_body" ~ paste0(gene_name, " (", genomic_region, ")"),
      genomic_region == "intergenic" ~ "intergenic",
      TRUE ~ genomic_region
    )
  )

# Reorder columns for better readability
sample_cols <- colnames(modtect_data)[!colnames(modtect_data) %in% c("chrom", "position", "reference_nt")]

output_cols <- c(
  "chrom", "position", "reference_nt", "site_id",
  "gene_name", "gene_id", "genomic_region", "annotation", "n_genes",
  sample_cols
)

final_annotated <- final_annotated %>%
  select(all_of(output_cols))

# Write annotated results
write_tsv(final_annotated, output_annotated)

cat(sprintf("Annotated results written to: %s\n", output_annotated))

# Generate summary statistics
summary_stats <- list(
  n_total_sites = nrow(final_annotated),
  n_sites_with_genes = sum(final_annotated$n_genes > 0, na.rm = TRUE),
  n_intergenic_sites = sum(final_annotated$genomic_region == "intergenic", na.rm = TRUE),
  n_genes_detected = sum(final_annotated$n_genes, na.rm = TRUE),
  n_samples_in_data = length(sample_cols),
  chromosomes = paste(unique(final_annotated$chrom), collapse = ", ")
)

# Write summary
summary_text <- c(
  "=== ModTect Annotation Summary ===",
  "",
  sprintf("Analysis Date: %s", Sys.time()),
  "",
  "--- Sites ---",
  sprintf("Total modification sites: %d", summary_stats$n_total_sites),
  sprintf("Sites overlapping genes: %d (%.1f%%)",
          summary_stats$n_sites_with_genes,
          100 * summary_stats$n_sites_with_genes / summary_stats$n_total_sites),
  sprintf("Intergenic sites: %d (%.1f%%)",
          summary_stats$n_intergenic_sites,
          100 * summary_stats$n_intergenic_sites / summary_stats$n_total_sites),
  sprintf("Total gene associations: %d", summary_stats$n_genes_detected),
  "",
  "--- Samples ---",
  sprintf("Number of samples: %d", summary_stats$n_samples_in_data),
  "",
  "--- Genomic Distribution ---",
  sprintf("Chromosomes: %s", substr(summary_stats$chromosomes, 1, 200)),
  "",
  "--- Top Genes with Most Modifications ---"
)

# Get top genes by modification count
if (summary_stats$n_sites_with_genes > 0) {
  top_genes <- final_annotated %>%
    filter(n_genes > 0) %>%
    separate_rows(gene_name, sep = ";") %>%
    group_by(gene_name) %>%
    summarise(n_sites = n(), .groups = "drop") %>%
    arrange(desc(n_sites)) %>%
    head(20) %>%
    transmute(
      rank = row_number(),
      gene = gene_name,
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
