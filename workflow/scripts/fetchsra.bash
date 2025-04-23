#!/usr/bin/env bash

# Get command line arguments from Snakemake
log=${snakemake_log[0]}
threads=${snakemake[threads]:6}
mem_mb=${snakemake_resources[mem_mb]:-2048}
accession=${snakemake_wildcards[accession]}
outdir=$(dirname "${snakemake_output[fq1]}")

# Create temporary directory
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Ensure output directory exists
mkdir -p "$outdir"

# Run fasterq-dump with error handling
(fasterq-dump --skip-technical \
    --temp $tmpdir \
    --threads $threads \
    --mem ${mem_mb}M \
    --outdir $outdir \
    $accession

# Compress output files if they exist
for file in "$outdir"/*.fastq; do
    if [ -f "$file" ]; then
        pigz -p "$threads" "$file"
    fi
done) 2>&1 | tee "$log"