#!/usr/bin/env bash

# Get command line arguments from Snakemake
log=${snakemake_log[0]}
threads=${snakemake[threads]:6}
mem_mb="-m${snakemake_resources[mem_mb]:-2048}M"
accession=${snakemake_wildcards[sample]}

reads=(${snakemake_output[reads]})
read_1=${reads[0]}

outdir=$(dirname "${read_1}")
mkdir -p "${outdir}"

# Create temporary directory
tmpdir=$(mktemp -d -p ${outdir})
trap 'rm -rf "$tmpdir"' EXIT
echo ${tmpdir}
# Ensure output directory exists


# Run fasterq-dump with error handling
(fasterq-dump --skip-technical \
    --temp ${tmpdir} \
    --threads ${threads} \
    --mem "${mem_mb}" \
    --outdir ${outdir} \
    ${accession}

# Compress output files if they exist
for file in "$outdir"/*.fastq; do
    if [ -f "$file" ]; then
        pigz -p "$threads" "$file"
    fi
done) 2>&1 | tee "$log"