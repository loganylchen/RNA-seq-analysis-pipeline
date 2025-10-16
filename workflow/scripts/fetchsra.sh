#!/usr/bin/env bash
set -x
set -e
# Get command line arguments from Snakemake
log=${snakemake_log[0]}
exec > "$log" 2>&1

threads=${snakemake[threads]:-6}
sra_or_fastq=${snakemake_params[sra]}
fq1==${snakemake_params[fq1]}
fq2==${snakemake_params[fq2]}
mem_mb="-m${snakemake_resources[mem_mb]:-2048}M"




reads=(${snakemake_output[reads]})
read_1=${reads[0]}

outdir=$(dirname "${read_1}")
mkdir -p "${outdir}"

# Create temporary directory
tmpdir=$(mktemp -d -p ${outdir})
trap 'rm -rf "$tmpdir"' EXIT
echo ${tmpdir}
# Ensure output directory exists



if [[ -f $fq1 ]];
then
    ln -s $fq1 ${reads[0]}
    if [[ -f $fq2 ]];
    then
        ln -s $fq2 ${reads[1]}
    fi
else
    echo "Downloading data from SRA"
    echo `date`
    (fasterq-dump --skip-technical \
        --temp ${tmpdir} \
        --threads ${threads} \
        --mem "${mem_mb}" \
        --outdir ${outdir} \
        ${sra_or_fastq}
        
        # Compress output files if they exist
        for file in "$outdir"/*.fastq; do
            if [ -f "$file" ]; then
                pigz -p "$threads" "$file"
            fi
    done) 2>&1 | tee "$log"
