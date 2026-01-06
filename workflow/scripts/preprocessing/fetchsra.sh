#!/usr/bin/env bash

# Get command line arguments from Snakemake
log=${snakemake_log[0]}
exec > "$log" 2>&1
set -x
set -e

# Get parameters
threads=${snakemake[threads]:-6}
sra_id=${snakemake_params[sra]}
fq1=${snakemake_params[fq1]}
fq2=${snakemake_params[fq2]}
sample_name=${snakemake_params[sample]}
project=${snakemake_params[project]}

# Get output paths
reads=(${snakemake_output[reads]})
read_1=${reads[0]}
read_2=${reads[1]}
outdir=$(dirname "${read_1}")
mkdir -p "${outdir}"

echo "=== Processing sample: ${sample_name} ==="
echo "Project: ${project}"
echo "SRA ID: ${sra_id}"
echo "Output directory: ${outdir}"
echo "Read 1 output: ${read_1}"
echo "Read 2 output: ${read_2}"

# Check if fastq files are provided
if [[ -n "$fq1" && -f "$fq1" ]]; then
    echo "Fastq files provided, creating symlinks..."
    ln -sf "$(readlink -f "$fq1")" "${read_1}"
    if [[ -n "$fq2" && -f "$fq2" ]]; then
        ln -sf "$(readlink -f "$fq2")" "${read_2}"
    fi
    echo "Symlinks created successfully"
else
    # No fastq files provided, fetch from SRA
    echo "No fastq files provided, fetching from SRA..."

    if [[ -z "$sra_id" ]]; then
        echo "Error: Neither fastq files nor SRA ID provided!"
        exit 1
    fi

    echo "Fetching SRA: ${sra_id}"

    # Create temp directory for SRA download
    tmp_dir="${outdir}/sra_tmp"
    mkdir -p "${tmp_dir}"

    # Download SRA file using fasterq-dump with gzip compression by default
    cd "${tmp_dir}"
    fasterq-dump --threads ${threads} --split-files --gzip --progress ${sra_id}

    # Check if files were downloaded
    if [[ ! -f "${sra_id}_1.fastq.gz" ]]; then
        echo "Warning: No paired-end files found, checking for single-end..."
        if [[ -f "${sra_id}.fastq.gz" ]]; then
            # Single-end data
            mv "${sra_id}.fastq.gz" "${read_1}"
            echo "Single-end data downloaded and moved to ${read_1}"
        else
            echo "Error: Failed to download SRA data!"
            exit 1
        fi
    else
        # Paired-end data
        mv "${sra_id}_1.fastq.gz" "${read_1}"
        if [[ -f "${sra_id}_2.fastq.gz" ]]; then
            mv "${sra_id}_2.fastq.gz" "${read_2}"
        fi
        echo "Paired-end data downloaded and moved"
    fi

    # Clean up SRA temp directory
    cd -
    rm -rf "${tmp_dir}"

    echo "SRA data download completed"
fi

echo "=== Sample ${sample_name} processing complete ==="
