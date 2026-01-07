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


download_sra_pair() {
    threads=$1
    srr=$2
    read1=$3  # e.g., "sample1" -> becomes sample1.read1.fq.gz
    read2=$4
    
    # Get URLs from ENA
    response=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=$srr&result=read_run&fields=fastq_ftp")
    urls=$(echo "$response" | tail -n +2 | cut -f2)
    
    # Convert to array (semicolon separated)
    IFS=';' read -r -a url_array <<< "$urls"
    
    if [ ${#url_array[@]} -eq 1 ]; then
        # Single-end
        echo "Single-end data detected"
        lftp -c "pget -n ${threads} ftp://${url_array[0]} -o ${read1}"
    elif [ ${#url_array[@]} -eq 2 ]; then
        # Paired-end
        echo "Paired-end data detected"
        # Download both in parallel using background processes
        lftp -c "pget -n ${threads} ftp://${url_array[0]} -o ${read1}";
        lftp -c "pget -n ${threads} ftp://${url_array[1]} -o ${read2}";
    else
        echo "Error: Unexpected number of files: ${#url_array[@]}"
        return 1
    fi
    
    echo "Download complete!"
}





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

  

    # Download SRA file using fasterq-dump with gzip compression by default
    
    
 
    download_sra_pair "${threads}" "${sra_id}" "${read_1}" "${read_2}"
   

    

    echo "SRA data download completed"
fi

echo "=== Sample ${sample_name} processing complete ==="
