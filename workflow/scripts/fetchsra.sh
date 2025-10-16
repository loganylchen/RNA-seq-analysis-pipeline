#!/usr/bin/env bash

# Get command line arguments from Snakemake
log=${snakemake_log[0]}
exec > "$log" 2>&1
set -x
set -e
threads=${snakemake[threads]:-6}
sra_or_fastq=${snakemake_params[sra]}
fq1=${snakemake_params[fq1]}
fq2=${snakemake_params[fq2]}
mem_mb="-m${snakemake_resources[mem_mb]:-2048}M"




reads=(${snakemake_output[reads]})
read_1=${reads[0]}
read_2=${reads[1]}
outdir=$(dirname "${read_1}")
mkdir -p "${outdir}"


if [[ -f $fq1 ]];
then
    ln -s $fq1 ${read_1}
    if [[ -f $fq2 ]];
    then
        ln -s $fq2 ${read_2}
    fi
fi
