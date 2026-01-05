#!/usr/bin/env bash


species=$(echo "${snakemake_params[species]}" | tr '[:upper:]' '[:lower:]')
build=${snakemake_params[build]}
release=${snakemake_params[release]}
gtf_release=$release
output_file=${snakemake_output[0]}
log=${snakemake_log[0]}


if [[ "$output_file" == *.gz ]]; then
    out_gz=true
    out_fmt=$(basename "$output_file" | sed 's/\.[^.]*$//' | sed 's/.*\.//')
else
    out_gz=false
    out_fmt=$(basename "$output_file" | sed 's/.*\.//')
fi


branch=""
if [ "$build" = "GRCh37" ]; then
    if [ "$release" -ge 81 ]; then
        branch="grch37/"
    fi
    if [ "$release" -gt 87 ]; then
        gtf_release=87
    fi
    elif [ -n "${snakemake_params[branch]}" ]; then
    branch="${snakemake_params[branch]}/"
fi


flavor=${snakemake_params[flavor]:-""}
if [ -n "$flavor" ]; then
    flavor="${flavor}."
fi


if [ "$out_fmt" = "gtf" ]; then
    suffix="gtf.gz"
    elif [ "$out_fmt" = "gff3" ]; then
    suffix="gff3.gz"
else
    echo "ERROR :'gtf[.gz]' and 'gff3[.gz]'。" >&2
    exit 1
fi


url=${snakemake_params[url]:-"ftp://ftp.ensembl.org/pub"}
url="${url}/${branch}release-${release}/${out_fmt}/${species}/${species^}.${build}.${gtf_release}.${flavor}${suffix}"


if $out_gz; then
    if ! curl -L "$url" > "$output_file" 2>> "$log"; then
        echo "Ensembl" >&2
        echo "" >&2
        exit 1
    fi
else
    if ! (curl -L "$url" | gzip -d > "$output_file") 2>> "$log"; then
        echo "Ensembl" >&2
        echo "" >&2
        exit 1
    fi
fi