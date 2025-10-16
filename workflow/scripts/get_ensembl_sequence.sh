#!/usr/bin/env bash

#
species=$(echo "${snakemake_params[species]}" | tr '[:upper:]' '[:lower:]')
release=${snakemake_params[release]}
build=${snakemake_params[build]}
datatype=${snakemake_params[datatype]:-""}
threads=${snakemake[threads]}
output_file=${snakemake_output[0]}
log=${snakemake_log[0]}


branch=""
if [ "$release" -ge 81 ] && [ "$build" = "GRCh37" ]; then
    branch="grch37/"
    elif [ -n "${snakemake_params[branch]}" ]; then
    branch="${snakemake_params[branch]}/"
fi


if [ "$release" -gt 75 ]; then
    spec="$build"
else
    spec="${build}.${release}"
fi


if [[ "$output_file" == *.gz ]]; then
    decompress=""
else
    decompress="gzip -dc"
fi


url="https://ftp.ensembl.org/pub"
url_prefix="${url}/${branch}release-${release}/fasta/${species}/${datatype}/${species^}.${spec}"


tmpdir=$(mktemp -d -p ./)
trap 'rm -rf "$tmpdir"' EXIT




suffixes=("dna.primary_assembly.fa.gz" "dna.toplevel.fa.gz")


fail=true
for suffix in "${suffixes[@]}"; do
    url_https="${url_prefix}.${suffix}"
    url_ftp="${url_https/https:\/\//ftp:\/\/}"
    
    if curl --location --head "$url_https" 2>/dev/null | grep -q 'Content-Length'; then
        (lftp -c "pget -n ${threads} ${url_https}" ) 2>&1 | tee -a "$log"
        ($decompress ${species^}.${spec}.${suffix} > "$output_file") 2>&1 | tee -a "$log"
        fail=false
        
        elif curl --location --head "$url_ftp" 2>/dev/null | grep -q 'Content-Length'; then
        (lftp -c "pget -n ${threads} ${url_ftp} ")  2>&1 | tee -a "$log"
        ($decompress ${species^}.${spec}.${suffix} > "$output_file") 2>&1 | tee -a "$log"
        fail=false
    fi
done



echo "Checking URL" >&2
echo "Species" >&2

if $fail; then
    exit 1
fi
