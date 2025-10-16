#!/usr/bin/env bash
log=${snakemake_log[0]}

exec > "${log}" 2>&1

extra=${snakemake_params[extra]:-""}
threads=${snakemake[threads]:-6}


input_files="${snakemake_input[sample]}"
n=${#input_files[@]}

echo ${input_files}
echo $n

if [ $n -ne 1 ] && [ $n -ne 2 ]; then
    echo "ERROR: input->sample has 1 or 2"
    exit 1
fi


if [ $n -eq 1 ]; then
    reads="--in1 ${input_files[0]}"
else
    reads="--in1 ${input_files[0]} --in2 ${input_files[1]}"
fi


trimmed=""
if [ -n "${snakemake_output[trimmed]}" ]; then
    if [ $n -eq 1 ]; then
        trimmed="--out1 ${snakemake_output[trimmed]}"
    else
        trimmed_files=("${snakemake_output[trimmed]}")
        trimmed="--out1 ${trimmed_files[0]} --out2 ${trimmed_files[1]}"
    fi
fi




html="--html ${snakemake_output[html]}"
json="--json ${snakemake_output[json]}"


(fastp --thread "$threads" \
    $extra \
    $reads \
    $trimmed \
    $json \
$html)