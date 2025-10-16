#!/usr/bin/env bash
log=${snakemake_log[0]}

exec > "${log}" 2>&1

extra=${snakemake_params[extra]:-""}
threads=${snakemake[threads]:-6}
pe=${snakemake_params['pe']}

input_files=("${snakemake_input[sample]}")
n=${#input_files[@]}

echo ${input_files}
echo $n
echo $pe



if [ "$pe" == "True" ]; then
    reads="--in1 ${input_files[0]} --in2 ${input_files[1]}"
    echo 'PE reads'
else
    reads="--in1 ${input_files[0]}"
    echo 'SE reads'
fi

set -x
set -e
trimmed=""
if [ -n "${snakemake_output[trimmed]}" ]; then
    if [ "$pe" == "True" ]; then
        trimmed_files=("${snakemake_output[trimmed]}")
        trimmed="--out1 ${trimmed_files[0]} --out2 ${trimmed_files[1]}"
        echo 'PE reads'
    else
        trimmed="--out1 ${snakemake_output[trimmed]}"
        echo 'SE reads'
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