#!/usr/bin/env bash
log=${snakemake_log[0]}

exec > "${log}" 2>&1

extra=${snakemake_params[extra]:-""}
threads=${snakemake[threads]:-6}
pe=${snakemake_params['pe']}

fq1="${snakemake_input[fq1]}"
fq2="${snakemake_input[fq2]}"

clean_fq1="${snakemake_output[fq1]}"
clean_fq2="${snakemake_output[fq2]}"






if [ "$pe" == "True" ]; then
    reads="--in1 ${fq1} --in2 ${fq2}"
    echo 'PE reads'
else
    reads="--in1 ${fq1}"
    echo 'SE reads'
fi

set -x
set -e
trimmed=""

if [ "$pe" == "True" ]; then

    trimmed="--out1 ${clean_fq1} --out2 ${clean_fq2}"
    echo 'PE reads'
else
    trimmed="--out1 ${clean_fq1}"
    echo 'SE reads'
fi





html="--html ${snakemake_output[html]}"
json="--json ${snakemake_output[json]}"


(fastp --thread "$threads" \
    $extra \
    $reads \
    $trimmed \
    $json \
$html)