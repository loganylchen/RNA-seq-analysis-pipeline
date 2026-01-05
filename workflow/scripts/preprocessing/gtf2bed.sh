#!/usr/bin/env bash
log=${snakemake_log[0]}

exec > "${log}" 2>&1


threads=${snakemake[threads]:-6}


gtf="${snakemake_input[gtf]}"



bed="${snakemake_output[bed]}"

awk '$3 == "exon" {print $1 "\t" ($4-1) "\t" $5}' ${gtf} > ${bed}.tmp;
bedtools sort -i ${bed}.tmp > ${bed}.sorted;
bedtools merge -i ${bed}.sorted  > ${bed};
rm ${bed}.tmp ${bed}.sorted;
    