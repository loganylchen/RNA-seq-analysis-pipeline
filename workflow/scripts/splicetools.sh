#!/usr/bin/env bash
log=${snakemake_log[0]}

exec > "${log}" 2>&1




ri_jcec="${snakemake_input[ri_jcec]}"
se_jcec="${snakemake_input[se_jcec]}"
input_dir="${snakemake_input[input_dir]}"

annotation_bed12="${snakemake_input[annotation_bed12]}"
genome_fasta="${snakemake_input[genome_fasta]}"
tpm_file="${snakemake_input[tpm_file]}"

control_tpm_threshold="${snakemake_params[control_tpm_threshold]}"
case_tpm_threshold="${snakemake_params[case_tpm_threshold]}"
control_n="${snakemake_params[control_n]}"
case_n="${snakemake_params[case_n]}"
fdr="${snakemake_params[fdr]}"

output_dir="${snakemake_output[output_dir]}"

perl /opt/SpliceTools/bin/RIMedley.pl \
	-r ${ri_jcec} \
	-a ${annotation_bed12} \
	-g ${genome_fasta} \
	-e ${tpm_file} \
	-TPM ${control_tpm_threshold},${case_tpm_threshold} \
	-SN ${control_n},${case_n} \
	-f ${fdr}

perl /opt/SpliceTools/bin/SEMedley.pl \
	-s ${se_jcec} \
	-a ${annotation_bed12} \
	-g ${genome_fasta} \
	-e ${tpm_file} \
	-TPM ${control_tpm_threshold},${case_tpm_threshold} \
	-SN ${control_n},${case_n} \
	-f ${fdr}

perl /opt/SpliceTools/bin/SpliceCompare.pl \
	-i ${input_dir} \
	-o ${output_dir} \
	-f ${fdr}