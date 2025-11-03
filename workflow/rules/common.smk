import pandas as pd
import glob
import re
from snakemake.utils import validate


READ_STRAND_INFER = re.compile(r"SSP estimation \(fwd/rev\) = (\d+\.\d+) / (\d+\.\d+)")


validate(config, schema="../schemas/config.schema.yaml")
project = config["project"]
case_condition = config["case_condition"]
control_condition = config["control_condition"]
discovery_sample_type = config["discovery_sample_type"]


samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

samples = samples.loc[samples["project"] == project]
case_samples = samples.loc[samples["condition"] == case_condition]
control_samples = samples.loc[samples["condition"] == control_condition]
discovery_samples = samples.loc[samples["sample_type"] == discovery_sample_type]

# project = samples["project"].unique().tolist()
# assert len(project) == 1, "Only one project is allowed!"
# project = project[0]


validate(samples, schema="../schemas/samples.schema.yaml")


def is_pe(wildcards):
    if samples.loc[wildcards.sample].loc["seq_type"] == "pe":
        return True
    else:
        return False


def get_sequence_type(qc_file):

    # SSP estimation (fwd/rev) = 0.44 / 0.56
    with open(qc_file) as fin:
        for line in fin:
            match = READ_STRAND_INFER.search(line)
            if match:
                fwd = float(match.group(1))
                rev = float(match.group(2))
                if fwd >= 0.8:
                    return "FWD"
                elif rev >= 0.8:
                    return "REV"
                else:
                    return "UNSTRAND"

    return "UNSTRAND"


def stringtie_strand_infer(qc_file):
    strand = get_sequence_type(qc_file)
    if strand == "FWD":
        return " --fr "
    elif strand == "REV":
        return " --rf "
    else:
        return ""


def featurecounts_strand_infer(qc_file):
    strand = get_sequence_type(qc_file)
    if strand == "FWD":
        return " -s 1 "
    elif strand == "REV":
        return " -s 2 "
    else:
        return " -s 0 "


def salmon_strand_infer(qc_file):
    strand = get_sequence_type(qc_file)
    if strand == "FWD":
        return " -l ISF "
    elif strand == "REV":
        return " -l ISR "
    else:
        return " -l IU "


def kallisto_strand_infer(qc_file):
    strand = get_sequence_type(qc_file)
    if strand == "FWD":
        return " --fr-stranded "
    elif strand == "REV":
        return " --rf-stranded "
    else:
        return "  "


def get_sra(wildcards):
    return samples.loc[wildcards.sample].loc["raw_data"]


def get_fq1(wildcards):
    return samples.loc[wildcards.sample].loc["fq1"]


def get_fq2(wildcards):
    return samples.loc[wildcards.sample].loc["fq2"]


def get_raw_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc["raw_data"]
    seq_type = samples.loc[wildcards.sample].loc["seq_type"]
    if seq_type == "se":
        return {
            "fq1": f"{wildcards.project}/data/{wildcards.sample}/{wildcards.sample}.fastq.gz",
            "fq2": "",
        }
    elif seq_type == "pe":
        return {
            "fq1": f"{wildcards.project}/data/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
            "fq2": f"{wildcards.project}/data/{wildcards.sample}/{wildcards.sample}_2.fastq.gz",
        }


def get_clean_data(wildcards):
    if samples.loc[wildcards.sample].loc["seq_type"] == "pe":
        return {
            "fq1": f"{wildcards.project}/clean_data/{wildcards.sample}_1.fastq.gz",
            "fq2": f"{wildcards.project}/clean_data/{wildcards.sample}_2.fastq.gz",
            "reads": [
                f"{wildcards.project}/clean_data/{wildcards.sample}_1.fastq.gz",
                f"{wildcards.project}/clean_data/{wildcards.sample}_2.fastq.gz",
            ],
        }
    elif samples.loc[wildcards.sample].loc["seq_type"] == "se":
        return {
            "fq1": f"{wildcards.project}/clean_data/{wildcards.sample}.fastq.gz",
            "fq2": "",
            "reads": [f"{wildcards.project}/clean_data/{wildcards.sample}.fastq.gz"],
        }
    else:
        raise ValueError(f"{wildcards.sample} is a wired name!")


def get_qc_files():
    qc_files = []
    for sample in samples.index:
        sample_project = samples.loc[sample, "project"]
        qc_files += [
            f"{sample_project}/qc/fastp/{sample}/{sample}.fastp.json",
            f"{sample_project}/qc/STAR/{sample}/{sample}.Log.final.out",
            f"{sample_project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
            f"{sample_project}/quantification/salmon/{sample}/",
        ]
    return qc_files


def get_final_output():
    final_output = [
        "resources/star_genome",
        # "resources/hisat2_genome/genome.1.ht2",
        # f"{project}/quantification/STAR_fc_count_matrix.txt",
        # f"{project}/quantification/HISAT2_fc_count_matrix.txt",
    ]
    for sample in samples.index:
        sample_project = samples.loc[sample, "project"]

        final_output += [
            # f"{sample_project}/quantification/{sample}.hisat2_counts.txt",
            # f"{sample_project}/qc/{sample}/rnaseq_qc_results.txt",
            # f"{sample_project}/quantification/{sample}.star_counts.txt",
            # f"{sample_project}/alignment/{sample}/{sample}.star.bam",
            # f"{sample_project}/assembly/{sample}/{sample}.stringtie.gtf",
            # f"{sample_project}/quantification/{sample}.salmon",
            # f"{sample_project}/modification/{sample}/modtect/{sample}.modtect.combined.txt",
        ]

    final_output += [
        # f"{sample_project}/assembly/gffcompare.annotated.gtf",
        # f"{sample_project}/transcript_splicing/rmats",
        # f"{sample_project}/assembly/transcriptome.fasta",
        # f"{sample_project}/quantification/STAR_fc_count_matrix_PUREE.txt",
        # f"{sample_project}/deseq2/discovery_count_matrix.rds",
        # f"{sample_project}/deseq2/validation_deg.rds",
        f"{sample_project}/visualization/Volcano_validation.pdf",
        # f"{sample_project}/enrichment/validation_gsea_enrichment.tsv",
        f"{sample_project}/qc/multiqc/",
    ]

    return final_output
