import pandas as pd
import glob
import re
from snakemake.utils import validate


READ_STRAND_INFER = re.compile(r"SSP estimation \(fwd/rev\) = (\d+\.\d+) / (\d+\.\d+)")

QUANTIFICATION_TOOLS = ["STAR_FC", "salmon", "kallisto"]
validate(config, schema="../schemas/config.schema.yaml")
project = config["project"]

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str}, comment="#")
    .set_index("sample_name", drop=False)
    .sort_index()
)

samples = samples.loc[samples["project_id"] == project]
datasets = samples["dataset_id"].unique().tolist()

# validate(samples, schema="../schemas/samples.schema.yaml")


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


def get_samples_strandness(qc_files):
    strandness_dict = {}
    for qc_file in qc_files:
        sample_name = qc_file.split("/")[-2]
        strandness_dict[sample_name] = get_sequence_type(qc_file)
    return strandness_dict


def get_info(wildcards):
    project = wildcards.project
    dataset_samples = f"{project}/config/samples_with_strandness.tsv"
    return dataset_samples


def get_condition(wildcards, condition_type):
    dataset = wildcards.dataset
    if condition_type == "case":
        return config["datasets"][dataset]["case_condition"]
    elif condition_type == "control":
        return config["datasets"][dataset]["control_condition"]
    else:
        raise ValueError(
            f"condition_type should be 'case' or 'control', got {condition_type}"
        )


def get_wildcards_element(wildcards, element):
    return getattr(wildcards, element)


def get_dataset(wildcards):
    return get_wildcards_element(wildcards, "dataset")


def get_project(wildcards):
    return get_wildcards_element(wildcards, "project")


def get_case_condition(wildcards):
    return get_condition(wildcards, "case")


def get_control_condition(wildcards):
    return get_condition(wildcards, "control")


def get_dataset_samples(wildcards):
    project = wildcards.project
    dataset = wildcards.dataset
    dataset_samples = samples[
        (samples["dataset_id"] == dataset) & (samples["project_id"] == project)
    ].index.tolist()
    return dataset_samples


def get_deseq2_design(wildcards):
    dataset = wildcards.dataset
    design = config["datasets"][dataset]["deseq2_design_formula"]
    return design


def get_edeger_design(wildcards):
    dataset = wildcards.dataset
    design = config["datasets"][dataset]["edeger_design_formula"]
    return design


def get_batch_vars(wildcards):
    dataset = wildcards.dataset
    batch_vars = config["datasets"][dataset].get("batch_vars", [])
    return batch_vars


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


def hisat2_strand_infer(qc_file):
    strand = get_sequence_type(qc_file)
    if strand == "FWD":
        return " --rna-strandness FR"
    elif strand == "REV":
        return " --rna-strandness RF"
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
        sample_project = samples.loc[sample, "project_id"]
        qc_files += [
            f"{sample_project}/qc/fastp/{sample}/{sample}.fastp.json",
            f"{sample_project}/qc/STAR/{sample}/{sample}.Log.final.out",
            f"{sample_project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
            f"{sample_project}/quantification/salmon/{sample}/",
            f"{sample_project}/qc/kallisto/{sample}/kallisto.log",
            f"{sample_project}/qc/hisat2/{sample}/{sample}.log",
            f"{sample_project}/qc/picard/{sample}/{sample}.alignment_summary_metrics.txt",
            f"{sample_project}/qc/picard/{sample}/{sample}.rnaseq_metrics.txt",
            f"{sample_project}/qc/picard/{sample}/{sample}.insert_size_metrics.txt",
            f"{sample_project}/qc/picard/{sample}/{sample}.insert_size_histogram.pdf",
            f"{sample_project}/qc/picard/{sample}/{sample}.gc_bias_metrics.txt",
            f"{sample_project}/qc/picard/{sample}/{sample}.gc_bias_summary_metrics.txt",
            f"{sample_project}/qc/picard/{sample}/{sample}.gc_bias_metrics.pdf",
            f"{sample_project}/qc/rnaseqc2/{sample}/",
        ]
    return qc_files


def get_final_output():
    final_output = [
        "resources/star_genome",
    ]
    for sample in samples.index:
        sample_project = samples.loc[sample, "project_id"]

        final_output += [
            # f"{sample_project}/modification/modtect/{sample}/{sample}.modtect.combined.txt",
        ]
    for dataset in datasets:
        for tool in QUANTIFICATION_TOOLS:

            final_output += [
                f"{sample_project}/DEG/deseq2/{tool}/{dataset}_deg.tsv",
                f"{sample_project}/DEG/edger/{tool}/{dataset}_deg.tsv",
                f"{sample_project}/DEG/limma_trend/{tool}/{dataset}_deg.tsv",
                f"{sample_project}/DEG/limma_voom/{tool}/{dataset}_deg.tsv",
                f"{sample_project}/visualization/{tool}_{dataset}_pca.png",
                f"{sample_project}/visualization/{tool}_{dataset}_pca_batchcorrected.png",
                f"{sample_project}/DEG/{tool}_{dataset}_summary_thresholds.tsv",
                f"{sample_project}/visualization/{tool}_{dataset}_upset.pdf",
                f"{sample_project}/quantification/{tool}/{dataset}_TPM_matrix_corrected.txt",
                f"{sample_project}/mime/{tool}/{dataset}_response_dataset.rds",
            ]

        final_output += [
            f"{sample_project}/purity/PUREE_{dataset}_tumor_purities.tsv",
            f"{sample_project}/enrichment/clusterprofiler/{dataset}_go_enrichment.tsv",
            # f"{sample_project}/visualization/STAR_FC_{dataset}_pca.png",
        ]
    final_output += [
        # f"{sample_project}/visualization/Volcano_validation.pdf",
        # f"{sample_project}/visualization/pca.png",
        # f"{sample_project}/visualization/kallisto_pca.png",
        # f"{sample_project}/visualization/salmon_pca.png",
        f"{sample_project}/qc/multiqc/",
        f"{sample_project}/qc/qc_summary.tsv",
        f"{sample_project}/qc/qc_summary_figures/",
        # f"{sample_project}/enrichment/clusterprofiler/validation_gsea_enrichment.tsv",
        # f"{sample_project}/modification/modtect/merged.modtect.txt",
        # f"{sample_project}/modification/modtect/annotated_modifications.tsv",
        # f"{sample_project}/modification/modtect/annotation_summary.txt",
        # f"{sample_project}/modification/modtect/modtect_statistical_results.tsv",
        # f"{sample_project}/modification/modtect/modtect_significant_modifications.tsv",
        # f"{sample_project}/modification/modtect/modtect_analysis_summary.txt",
        # f"{sample_project}/transcript_splicing/rmats",
        # f"{sample_project}/transcript_splicing/splicetools/",
        # f"{sample_project}/transcript_splicing/rmats_analysis/summary_statistics.csv",
        # f"{sample_project}/transcript_splicing/splicetools_analysis/combined_summary.csv",
        # f"{sample_project}/quantification/STAR_FC/TPM_matrix.txt",
        # f"{sample_project}/dcb/discovery_dcb.tsv",
        # f"{sample_project}/dcb/validation_dcb.tsv",
        # f"{sample_project}/DEG/visualization/deg_heatmap.png",
        # f"{sample_project}/DEG/visualization/discovery_vs_validation_scatter.png",
        # f"{sample_project}/DEG/visualization/significant_genes_boxplot.png",
        # f"{sample_project}/DEG/classifier/lasso_signature_genes.tsv",
        # f"{sample_project}/DEG/classifier/lasso_coefficients.tsv",
        # f"{sample_project}/DEG/classifier/lasso_roc_curve.png",
    ]

    return final_output + get_qc_files()
