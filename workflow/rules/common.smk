import pandas as pd
import glob
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index())



validate(samples, schema="../schemas/samples.schema.yaml")



def is_pe(wildcards):
    if samples.loc[wildcards.sample].loc['seq_type']  == 'pe':
        return True
    else:
        return False

def get_sra(wildcards):
    return samples.loc[wildcards.sample].loc['raw_data']


def get_raw_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['raw_data']
    seq_type = samples.loc[wildcards.sample].loc['seq_type']
    if seq_type == 'se':
        return {
            'fq1':f"{wildcards.project}/data/{wildcards.sample}/{wildcards.sample}.fastq.gz",
        }
    elif seq_type == 'pe':
        return {
            'fq1':f"{wildcards.project}/data/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
            'fq2':f"{wildcards.project}/data/{wildcards.sample}/{wildcards.sample}_2.fastq.gz"
        }


def get_clean_data_star(wildcards):
    if samples.loc[wildcards.sample].loc['seq_type'] == 'pe':
        return {
            'fq1': f"{wildcards.project}/clean_data/{wildcards.sample}_1.fastq.gz",
            'fq2': f"{wildcards.project}/clean_data/{wildcards.sample}_2.fastq.gz",
            'index':"resources/star_genome",
        }
    elif samples.loc[wildcards.sample].loc['seq_type'] == 'se':
        return {
            'fq1': f"{wildcards.project}/clean_data/{wildcards.sample}.fastq.gz",
            'index': "resources/star_genome",
        }
    else:
        raise ValueError(f'{wildcards.sample} is a wired name!')

def get_clean_data_hisat2(wildcards):
    if samples.loc[wildcards.sample].loc['seq_type'] == 'pe':
        return {
            'reads': [f"{wildcards.project}/clean_data/{wildcards.sample}_1.fastq.gz",
             f"{wildcards.project}/clean_data/{wildcards.sample}_2.fastq.gz"],
             'idx':multiext(
            "resources/hisat2_genome/genome",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
        }
    elif samples.loc[wildcards.sample].loc['seq_type'] == 'se':
        return {
            'reads': [f"{wildcards.project}/clean_data/{wildcards.sample}.fastq.gz"],
             'idx':multiext(
            "resources/hisat2_genome/genome",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
        }
    else:
        raise ValueError(f'{wildcards.sample} is a wired name!')

def get_final_output():
    final_output = []
    for sample in samples.index:
        project = samples.loc[sample,'project']
        final_output+=[
            f"{project}/quantification/{sample}.hisat2_counts.txt",,
            f"{project}/quantification/{sample}.star_counts.txt",,
        ]
    return final_output




