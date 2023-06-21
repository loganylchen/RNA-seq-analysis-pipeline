import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")


def check_raw_data(raw_data_string:str):
    if raw_data_string.endswith('.fq') or raw_data_string.endswith('.fq.gz') or raw_data_string.endswith('.fastq') or raw_data_string.endswith('.fastq.gz'):
        fq1, fq2 = raw_data_string.split(',')
        return 'fastq',fq1,fq2
    elif raw_data_string.startswith('SRR'):
        return 'sra',raw_data_string, None
    else:
        raise ValueError(f'{raw_data_string} is not a valide datatype')
def get_fq(wildcards):
    raw_data = samples.loc[wildcards.sample].loc['raw_data']
    data_type , *data = check_raw_data(raw_data)
    if data_type == 'fastq':
        return {
            'fq1':data[0],
            'fq2':data[1]
        }
    else:
        return {
            'fq1':f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_1.fastq.gz",
            'fq2':f"results/raw_fastq/{wildcards.sample}/{wildcards.sample}_2.fastq.gz"
        }

def get_sra(wildcards):
    return samples.loc[wildcards.sample].loc['raw_data']

def get_final_output():
    contrasts = config['diffexp']['contrasts']
    final_output = expand("results/star/{sample.sample_name}/ReadsPerGene.out.tab",sample=samples.itertuples())
    final_output.append("results/deseq2/count_matrix.rds")
    for key in contrasts:
        final_output.append(f"results/diffexp/{key}.diffexp.tsv")
        final_output.append(directory(f"results/enrichment/{key}"))
        final_output.append(directory(f"results/visualization/Volcano.{key}.diffexp.pdf"))
        final_output.append(directory(f"results/visualization/Volcano.{key}.diffexp.png"))
    # final_output.append("results/counts/all.symbol.tsv")
    return final_output



def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]
# units = (
#     pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
#     .set_index(["sample_name", "unit_name"], drop=False)
#     .sort_index()
# )
# validate(units, schema="../schemas/units.schema.yaml")
#
#
# def get_cutadapt_input(wildcards):
#     unit = units.loc[wildcards.sample].loc[wildcards.unit]
#
#     if pd.isna(unit["fq1"]):
#         # SRA sample (always paired-end for now)
#         accession = unit["sra"]
#         return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])
#
#     if unit["fq1"].endswith("gz"):
#         ending = ".gz"
#     else:
#         ending = ""
#
#     if pd.isna(unit["fq2"]):
#         # single end local sample
#         return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
#             S=unit.sample_name, U=unit.unit_name, E=ending
#         )
#     else:
#         # paired end local sample
#         return expand(
#             "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
#                 S=unit.sample_name, U=unit.unit_name, E=ending
#             ),
#             read=["fq1", "fq2"],
#         )
#
#
# def get_cutadapt_pipe_input(wildcards):
#     files = list(
#         sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
#     )
#     assert len(files) > 0
#     return files
#
#
# def is_paired_end(sample):
#     sample_units = units.loc[sample]
#     fq2_null = sample_units["fq2"].isnull()
#     sra_null = sample_units["sra"].isnull()
#     paired = ~fq2_null | ~sra_null
#     all_paired = paired.all()
#     all_single = (~paired).all()
#     assert (
#         all_single or all_paired
#     ), "invalid units for sample {}, must be all paired end or all single end".format(
#         sample
#     )
#     return all_paired
#
#

#
#
# def get_strandedness(units):
#     if "strandedness" in units.columns:
#         return units["strandedness"].tolist()
#     else:
#         strand_list = ["none"]
#         return strand_list * units.shape[0]
#
#
# def get_deseq2_threads(wildcards=None):
#     # https://twitter.com/mikelove/status/918770188568363008
#     few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
#     return 1 if len(samples) < 100 or few_coeffs else 6
#
#
# def is_activated(xpath):
#     c = config
#     for entry in xpath.split("/"):
#         c = c.get(entry, {})
#     return bool(c.get("activate", False))
#
#
# def get_bioc_species_name():
#     first_letter = config["ref"]["species"][0]
#     subspecies = config["ref"]["species"].split("_")[1]
#     return first_letter + subspecies
#
#
# def get_fastqs(wc):
#     if config["trimming"]["activate"]:
#         return expand(
#             "results/trimmed/{sample}/{unit}_{read}.fastq.gz",
#             unit=units.loc[wc.sample, "unit_name"],
#             sample=wc.sample,
#             read=wc.read,
#         )
#     unit = units.loc[wc.sample]
#     if all(pd.isna(unit["fq1"])):
#         # SRA sample (always paired-end for now)
#         accession = unit["sra"]
#         return expand(
#             "sra/{accession}_{read}.fastq", accession=accession, read=wc.read[-1]
#         )
#     fq = "fq{}".format(wc.read[-1])
#     return units.loc[wc.sample, fq].tolist()
#
#
# def get_contrast(wildcards):
#     return config["diffexp"]["contrasts"][wildcards.contrast]