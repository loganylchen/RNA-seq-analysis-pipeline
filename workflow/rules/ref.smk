rule get_genome:
    output:
        "resources/raw_genome.fasta",
    log:
        "logs/ref/get-genome.log",
    params:
        species=config["reference"]["species"],
        datatype="dna",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
    benchmark:
        "benchmarks/get_genome.benchmark.txt"
    wrapper:
        "v1.21.4/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/raw_genome.gtf",
    params:
        species=config["reference"]["species"],
        fmt="gtf",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        flavor="",
    log:
        "logs/ref/get_annotation.log",
    benchmark:
        "benchmarks/get_annotation.benchmark.txt"
    wrapper:
        "v1.21.4/bio/reference/ensembl-annotation"


rule filtering_genome_and_annotation:
    input:
        fasta="resources/raw_genome.fasta",
        gtf="resources/raw_genome.gtf",
    output:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    log:
        "logs/ref/filtering_references.log",
    params:
        select_contigs=config["reference"]["select_contigs"],
    container:
        "docker://btrspg/biopython:1.85"
    conda:
        "../envs/biopython.yaml"
    benchmark:
        "benchmarks/filtering_references.benchmark.txt"
    script:
        "../scripts/reference_filtering.py"

rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/ref/genome-faidx.log",
    wrapper:
        "v1.21.4/bio/samtools/faidx"

rule create_dict:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.dict",
    log:
        "logs/picard/create_dict.log",
    params:
        extra="",  # optional: extra arguments for picard.
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=10240,
    wrapper:
        "v2.3.1/bio/picard/createsequencedictionary"



rule star_index:
    input:
        fasta='resources/genome.fasta',
        gtf='resources/genome.gtf'
    output:
        directory("resources/star_genome")
    threads: config['threads']['star']
    params:
        sjdb_overhang=100,
        extra=config['star']['index_extra'],
    log:
        "logs/star_index_genome.log",
    benchmark:
        "benchmarks/star_index.benchmark.txt"
    wrapper:
        "v1.21.4/bio/star/index"

rule hisat2_index:
    input:
        fasta='resources/genome.fasta'
    output:
        multiext(
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
    params:
        extra=config['hisat2']['index_extra'],
    log:
        "logs/hisat2_index_genome.log",
    threads: config['threads']['hisat2'],
    wrapper:
        "v6.0.0/bio/hisat2/index"


