rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/ref/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    benchmark:
        "benchmarks/get_genome.benchmark.txt"
    wrapper:
        "v1.21.4/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/ref/get_annotation.log",
    benchmark:
        "benchmarks/get_annotation.benchmark.txt"
    wrapper:
        "v1.21.4/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/ref/genome-faidx.log",
    cache: True
    wrapper:
        "v1.21.4/bio/samtools/faidx"



rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: 4
    params:
        extra="--sjdbGTFfile resources/genome.gtf --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    benchmark:
        "benchmarks/star_index.benchmark.txt"
    wrapper:
        "v1.21.4/bio/star/index"

rule get_vep_cache:
    output:
        directory("resources/vep/cache"),
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    log:
        "logs/vep/cache.log",
    wrapper:
        "v1.25.0/bio/vep/cache"

rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=100
    log:
        "logs/download_vep_plugins.log"
    wrapper:
        "v1.25.0/bio/vep/plugins"

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

rule gtf_to_bed:
    input:
        "resources/genome.gtf"
    output:
        "resources/genome.bed"
    log:
        "logs/gtf2bed/gtf2bed.log"
    conda:
        "../envs/pygtftk.yaml"
    shell:
        "gtftk select_by_key -e -i {input} | "
        "gtftk convert --outputfile {output} -f bed6"

rule bedtools_merge:
    input:
        "resources/genome.sorted.bed"
    output:
        "resources/genome.merged.bed"
    params:
        ## Add optional parameters
        extra=" " ## In this example, we want to count how many input lines we merged per output line
    log:
        "logs/merge.bedtools.log"
    wrapper:
        "v2.3.1/bio/bedtools/merge"

rule bedtools_sort_bed:
    input:
        in_file="resources/genome.bed",
        # an optional sort file can be set as genomefile by the variable genome or
        # as fasta index file by the variable faidx
        genome="resources/genome.fasta.fai"
    output:
        "resources/genome.sorted.bed"
    params:
        ## Add optional parameters
        extra=""
    log:
        "logs/sorted.bed.log"
    wrapper:
        "v2.3.1/bio/bedtools/sort"

