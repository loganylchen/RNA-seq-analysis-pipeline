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