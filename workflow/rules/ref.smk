rule get_genome:
    output:
        "resources/raw_genome.fasta",
    log:
        "logs/ref/get-genome.log",
    container:
        "docker://btrspg/lftp:latest"
    params:
        species=config["reference"]["species"],
        datatype="dna",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
    resources:
        mem_mb=1024 * 10,
    benchmark:
        "benchmarks/get_genome.benchmark.txt"
    threads: config["threads"]["lftp"]
    script:
        "../scripts/get_ensembl_sequence.sh"


rule get_annotation:
    output:
        "resources/raw_genome.gtf",
    params:
        species=config["reference"]["species"],
        fmt="gtf",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        flavor="",
    container:
        "docker://btrspg/lftp:latest"
    threads: 1
    resources:
        mem_mb=1024 * 2,
    log:
        "logs/ref/get_annotation.log",
    benchmark:
        "benchmarks/get_annotation.benchmark.txt"
    script:
        "../scripts/get_ensembl_annotation.sh"


rule filtering_genome_and_annotation:
    input:
        fasta="resources/raw_genome.fasta",
        gtf="resources/raw_genome.gtf",
    output:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    log:
        "logs/ref/filtering_references.log",
    resources:
        mem_mb=1024 * 10,
    params:
        select_contigs=config["reference"]["select_contigs"],
    container:
        "docker://btrspg/biopython:1.85"
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
    resources:
        mem_mb=1024 * 10,
    wrapper:
        "v1.21.4/bio/samtools/faidx"


rule star_index:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: config["threads"]["star"]
    resources:
        mem_mb=1024 * 100,
    params:
        extra=config["star"]["index_extra"],
    log:
        "logs/star_index_genome.log",
    benchmark:
        "benchmarks/star_index.benchmark.txt"
    container:
        (
            "docker://btrspg/star:2.7.11b"
            if config["container"].get("star", None) is None
            else config["container"].get("star", None)
        )
    shell:
        "STAR --runMode genomeGenerate "
        "--genomeDir {output} "
        "--genomeFastaFiles {input.fasta} "
        "--sjdbGTFfile {input.gtf} "
        "{params.extra} "
        "--runThreadN {threads} &>{log} "


rule hisat2_index:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
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
        prefix="resources/hisat2_genome/genome",
        extra=config["hisat2"]["index_extra"],
    resources:
        mem_mb=1024 * 40,
    log:
        "logs/hisat2_index_genome.log",
    container:
        (
            "docker://btrspg/hisat2:2.2.1"
            if config["container"].get("hisat2", None) is None
            else config["container"].get("hisat2", None)
        )
    threads: config["threads"]["hisat2"]
    benchmark:
        "benchmarks/hisat2_index.benchmark.txt"
    shell:
        "hisat2_extract_splice_sites.py {input.gtf} > {params.prefix}.ss 2>{log};"
        "hisat2_extract_exons.py {input.gtf} > {params.prefix}.exon 2>>{log}; "
        "hisat2-build --threads {threads} "
        "--exon {params.prefix}.exon "
        "--ss {params.prefix}.ss {input.fasta} {params.prefix} &>>{log} "
