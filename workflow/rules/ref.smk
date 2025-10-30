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


rule get_transcript_sequence:
    input:
        fasta="resources/genome.fasta",
        gtf="{project}/assembly/gffcompare.annotated.gtf",
    output:
        fasta="{project}/assembly/transcriptome.fasta",
    params:
        extra=config.get("gffread", {}).get("extra", ""),
    threads: 1
    resources:
        mem_mb=1024 * 4,
    container:
        (
            "docker://btrspg/gffread:0.12.7"
            if config["container"].get("gffread", None) is None
            else config["container"].get("gffread", None)
        )
    log:
        "logs/{project}/get_transcript_sequence.log",
    shell:
        "gffread -w {output} -g {input.fasta} {input.gtf} {params.extra} &> {log} "


rule salmon_index:
    input:
        fasta="{project}/assembly/transcriptome.fasta",
    output:
        index=directory("{project}/assembly/transcriptome_salmon_index"),
    params:
        extra=config.get("salmon", {}).get("extra_index", ""),
    threads: config["threads"]["salmon"]
    resources:
        mem_mb=1024 * 10,
    container:
        (
            "docker://btrspg/salmon:1.10.3"
            if config["container"].get("salmon", None) is None
            else config["container"].get("salmon", None)
        )
    log:
        "logs/{project}/salmon_index.log",
    shell:
        "salmon index "
        "{params.extra} "
        "-p {threads} "
        "-t {input.fasta} "
        "-i {output.index} "
        "&> {log} "


rule geneid_to_genename:
    input:
        gtf="resources/genome.gtf",
    output:
        tsv="resources/gene_id_to_gene_name.tsv",
    log:
        "logs/geneid2genename.log",
    benchmark:
        "benchmarks/geneid_to_genename.benchmark.txt"
    script:
        "../scripts/get_genename.py"


rule gtf_to_bed:
    input:
        gtf="resources/genome.gtf",
    output:
        bed="resources/transcriptome.bed",
    log:
        "logs/gtf2bed.log",
    benchmark:
        "benchmarks/gtf2bed.benchmark.txt"
    container:
        (
            "docker://btrspg/bedtools:2.31.1"
            if config["container"].get("bedtools", None) is None
            else config["container"].get("bedtools", None)
        )
    script:
        "../scripts/gtf2bed.sh"


rule ref_dict:
    input:
        fasta="resources/genome.fasta",
    output:
        ref_dict="resources/genome.dict",
    log:
        "logs/ref_dict.log",
    container:
        (
            "docker://btrspg/picard:3.4.0"
            if config["container"].get("picard", None) is None
            else config["container"].get("picard", None)
        )
    threads: 1
    shell:
        "picard CreateSequenceDictionary "
        "R={input.fasta} "
        "O={output.ref_dict} &>{log}"
