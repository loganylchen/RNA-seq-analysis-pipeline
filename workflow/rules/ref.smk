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
        mem_mb=config["resources"]["mem_mb"].get("lftp", 10240),
    benchmark:
        "benchmarks/get_genome.benchmark.txt"
    threads: config["threads"].get("lftp", 4)
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
    threads: config["threads"].get("lftp", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("annotation", 2048),
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
        mem_mb=config["resources"]["mem_mb"].get("filtering", 10240),
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
        mem_mb=config["resources"]["mem_mb"].get("faidx", 10240),
    wrapper:
        "v1.21.4/bio/samtools/faidx"


rule star_index:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    threads: config["threads"].get("star", 8)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("star_index", 102400),
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
        mem_mb=config["resources"]["mem_mb"].get("hisat2_index", 40960),
    log:
        "logs/hisat2_index_genome.log",
    container:
        (
            "docker://btrspg/hisat2:2.2.1"
            if config["container"].get("hisat2", None) is None
            else config["container"].get("hisat2", None)
        )
    threads: config["threads"].get("hisat2", 8)
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
        gtf="{project}/assembly/stringtie/gffcompare.annotated.gtf",
    output:
        fasta="{project}/assembly/stringtie/transcriptome.fasta",
    params:
        extra=config.get("gffread", {}).get("extra", ""),
    threads: config["threads"].get("gffread", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
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
        fasta="{project}/assembly/stringtie/transcriptome.fasta",
    output:
        index=directory("{project}/assembly/stringtie/transcriptome_salmon_index"),
    params:
        extra=config.get("salmon", {}).get("extra_index", ""),
    threads: config["threads"].get("salmon", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("salmon", 10240),
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


rule kallisto_index:
    input:
        fasta="{project}/assembly/stringtie/transcriptome.fasta",
    output:
        index="{project}/assembly/stringtie/transcriptome_kallisto.idx",
    params:
        extra=config.get("kallisto", {}).get("extra_index", ""),
    threads: config["threads"].get("kallisto", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("kallisto", 10240),
    container:
        (
            "docker://btrspg/kallisto:0.51.1"
            if config["container"].get("kallisto", None) is None
            else config["container"].get("kallisto", None)
        )
    log:
        "logs/{project}/kallisto_index.log",
    shell:
        "kallisto index "
        "-i {output.index} "
        "{params.extra} "
        "-t {threads} "
        "{input.fasta} "
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
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
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
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    script:
        "../scripts/gtf2bed.sh"


rule gtf_to_refflat:
    input:
        gtf="resources/genome.gtf",
    output:
        refflat="resources/genome.refFlat",
    log:
        "logs/gtf2refflat.log",
    container:
        (
            "docker://btrspg/python3:20251024"
            if config["container"].get("python3", None) is None
            else config["container"].get("python3", None)
        )
    resources:
        mem_mb=config["resources"]["mem_mb"].get("default", 4096),
    script:
        "../scripts/gtf2refflat.py"


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
    threads: config["threads"].get("picard", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("picard", 8192),
    shell:
        "picard CreateSequenceDictionary "
        "-R {input.fasta} "
        "-O {output.ref_dict} &>{log}"
