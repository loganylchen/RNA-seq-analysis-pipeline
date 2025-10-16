rule get_fastq:
    output:
        reads=(
            [
                "{project}/data/{sample}/{sample}_1.fastq.gz",
                "{project}/data/{sample}/{sample}_2.fastq.gz",
            ]
            if is_pe
            else ["{project}/data/{sample}/{sample}.fastq.gz"]
        ),
    log:
        "logs/{project}/get_{sample}.fastq.log",
    params:
        sra=get_sra,
        fq1=get_fq1,
        fq2=get_fq2,
    container:
        "docker://btrspg/sratools:3.2.1"
    resources:
        mem_mb=10 * 1024,
    threads: config["threads"]["sra"]
    script:
        "../scripts/fetchsra.sh"
