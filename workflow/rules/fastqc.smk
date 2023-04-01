rule fastp_fastqc:
    input:
        unpack(get_fq)
    output:
        fq1="results/clean_fastq/{sample}/{sample}_1.fastq.gz",
        fq2="results/clean_fastq/{sample}/{sample}_2.fastq.gz",
        qc_html="results/qc/{sample}/{sample}.fastp.html",
        qc_json="results/qc/{sample}/{sample}.fastp.json",
    log:
        "logs/clean_data/{sample}.log"
    benchmark:
        "benchmarks/{sample}.fastp_qc.benchmark.txt"
    params:
        extra=config['params']["fastp"]
    conda:
        "../envs/fastp.yaml"
    threads:
        config['threads']["fastp"]
    shell:
        "(fastp {params.extra} "
        "-i {input.fq1} "
        "-I {input.fq2} "
        "-o {output.fq1} "
        "-O {output.fq2} "
        "--html {output.qc_html} "
        "--json {output.qc_json}) 2>{log} "
