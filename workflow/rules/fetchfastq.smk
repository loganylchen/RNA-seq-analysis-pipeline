rule sratools_fetchfastq_pe:
    output:
        fq1="results/raw_fastq/{sample}/{sample}_1.fastq.gz",
        fq2="results/raw_fastq/{sample}/{sample}_2.fastq.gz",
        outdir=directory("results/raw_fastq/{sample}")
    log:
        "logs/raw_fastq/{sample}_fetchfastq.log"
    params:
        extra=config['params']['sratools_fetchfastq'],
        sra=get_sra
    benchmark:
        "benchmarks/{sample}.sratools_fetchfastq.benchmark.txt"
    conda:
        "../envs/sratools.yaml"
    shell:
        "fastq-dump "
        "{params.extra} "
        "--outdir {output.outdir} {params.sra} 2>{log} && "
        "mv {output.outdir}/{params.sra}_1.fastq.gz {output.fq1} && "
        "mv {output.outdir}/{params.sra}_2.fastq.gz {output.fq2} "


rule sratools_fetchfastq_se:
    output:
        fq="results/raw_fastq/{sample}/{sample}.fastq.gz",
        outdir=directory("results/raw_fastq/{sample}")
    log:
        "logs/raw_fastq/{sample}_fetchfastq.log"
    params:
        extra=config['params']['sratools_fetchfastq'],
        sra=get_sra
    benchmark:
        "benchmarks/{sample}.sratools_fetchfastq.benchmark.txt"
    conda:
        "../envs/sratools.yaml"
    shell:
        "fastq-dump "
        "{params.extra} "
        "--outdir {output.outdir} {params.sra} 2>{log} && "
        "mv {output.outdir}/{params.sra}.fastq.gz {output.fq} "
