rule sratools_fetchfastq:
    output:
        fq1="results/raw_fastq/{sample}/{sra}_1.fastq.gz",
        fq2="results/raw_fastq/{sample}/{sra}_2.fastq.gz",
        outdir=directory("results/raw_fastq/{sample}")
    log:
        "logs/raw_fastq/{sample}_{sra}_fetchfastq.log"
    params:
        extra=config['params']['sratools_fetchfastq'],
        sample_id="{sra}"
    benchmark:
        "benchmarks/{sample}_{sra}.sratools_fetchfastq.benchmark.txt"
    conda:
        "../envs/sratools.yaml"
    shell:
        "fastq-dump "
        "{params.extra} "
        "--outdir {output.outdir} {params.sample_id} 2>{log}"
