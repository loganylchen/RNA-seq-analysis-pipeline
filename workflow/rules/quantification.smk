rule featurecounts_quantification_star:
    input:
        bam="{project}/alignment/{sample}/{sample}.star.bam",
        rnaseq_qc="{project}/qc/{sample}/rnaseq_qc_results.txt",
        gtf="resources/genome.gtf",
    output:
        counts="{project}/quantification/{sample}.star_counts.txt",
    params:
        extra=config["featurecounts"]["extra"],
        strand_param=lambda wildcards, input: featurecounts_strand_infer(
            input.rnaseq_qc
        ),
    container:
        (
            "docker://btrspg/subread:2.1.1"
            if config["container"].get("subread", None) is None
            else config["container"].get("subread", None)
        )
    threads: config["threads"]["featurecounts"]
    log:
        "logs/quantification/{project}_{sample}_star.log",
    shell:
        "featureCounts -T {threads} "
        "{params.extra} "
        "{params.strand_param} "
        "-a {input.gtf} "
        "-o {output.counts} "
        "{input.bam} &>{log}  "


rule salmon_quantification:
    input:
        unpack(get_clean_data),
        idx="{project}/assembly/transcriptome_salmon_index",
        rnaseq_qc="{project}/qc/{sample}/rnaseq_qc_results.txt",
    output:
        outdir=directory("{project}/quantification/{sample}.salmon"),
    params:
        extra=config.get("salmon", {}).get("extra", ""),
        strand_param=lambda wildcards, input: salmon_strand_infer(input.rnaseq_qc),
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
        "logs/{project}/{sample}_salmon_quantify.log",
    shell:
        "salmon quant "
        "-i {input.idx} "
        "{params.strand_param} "
        "{params.extra} "
        "-1 {input.fq1} "
        "-2 {input.fq2} "
        "-p {threads} "
        "-o {output.outdir} "
        "&> {log} "
