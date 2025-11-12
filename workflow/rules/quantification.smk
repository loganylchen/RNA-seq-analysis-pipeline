rule featurecounts_quantification_star:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
        rnaseq_qc="{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
        gtf="resources/genome.gtf",
    output:
        counts="{project}/quantification/featurecounts/{sample}.txt",
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
        "logs/{project}/{sample}_featurecounts_star.log",
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
        idx="{project}/assembly/stringtie/transcriptome_salmon_index",
        rnaseq_qc="{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
    output:
        outdir=directory("{project}/quantification/salmon/{sample}/"),
        quantification_file="{project}/quantification/salmon/{sample}/quant.sf",
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


rule kallisto_quantification:
    input:
        unpack(get_clean_data),
        idx="{project}/assembly/stringtie/transcriptome_kallisto.idx",
        rnaseq_qc="{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
    output:
        quantification_file="{project}/quantification/kallisto/{sample}/abundance.tsv",
        qc_log="{project}/qc/kallisto/{sample}/kallisto.log",
    params:
        extra=config.get("kallisto", {}).get("extra", ""),
        outdir=lambda wc, output: os.path.dirname(output),
        strand_param=lambda wildcards, input: kallisto_strand_infer(input.rnaseq_qc),
    threads: config["threads"].get("kallisto", 1)
    resources:
        mem_mb=1024 * 10,
    container:
        (
            "docker://btrspg/kallisto:0.51.1"
            if config["container"].get("kallisto", None) is None
            else config["container"].get("kallisto", None)
        )
    log:
        "logs/{project}/{sample}_kallisto_quantify.log",
    shell:
        "kallisto quant "
        "-i {input.idx} "
        "{params.strand_param} "
        "{params.extra} "
        "-t {threads} "
        "-o {params.outdir} "
        "{input.fq1} "
        "{input.fq2} "
        "&> {log}; "
        "cp {log} {output.qc_log};"
