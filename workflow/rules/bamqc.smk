rule qualimap_rnaseq_qc:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
        gtf="resources/genome.gtf",
    output:
        rnaseq_qc="{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
        temp_dir=temp(directory("{project}/qc/temp/{sample}")),
    params:
        pe="-pe" if is_pe else "",
        outdir=lambda wc, output: os.path.dirname(output.rnaseq_qc),
    container:
        (
            "docker://btrspg/qualimap:2.3"
            if config["container"].get("qualimap", None) is None
            else config["container"].get("qualimap", None)
        )
    threads: config["threads"].get("qualimap", 10)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("qualimap", 20480),
    log:
        "logs/{project}/{sample}_qualimap-rnaseq-qc.log",
    shell:
        "mkdir -p {output.temp_dir};"
        "JAVA_OPTS='-Djava.io.tmpdir={output.temp_dir}' qualimap rnaseq "
        "-bam {input.bam} "
        "-gtf {input.gtf} "
        "-outdir {params.outdir} "
        "{params.pe} "
        "--java-mem-size=20G &>{log}"


rule picard_alignment_metrics:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
        ref="resources/genome.fasta",
    output:
        metrics="{project}/qc/picard/{sample}/{sample}.alignment_summary_metrics.txt",
    container:
        (
            "docker://btrspg/picard:3.4.0"
            if config["container"].get("picard", None) is None
            else config["container"].get("picard", None)
        )
    threads: config["threads"].get("picard", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("picard", 8192),
    log:
        "logs/{project}/{sample}_picard_alignment_metrics.log",
    shell:
        "picard CollectAlignmentSummaryMetrics "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.metrics} "
        "&>{log}"


rule picard_rnaseq_metrics:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
        ref="resources/genome.fasta",
        refflat="resources/genome.refFlat",
    output:
        metrics="{project}/qc/picard/{sample}/{sample}.rnaseq_metrics.txt",
    container:
        (
            "docker://btrspg/picard:3.4.0"
            if config["container"].get("picard", None) is None
            else config["container"].get("picard", None)
        )
    threads: config["threads"].get("picard", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("picard", 8192),
    log:
        "logs/{project}/{sample}_picard_rnaseq_metrics.log",
    shell:
        "picard CollectRnaSeqMetrics "
        "-I {input.bam} "
        "-O {output.metrics} "
        "-REF_FLAT {input.refflat} "
        "-STRAND NONE "
        "&>{log}"


rule picard_insert_size_metrics:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
    output:
        metrics="{project}/qc/picard/{sample}/{sample}.insert_size_metrics.txt",
        plot="{project}/qc/picard/{sample}/{sample}.insert_size_histogram.pdf",
    container:
        (
            "docker://btrspg/picard:3.4.0"
            if config["container"].get("picard", None) is None
            else config["container"].get("picard", None)
        )
    threads: config["threads"].get("picard", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("picard", 8192),
    log:
        "logs/{project}/{sample}_picard_insert_size_metrics.log",
    shell:
        "picard CollectInsertSizeMetrics "
        "-I {input.bam} "
        "-O {output.metrics} "
        "-H {output.plot} "
        "&>{log}"


rule picard_gc_bias_metrics:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
        ref="resources/genome.fasta",
    output:
        metrics="{project}/qc/picard/{sample}/{sample}.gc_bias_metrics.txt",
        chart="{project}/qc/picard/{sample}/{sample}.gc_bias_metrics.pdf",
        summary="{project}/qc/picard/{sample}/{sample}.gc_bias_summary_metrics.txt",
    container:
        (
            "docker://btrspg/picard:3.4.0"
            if config["container"].get("picard", None) is None
            else config["container"].get("picard", None)
        )
    threads: config["threads"].get("picard", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("picard", 8192),
    log:
        "logs/{project}/{sample}_picard_gc_bias_metrics.log",
    shell:
        "picard CollectGcBiasMetrics "
        "-I {input.bam} "
        "-O {output.metrics} "
        "-CHART {output.chart} "
        "-S {output.summary} "
        "-R {input.ref} "
        "&>{log}"


rule picard_duplicate_metrics:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
    output:
        metrics="{project}/qc/picard/{sample}/{sample}.duplicate_metrics.txt",
        marked_bam=temp("{project}/qc/picard/{sample}/{sample}.marked_duplicates.bam"),
    container:
        (
            "docker://btrspg/picard:3.4.0"
            if config["container"].get("picard", None) is None
            else config["container"].get("picard", None)
        )
    threads: config["threads"].get("picard", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("picard", 8192),
    log:
        "logs/{project}/{sample}_picard_duplicate_metrics.log",
    shell:
        "picard MarkDuplicates "
        "-I {input.bam} "
        "-O {output.marked_bam} "
        "-M {output.metrics} "
        "-REMOVE_DUPLICATES false "
        "&>{log}"
