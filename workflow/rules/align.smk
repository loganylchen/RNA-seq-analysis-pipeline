rule star_align:
    input:
        unpack(get_clean_data),
        idx="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln="{project}/alignment/STAR/{sample}/{sample}.bam",
        temp_dir=directory("{project}/alignment/STAR_TMP/{sample}"),
        reads_per_gene="{project}/quantification/STAR/{sample}/{sample}.ReadsPerGene.out.tab",
        qc_log="{project}/qc/STAR/{sample}/{sample}.Log.final.out",
    log:
        "logs/{project}/{sample}_star.log",
    container:
        (
            "docker://btrspg/star:2.7.11b"
            if config["container"].get("star", None) is None
            else config["container"].get("star", None)
        )
    params:
        extra=lambda wc, input: f' --quantMode GeneCounts --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {input.gtf} {config["star"]["extra"]}',
    threads: config["threads"]["star"]
    shell:
        "STAR --genomeDir {input.idx} "
        "--readFilesIn {input.reads} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {output.temp_dir}/ "
        "--outTmpDir {output.temp_dir}/_STARTMP "
        "--outStd BAM_SortedByCoordinate "
        "--outReadsUnmapped Fastx"
        "{params.extra} "
        "--runThreadN {threads} "
        "> {output.aln} 2>{log};"
        "mv {output.temp_dir}/ReadsPerGene.out.tab {output.reads_per_gene}; "
        "mv {output.temp_dir}/Log.final.out {output.qc_log}; "


rule hisat2_align:
    input:
        unpack(get_clean_data),
        idx="resources/hisat2_genome/genome.1.ht2",
        rnaseq_qc="{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
    output:
        aln="{project}/alignment/hisat2/{sample}/{sample}.bam",
        qc_log="{project}/qc/hisat2/{sample}/{sample}.log",
    log:
        "logs/{project}/{sample}_hisat2.log",
    container:
        (
            "docker://btrspg/hisat2:2.2.1"
            if config["container"].get("hisat2", None) is None
            else config["container"].get("hisat2", None)
        )
    params:
        extra=config["hisat2"]["extra"],
        prefix="resources/hisat2_genome/genome",
        strand_param=lambda wildcards, input: hisat2_strand_infer(input.rnaseq_qc),
    threads: config["threads"].get("hisat2", 10)
    shell:
        "hisat2 -x {params.prefix} "
        "-1 {input.fq1} "
        "-2 {input.fq2} "
        "--new-summary {output.qc_log} "
        "{params.strand_param} "
        "{params.extra} "
        "-p {threads} "
        "-S {output.aln} "
        "&>{log};"
