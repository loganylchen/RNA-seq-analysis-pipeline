rule star_align:
    input:
        unpack(get_clean_data),
        idx="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln="{project}/alignment/{sample}/{sample}.star.bam",
        temp_dir=temp(directory("{project}/alignment/{sample}/STAR_TMP")),
        reads_per_gene="{project}/quantification/{sample}/{sample}.ReadsPerGene.out.tab",
        qc_log="{project}/qc/{sample}/{sample}.Log.final.out",
    log:
        "logs/{project}/{sample}_star.log",
    container:
        (
            "docker://btrspg/star:2.7.11b"
            if config["container"].get("star", None) is None
            else config["container"].get("star", None)
        )
    params:
        extra=lambda wc, input: f'--quantMode GeneCounts --chimOutJunctionFormat 1 --outSAMtype BAM SortedByCoordinate --sjdbGTFfile {input.gtf} {config["star"]["extra"]}',
    threads: config["threads"]["star"]
    shell:
        "STAR --genomeDir {input.idx} "
        "--readFilesIn {input.reads} "
        "--readFilesCommand zcat "
        "--outFileNamePrefix {output.temp_dir}/ "
        "--outTmpDir {output.temp_dir}/_STARTMP "
        "--outStd BAM_SortedByCoordinate "
        "{params.extra} "
        "--runThreadN {threads} "
        "> {output.aln} 2>{log};"
        "mv {output.temp_dir}/ReadsPerGene.out.tab {output.reads_per_gene}; "
        "mv {output.temp_dir}/Log.final.out {output.qc_log}; "



