rule add_read_group:
    input:
        aln="{project}/alignment/{sample}/{sample}.star.bam",
    output:
        withrg_bam=temp("{project}/alignment/{sample}/{sample}.withrg.bam"),
    log:
        "logs/{project}/{sample}_add_readgroup.log",
    container:
        (
            "docker://btrspg/picard:3.4.0"
            if config["container"].get("picard", None) is None
            else config["container"].get("picard", None)
        )
    threads: 1
    shell:
        "picard AddOrReplaceReadGroups "
        "I={input.aln} "
        "O={output.withrg_bam} "
        "RGPL=illumina RGLB={wildcards.project} RGPU=NONE RGSM={wildcards.sample} "
        "&>{log};"
        "picard BuildBamIndex I={output.withrg_bam};"


rule rnaseq_bam_split:
    input:
        aln="{project}/alignment/{sample}/{sample}.star.bam",
        fasta="resources/genome.fasta",
        ref_dict="resources/genome.dict",
    output:
        split_bam="{project}/alignment/{sample}/{sample}.split.bam",
    log:
        "logs/{project}/{sample}_bam_split.log",
    container:
        (
            "docker://btrspg/gatk:3.8"
            if config["container"].get("gatk", None) is None
            else config["container"].get("gatk", None)
        )
    threads: config["threads"]["gatk"]
    shell:
        "gatk3 -T SplitNCigarReads "
        "-R {input.fasta} "
        "-I {input.aln} "
        "--num_threads {threads} "
        "-o {output.split_bam} &>{log}"
