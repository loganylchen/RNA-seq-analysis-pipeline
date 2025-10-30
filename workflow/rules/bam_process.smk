rule rnaseq_bam_split:
    input:
        aln="{project}/alignment/{sample}/{sample}.star.bam",
        fasta="resources/genome.fasta",
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
