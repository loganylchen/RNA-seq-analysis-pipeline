rule add_read_group:
    input:
        aln="{project}/alignment/STAR/{sample}/{sample}.bam",
    output:
        withrg_bam=temp("{project}/alignment/STAR/{sample}/{sample}.withrg.bam"),
        withrg_bai=temp("{project}/alignment/STAR/{sample}/{sample}.withrg.bai"),
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
        "-I {input.aln} "
        "-O {output.withrg_bam} "
        "--RGPL illumina "
        "--RGLB {wildcards.project} "
        "--RGPU NONE "
        "--RGSM {wildcards.sample} "
        "--CREATE_INDEX true "
        "&>{log};"


rule rnaseq_bam_split:
    input:
        aln="{project}/alignment/STAR/{sample}/{sample}.withrg.bam",
        bai="{project}/alignment/STAR/{sample}/{sample}.withrg.bai",
        fasta="resources/genome.fasta",
        ref_dict="resources/genome.dict",
    output:
        split_bam=temp("{project}/alignment/STAR/{sample}/{sample}.split.bam"),
    resources:
        tmpdir="tmp_dir/",
    log:
        "logs/{project}/{sample}_bam_split.log",
    container:
        (
            "docker://broadinstitute/gatk:4.6.2.0"
            if config["container"].get("gatk4", None) is None
            else config["container"].get("gatk4", None)
        )
    threads: 1
    shell:
        "gatk --java-options '-Djava.io.tmpdir={resources.tmpdir}' "
        "SplitNCigarReads "
        "-R {input.fasta} "
        "-I {input.aln} "
        "-O {output.split_bam} &>{log}"
