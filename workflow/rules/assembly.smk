rule assembly_stringtie2:
    input:
        bam="{project}/alignment/STAR/{sample}/{sample}.bam",
        rnaseq_qc="{project}/qc/qualimap-rnaseq/{sample}/rnaseq_qc_results.txt",
        gtf="resources/genome.gtf",
    output:
        gtf="{project}/assembly/stringtie/{sample}/{sample}.stringtie.gtf",
    params:
        extra=config["stringtie"]["extra"],
        strand_param=lambda wildcards, input: stringtie_strand_infer(input.rnaseq_qc),
    container:
        (
            "docker://btrspg/stringtie:2.2.3"
            if config["container"].get("stringtie", None) is None
            else config["container"].get("stringtie", None)
        )
    threads: config["threads"].get("stringtie", 4)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("stringtie", 10240),
        tmpdir="./temp",
    log:
        "logs/{project}/{sample}_stringtie.log",
    shell:
        "stringtie "
        "-G {input.gtf} "
        "-p {threads} "
        "{params.strand_param} {params.extra} "
        "-o {output.gtf} "
        "{input.bam} "
        "&>{log}"


rule assembly_merge:
    input:
        gtfs=expand(
            "{project}/assembly/stringtie/{sample}/{sample}.stringtie.gtf",
            project=project,
            sample=discovery_samples.index.tolist(),
        ),
        ref_gtf="resources/genome.gtf",
    output:
        gtf="{project}/assembly/stringtie/merged.gtf",
    params:
        extra=config["stringtie"]["merge_extra"],
    container:
        (
            "docker://btrspg/stringtie:2.2.3"
            if config["container"].get("stringtie", None) is None
            else config["container"].get("stringtie", None)
        )
    threads: config["threads"].get("stringtie", 1)
    resources:
        mem_mb=config["resources"]["mem_mb"].get("stringtie", 10240),
        tmpdir="./temp",
    log:
        "logs/{project}/stringtie_merge.log",
    shell:
        "stringtie "
        "--merge "
        "-G {input.ref_gtf} "
        "{params.extra} "
        "-o {output.gtf} "
        "{input.gtfs} "
        "&>{log}"
