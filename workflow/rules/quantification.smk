rule featurecounts_quantification_star:
    input:
        bam="{project}/alignment/{sample}/{sample}.star.bam",
        gtf="resources/genome.gtf",
    output:
        counts="{project}/quantification/{sample}.star_counts.txt",
    params:
        extra=config["featurecounts"]["extra"],
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
        "featureCounts -T {threads} {params.extra} "
        "-a {input.gtf} "
        "-o {output.counts} "
        "{input.bam} &>{log}  "
