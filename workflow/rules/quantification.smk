rule featurecounts_quantification_star:
    input:
        bam = "{project}/alignment/{sample}/{sample}.star.bam",
        gtf = "resources/genome.gtf",
    output:
        counts = "{project}/quantification/{sample}.star_counts.txt",
    params:
        extra = config['featurecounts']['extra'],
    conda:
        "../envs/featurecounts.yaml",
    threads: config['threads']['featurecounts'],
    log:
        "logs/quantification/{project}_{sample}_star.log"
    shell:
        "featureCounts -T {threads} {params.extra} "
        "-a {input.gtf} "
        "-o {output.counts} "
        "{input.bam} > {log} 2>&1 "

rule featurecounts_quantification_hisat2:
    input:
        bam = "{project}/alignment/{sample}/{sample}.hisat2.bam",
        gtf = "resources/genome.gtf",
    output:
        counts = "{project}/quantification/{sample}.hisat2_counts.txt",
    params:
        extra = config['featurecounts']['extra'],
    conda:
        "../envs/featurecounts.yaml",
    threads: config['threads']['featurecounts'],
    log:
        "logs/quantification/{project}_{sample}_hisat2.log"
    shell:
        "featureCounts -T {threads} {params.extra} "
        "-a {input.gtf} "
        "-o {output.counts} "
        "{input.bam} > {log} 2>&1 "

