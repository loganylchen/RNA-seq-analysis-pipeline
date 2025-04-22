rule featurecounts_quantification_star:
    input:
        bam = "{project}/alignment/{sample}/{sample}.star.bam",
        gtf = "resources/genome.gtf",
    output:
        counts = "{project}/quantification/{sample}.star_counts.txt",
    params:
        extra = config['featurecounts']['extra']
    threads: config['threads']['featurecounts']
    conda:
        '../envs/featurecounts.yaml'
    log:
        "logs/{project}_{sample}_featurecounts_star.log"
    shell:
        "featureCounts -T {threads} "
        "{params.extra} "
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
        extra = config['featurecounts']['extra']
    threads: config['threads']['featurecounts']
    conda:
        '../envs/featurecounts.yaml'
    log:
        "logs/{project}_{sample}_featurecounts_hisat2.log"
    shell:
        "featureCounts -T {threads} "
        "{params.extra} "
        "-a {input.gtf} "
        "-o {output.counts} "
        "{input.bam} > {log} 2>&1 "
        