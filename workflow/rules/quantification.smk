rule featurecounts_quantification:
    input:
        bam = "{project}/alignment/{sample}/{sample}.star.bam",
        gtf = "resources/genome.gtf",
    output:
        counts = "{project}/quantification/{sample}.counts.txt",
        summary = "{project}/quantification/{sample}.counts.summary.txt"
    params:
        # Optional parameters can be added here
        extra = "-p -B -C"  # Example: count paired-end reads, both ends mapped, exclude chimeric reads
    threads: 8
    log:
        "logs/quantification/{sample}.log"
    shell:
        """
        featureCounts -T {threads} {params.extra} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam} > {log} 2>&1
        """