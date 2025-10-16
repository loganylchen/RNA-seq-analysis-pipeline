rule star_align:
    input:
        unpack(get_clean_data),
        idx="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln="{project}/alignment/{sample}/{sample}.star.bam",
        reads_per_gene="{project}/quantification/{sample}/{sample}.ReadsPerGene.out.tab",
    log:
        "logs/{project}_{sample}_star.log",
    params:
        extra=lambda wc, input: f'--sjdbGTFfile {input.gtf} {config["star"]["extra"]}', 
    threads: config['threads']['star']
    container:
        (
            "docker://btrspg/star:2.7.11b"
            if config["container"].get("star", None) is None
            else config["container"].get("star", None)
        )
    shell:
        



rule hisat2_align:
    input:
        unpack(get_clean_data),
        idx=multiext(
            "resources/hisat2_genome/genome",
            ".1.ht2",
            ".2.ht2",
            ".3.ht2",
            ".4.ht2",
            ".5.ht2",
            ".6.ht2",
            ".7.ht2",
            ".8.ht2",
        ),
    output:
        "{project}/alignment/{sample}/{sample}.hisat2.bam",
    log:
        "logs/{project}_{sample}_hisat2.log",
    params:
        extra=config['hisat2']['extra'],
    threads: config['threads']['hisat2']
    wrapper:
        "v6.0.0/bio/hisat2/align"