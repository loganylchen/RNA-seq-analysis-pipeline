rule star_align:
    input:
        get_clean_data_star,
    output:
        aln="{project}/alignment/{sample}/{sample}.star.bam",
        sj="{project}/alignment/{sample}/{sample}.SJ.out.tab",
    log:
        "logs/{project}_{sample}_star.log",
    params:
        extra=config['star']['extra'],
    threads: config['threads']['star']
    wrapper:
        "v6.0.0/bio/star/align"



rule hisat2_align:
    input:
        get_clean_data_hisat2,
    output:
        "{project}/alignment/{sample}/{sample}.hisat2.bam",
    log:
        "logs/{project}_{sample}_hisat2.log",
    params:
        extra=config['hisat2']['extra'],
    threads: config['threads']['hisat2']
    wrapper:
        "v6.0.0/bio/hisat2/align"