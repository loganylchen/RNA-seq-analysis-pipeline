rule star_align:
    input:
        get_clean_data_star,
    output:
        aln="{project}/alignment/{accession}/{accession}.star.bam",
        sj="{project}/alignment/{accession}/{accession}.SJ.out.tab",
    log:
        "logs/{project}_{accession}_star.log",
    params:
        extra=config['star']['extra'],
    threads: config['threads']['star']
    wrapper:
        "v6.0.0/bio/star/align"




rule hisat2_align:
    input:
        get_clean_data_hisat2
        
    output:
        "mapped/{sample}.bam",
    log:
        "logs/hisat2_align_{sample}.log",
    params:
        extra="",
    threads: 2
    wrapper:
        "v6.0.0/bio/hisat2/align"