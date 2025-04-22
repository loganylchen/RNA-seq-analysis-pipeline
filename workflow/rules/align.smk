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
        get_clean_data_hisat2,
    output:
        "{project}/alignment/{accession}/{accession}.hisat2.bam",
    log:
        "logs/{project}_{accession}_hisat2.log",
    params:
        extra=config['hisat2']['extra'],
    threads: config['threads']['hisat2']
    wrapper:
        "v6.0.0/bio/hisat2/align"