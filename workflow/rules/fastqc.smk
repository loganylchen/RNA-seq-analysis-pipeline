rule fastp_se:
    input:
        sample=["{project}/data/{accession}.fastq.gz"]
    output:
        trimmed="{project}/clean_data/{accession}.fastq.gz",
        html="{project}/report/{accession}.fastp.html",
        json="{project}/report/{accession}.fastp.json",
    log:
        "logs/{project}/fastp_{accession}.log"
    params:
        extra=config['fastp']['se_extra']
    threads: config['threads']['fastp']
    wrapper:
        "v6.0.0/bio/fastp"

rule fastp_pe:
    input:
        sample=["{project}/data/{accession}_1.fastq.gz", "{project}/data/{accession}_2.fastq.gz"]
    output:
        trimmed=["{project}/clean_data/{accession}_1.fastq.gz", "{project}/clean_data/{accession}_1.fastq.gz",],
        html="{project}/report/{accession}.fastp.html",
        json="{project}/report/{accession}.fastp.json"
    log:
        "logs/{project}/fastp_{accession}.log"
    params:
        extra=config['fastp']['pe_extra']
    threads: config['threads']['fastp']
    wrapper:
        "v6.0.0/bio/fastp"


