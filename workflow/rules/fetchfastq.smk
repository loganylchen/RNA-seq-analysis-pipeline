rule get_fastq_pe:
    output:
        "{project}/data/{accession}_1.fastq.gz",
        "{project}/data/{accession}_2.fastq.gz",
    log:
        "logs/{project}/get_{accession}.fastq.log"
    params:
        extra="--skip-technical"
    threads: config['threads']['sra']
    wrapper:
        "v6.0.0/bio/sra-tools/fasterq-dump"

rule get_fastq_se:
    output:
        "{project}/data/{accession}.fastq.gz"
    log:
        "logs/{project}/get_{accession}.fastq.log"
    params:
        extra="--skip-technical"
    threads: config['threads']['sra']
    wrapper:
        "v6.0.0/bio/sra-tools/fasterq-dump"


