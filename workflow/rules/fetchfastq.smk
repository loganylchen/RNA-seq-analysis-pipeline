rule get_fastq:
    output:
        unpack(get_raw_fq)
    log:
        "logs/{project}/get_{accession}.fastq.log"
    contaniner:
        'docker://btrspg/sratools:3.2.1'
    params:
        extra="--skip-technical"
    resources:
        mem_mb=10* 1024,
    threads: config['threads']['sra']
    script:
        "../scripts/fetchsra.bash"



