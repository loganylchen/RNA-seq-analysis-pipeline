rule fastp_se:
    input:
        sample=["{project}/data/{sample}.fastq.gz"]
    output:
        unpack(get_trimmed_data),
        html="{project}/report/{sample}.fastp.html",
        json="{project}/report/{sample}.fastp.json",
    log:
        "logs/{project}/fastp_{sample}.log"
    params:
        extra=config['fastp']['se_extra']
    threads: config['threads']['fastp']
    wrapper:
        "v6.0.0/bio/fastp"

rule fastp_pe:
    input:
        sample=["{project}/data/{sample}_1.fastq.gz", "{project}/data/{sample}_2.fastq.gz"]
    output:
        unpack(get_trimmed_data),
        html="{project}/report/{sample}.fastp.html",
        json="{project}/report/{sample}.fastp.json"
    log:
        "logs/{project}/fastp_{sample}.log"
    params:
        extra=config['fastp']['pe_extra']
    threads: config['threads']['fastp']
    wrapper:
        "v6.0.0/bio/fastp"


